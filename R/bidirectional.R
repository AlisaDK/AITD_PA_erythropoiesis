library(TwoSampleMR)
library(tidyverse)
library(RadialMR)
library(MendelianRandomization)


#EXPOSURE: PA
pernicious_anemia <- read_delim("data-raw/b12_meta_ukbb_estbb_finngen_061020.out",
                              "\t", escape_double = FALSE, trim_ws = TRUE)

PA <- pernicious_anemia %>%
    filter(`rs_number` %in% c("rs6679677", "rs12616502", "rs28414666", "rs2476491", "rs74203920"))

PA <- data.frame(SNP = PA$`rs_number`, beta.exposure= PA$beta, se.exposure= PA$se, id.exposure=1, exposure="Pernicious anemia",
           units.exposure="log odds", pval.exposure=PA$`p-value`, prevalence.exposure=0.0033, ncontrol.exposure=PA$n_samples,  ncase.exposure=round(0.0033*PA$n_samples),
              effect_allele.exposure= PA$reference_allele, other_allele.exposure= PA$other_allele,  eaf.exposure= PA$eaf)
dat_exposure <- PA

#OUTCOME: AITD
AITD <- read_table2("data-raw/Autoimmune_thyroid_disease_Meta_results_200313.cl.txt")
head(AITD)
AITD <- AITD %>%
    rename(SNP=rsID) %>%
 filter(SNP %in% PA$SNP)
AITD$`UKB-frq` <-as.numeric(AITD$`UKB-frq`)
AITD <- AITD %>%
 mutate(beta=log(`OR-A1`)) %>%
 mutate(z=qnorm(P)) %>%
 mutate(se=abs(beta/z)) %>%
 mutate(eaf=`UKB-frq`/100)
dat_outcome <-data.frame(SNP = AITD$SNP, beta.outcome= AITD$beta, se.outcome=AITD$se, id.outcome=1, outcome="AITD",
                   units.outcome="log odds", pval.outcome=AITD$P, ncase.outcome=30234,
                       ncontrol.outcome=725172, prevalence.outcome=0.038,
                       effect_allele.outcome=AITD$A1, other_allele.outcome=AITD$A0, eaf.outcome=AITD$eaf)
#HARMONIZING
harm_dat <- harmonise_data(dat_exposure, dat_outcome)

snps_reverse <- harm_dat %>%
    select(exposure, outcome, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, ncase.exposure, ncontrol.exposure,
           beta.outcome, se.outcome,  pval.outcome, ncase.outcome, ncontrol.outcome) %>%
    rename(effect_allele=effect_allele.exposure, other_allele=other_allele.exposure) %>%
    mutate(F.stat=round((beta.exposure*beta.exposure)/(se.exposure*se.exposure)))

snps_reverse <-snps_reverse[order(snps_reverse$SNP),]
save(snps_reverse, file = "data/snps_reverse.RData")

#Steiger filtering: error note, so I calculated r2 by get_R_from_lor function below
#harm_dat.steiger<-steiger_filtering(harm_dat) %>%
#    filter(steiger_dir=="TRUE")

binary.exposure <- harm_dat %>%
    filter(units.exposure=="log odds") %>%
    mutate(r.exposure=NA)

binary.exposure$r.exposure <- get_r_from_lor(
    lor=binary.exposure$beta.exposure,
    af=binary.exposure$eaf.exposure,
    ncase=binary.exposure$ncase.exposure,
    ncontrol=binary.exposure$ncontrol.exposure,
    prevalence=binary.exposure$prevalence.exposure,
    model = "logit",
    correction = FALSE
)

binary.outcome <- harm_dat %>%
    filter(units.outcome=="log odds") %>%
    mutate(r.outcome=NA)

binary.outcome$r.outcome <- get_r_from_lor(
    lor=binary.outcome$beta.outcome,
    af=binary.outcome$eaf.outcome,
    ncase=binary.outcome$ncase.outcome,
    ncontrol=binary.outcome$ncontrol.outcome,
    prevalence=binary.outcome$prevalence.outcome,
    model = "logit",
    correction = FALSE
)

harm_dat_new <- left_join(harm_dat, binary.exposure)
harm_dat <- left_join(harm_dat_new, binary.outcome)
names(harm_dat)
reverse_r2 <- harm_dat[,c(1, 23, 15, 32,33)] %>%
    mutate(r2_exposure=r.exposure^2, r2_outcome=r.outcome^2,) %>%
    select(, !contains("r."))
#reverse_r2 holds r2 for each of the 5 individual SNPs showing that r2exp>r2out

#Radial MR: note ivw_radial does not work every time, probably depends on iterations. However, every time, the result is the same
radial_data<-format_radial(harm_dat$beta.exposure,harm_dat$beta.outcome,
                            harm_dat$se.exposure,harm_dat$se.outcome,
                            harm_dat$SNP)

ivw_model<-ivw_radial(radial_data,0.05/nrow(radial_data),3,0.00001)
ivw_model$outliers
radial_outliers<-radial_data[radial_data$SNP %in% ivw_model$outliers$SNP,]
radial_data2<-radial_data[-c(as.numeric(row.names(radial_outliers))),]

harm_radial <- harm_dat %>%
 filter(SNP %in% radial_data2$SNP)
mrradial1 <- mr(harm_radial, method_list=c( "mr_ivw"))
radial_outliers1 <- radial_outliers %>%
 summarise(snp.outliers=n()) %>%
 mutate(outcome="AITD")
radial_outliers2 <- radial_outliers %>%
 summarise(snp.list = paste0(SNP, collapse = " ")) %>%
 mutate(outcome="AITD")

radial_outliers <- left_join(radial_outliers1, radial_outliers2)

#HETREOGENEITY
het<-mr_heterogeneity(harm_radial)
het <- het %>%
 mutate(I2=100*(Q-Q_df)/Q) %>%
 mutate(I2=if_else(I2<0, 0, I2)) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/
het <- subset(het, method=="Inverse variance weighted")
mrradial <- merge(mrradial1, het, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x = TRUE, all.y = T) %>%
 mutate(method="IVW, -RadialMR outliers")
# #Merge Radial MR results
mrradial <- merge(mrradial1, het, by=c("outcome"))
mrradial_all <- merge(mrradial, radial_outliers)
mrradial_all$snp.outliers <-as.numeric(mrradial_all$snp.outliers)
mrradial <- mrradial_all %>%
 mutate(snp.outliers = if_else(is.na(snp.outliers), 0, snp.outliers))

mrradial <- mrradial %>%
 rename(id.exposure=id.exposure.x, id.outcome=id.outcome.x, exposure=exposure.x, method=method.x) %>%
 select(, !contains(".y")) %>%
 mutate(method="IVW, -RadialMR outliers")
save(mrradial, file = "data/mrradial_bidirectional.RData")
#load("data/mrradial_bidirectional.RData")


##############################
#MR analyses
##############################
res <- mr(harm_dat, method_list=c( "mr_ivw", "mr_weighted_median", "mr_egger_regression"))

mrpresso<-run_mr_presso(harm_dat, NbDistribution = 1000)
mr_res <- generate_odds_ratios(res)

scatter_plot <- mr_scatter_plot(mr_res, harm_dat)
ggsave(scatter_plot[[1]], file="scatter_PA_AITD.pdf", width=7, height=7)

single_snp_analysis <- mr_singlesnp(harm_dat)
forest_plot <- mr_forest_plot(single_snp_analysis)
ggsave(forest_plot[[1]], file="forest_plot_PA_AITD.pdf", width=7, height=7)

ivw <- subset_on_method(res,
                        single_snp_method = "Wald ratio",
                        multi_snp_method = "Inverse variance weighted"
)
wm <- subset_on_method(
    res,
    single_snp_method = "Wald ratio",
    multi_snp_method = "Weighted median"
)
egger <- subset_on_method(
    res,
    single_snp_method = "Wald ratio",
    multi_snp_method = "MR Egger"
)

#Steiger: all OK

#HETREOGENEITY
het<-mr_heterogeneity(harm_dat)
het <- het %>%
    mutate(I2=100*(Q-Q_df)/Q) %>%
    mutate(I2=if_else(I2<0, 0, I2)) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/
het <- subset(het, method=="Inverse variance weighted")

#PLEIOTROPY
#Egger intercept
plt<-mr_pleiotropy_test(harm_dat)
plt <- plt %>%
    rename(se.egger=se, pval.egger=pval) %>%
    mutate(method="MR Egger")

#Egger I2GX
i2gxf <- function(x, y) {
    list(id.exposure=y$id.exposure, id.outcome=y$id.outcome, I2GX=mr_egger(mr_input(bx=x$beta.exposure, bxse=x$se.exposure, by=x$beta.outcome, byse=x$se.outcome))@I.sq)
}

I2GX <- harm_dat %>%
    group_by(id.exposure, id.outcome) %>%
    group_map(i2gxf)

df <- data.frame(matrix(unlist(I2GX), nrow=length(I2GX), byrow=T),stringsAsFactors=FALSE)
mregger.i2gx <- df %>%
    rename(id.exposure=X1, id.outcome=X2, I2gx=X3) %>%
    mutate(I2gx=100*I2gx) %>%
    mutate(method="MR Egger")

df <- data.frame(matrix(unlist(mrpresso), nrow=length(mrpresso), byrow=T),stringsAsFactors=FALSE)
mr.presso <- df %>%
    rename(b=X1, se=X2, t=X3, p=X4, no_b=X5, no_se=X6, no_t=X7, no_p=X8, global.p=X9, bias.p=X10, snp.total=X11, snp.outliers=X12, snp.list=X13)

d <- ivw %>%
    mutate(snp.total=nsnp) %>%
    select(b, snp.total, exposure, outcome, id.exposure, id.outcome)

mr.presso <- merge(d, mr.presso, all.x = TRUE, all.y = TRUE)
mr.presso <- mr.presso %>%
    mutate(pval=p) %>%
    mutate(method="MR-PRESSO") %>%
    select(id.exposure, id.outcome, exposure, outcome, method, b, se, pval, no_b, no_se, global.p, bias.p, snp.total, snp.outliers, snp.list )

mr.presso[,7:14] <- lapply(mr.presso[,7:14], as.numeric)
mr.presso <- mr.presso %>%
    mutate(nsnp=snp.total-snp.outliers)

df <- data.frame(matrix(unlist(mrpresso), nrow=length(mrpresso), byrow=T),stringsAsFactors=FALSE)
mrpresso.steiger <- df %>%
    rename(b=X1, se=X2, t=X3, p=X4, no_b=X5, no_se=X6, no_t=X7, no_p=X8, global.p=X9, bias.p=X10, snp.total=X11, snp.outliers=X12, snp.list=X13)

mr_res <- merge(res, plt, all.x=TRUE, all.y=TRUE)
mr_res <- merge(mr_res, mregger.i2gx, all.x=TRUE, all.y=TRUE)
mr_res <- merge(mr_res, het, all.x=TRUE, all.y=TRUE)
mr_res <- bind_rows(mr_res, mr.presso, mrradial)

mr_res <-mr_res[order(mr_res$id.exposure, mr_res$id.outcome, mr_res$method),]

mr_res <- mr_res[,c(5, 4, 3, 7:9, 14, 16, 17,  10, 13, 12, 6, 23, 24)]

harm_dat <- harm_dat %>%
    mutate(r.exposure=NA) %>%
    mutate(r.outcome=NA)

harm_dat$r.exposure <- get_r_from_lor(
    lor=harm_dat$beta.exposure,
    af=harm_dat$eaf.exposure,
    ncase=harm_dat$ncase.exposure,
    ncontrol=harm_dat$ncontrol.exposure,
    prevalence=harm_dat$prevalence.exposure,
    model = "logit",
    correction = FALSE
)

harm_dat$r.outcome <- get_r_from_lor(
    lor=harm_dat$beta.outcome,
    af=harm_dat$eaf.outcome,
    ncase=harm_dat$ncase.outcome,
    ncontrol=harm_dat$ncontrol.outcome,
    prevalence=harm_dat$prevalence.outcome,
    model = "logit",
    correction = FALSE
)

harm_dat <- harm_dat %>%
    select(beta.exposure, beta.outcome, eaf.exposure, eaf.outcome,
id.outcome, se.outcome,, outcome, units.outcome, pval.outcome, ncase.outcome,
ncontrol.outcome, prevalence.outcome, se.exposure, id.exposure, exposure, units.exposure, pval.exposure,
prevalence.exposure, ncontrol.exposure, ncase.exposure, r.exposure, r.outcome)

harm_dat <- harm_dat %>%
    mutate(samplesize.exposure=650000, samplesize.outcome=750000) %>%
    select(-starts_with("pval"))


out <- directionality_test(harm_dat) %>%
    select(exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction)
mr_res<- merge(out, mr_res, all.x=TRUE, all.y=TRUE)


reverse <- generate_odds_ratios(mr_res)
names(reverse)
reverse <- reverse[,c(1, 2, 6, 21:23, 9:13, 15, 14, 16:18, 3:5)]
save(reverse, file = "data/reverse.RData")
