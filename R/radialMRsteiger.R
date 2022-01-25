library(TwoSampleMR)
library(dplyr)
library(tidyverse)
library(RadialMR)
vignette("RadialMR")

load("data/harm_dat_sf.RData")

harm_radial.steiger <- harm_dat.steiger %>%
    filter(id.exposure!=4) %>% #id.exposure 4 is FT4 DIO3OS with 2 SNPs
    filter(!(id.outcome>20 & id.exposure==8)) #for some reason there are only 4!!! hyperthyroidism (exposure=8) SNPs in Neale lab data, and MRradial does not work

radial_result.steiger <- harm_dat.steiger %>% filter(id.exposure==-1)
radial_outliers.steiger <- harm_dat.steiger %>% filter(id.exposure==-1)
d <- subset(harm_radial.steiger, !duplicated(paste(id.exposure, " - ", id.outcome)),
            select = c(exposure, outcome, id.exposure, id.outcome))
for (j in 1:nrow(d)) {
    x <- subset(harm_radial.steiger, exposure == d$exposure[j] & outcome == d$outcome[j])
    message(x$exposure[1], " - ", x$outcome[1])
    radial_data.steiger <-format_radial(x$beta.exposure,x$beta.outcome,
                                x$se.exposure,x$se.outcome,
                                x$SNP)
    ivw_model <- ivw_radial(radial_data.steiger,0.05/nrow(radial_data.steiger),3,0.00001)
    #    message(typeof(ivw_model$outliers))
    if (is.character(ivw_model$outliers)) {
        keep <- subset(harm_radial.steiger, exposure == d$exposure[j] & outcome == d$outcome[j])
    } else {
        keep <- subset(harm_radial.steiger, exposure == d$exposure[j] & outcome == d$outcome[j] & !SNP %in% ivw_model$outliers$SNP)
        radial_outliers.steiger <- rbind(radial_outliers.steiger, subset(harm_radial.steiger, exposure == d$exposure[j] & outcome == d$outcome[j] & SNP %in% ivw_model$outliers$SNP))
    }
    radial_result.steiger <- rbind(radial_result.steiger, keep)
}

results.steiger <- mr(radial_result.steiger) %>%
    filter(method=="Wald ratio" | method=="Inverse variance weighted") %>%
    mutate(method="IVW, -RadialMR outliers")

radial_outliers1 <- radial_outliers.steiger %>%
    select(SNP, exposure, outcome) %>%
    group_by(exposure, outcome) %>%
    summarise(snp.outliers=n()) %>%
    ungroup()

radial_outliers2 <- radial_outliers.steiger %>%
    select(SNP, exposure, outcome) %>%
    group_by(exposure, outcome) %>%
    summarise(snp.list = paste0(SNP, collapse = " ")) %>%
    ungroup()

radial_outliers.steiger <- left_join(radial_outliers1, radial_outliers2)


#HETREOGENEITY
het.steiger<-mr_heterogeneity(radial_result.steiger)
het.steiger <- het.steiger %>%
    mutate(I2=100*(Q-Q_df)/Q) %>%
    mutate(I2=if_else(I2<0, 0, I2)) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/
het.steiger <- subset(het.steiger, method=="Inverse variance weighted") %>%
    mutate(method="IVW, -RadialMR outliers")


#Merge results
mrradial.steiger <- merge(results.steiger, het.steiger, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x = TRUE, all.y = T) #check for column names
names(mrradial.steiger)
mrradial_all <- merge(mrradial.steiger, radial_outliers.steiger, all.x=TRUE)

mrradial_all$snp.outliers <-as.numeric(mrradial_all$snp.outliers)

mrradial.steiger <- mrradial_all %>%
    mutate(snp.outliers = if_else(is.na(snp.outliers), 0, snp.outliers))

save(mrradial.steiger, file = "data/mrradial_s_thyroid_CBC.RData")
