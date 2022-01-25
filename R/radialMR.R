library(TwoSampleMR)
library(dplyr)
library(tidyverse)
library(RadialMR)
vignette("RadialMR")

load("data/harm_dat.RData")

harm_radial <- harm_dat %>%
    filter(id.exposure!=4) %>%
    filter(!(id.outcome>20 & id.exposure==8)) #for some reason there are only 4!!! hyperthyroidism (exposure=8) SNPs in Neale lab data, and MRradial does not work

radial_result <- harm_dat %>% filter(id.exposure==-1)
radial_outliers <- harm_dat %>% filter(id.exposure==-1)
d <- subset(harm_radial, !duplicated(paste(id.exposure, " - ", id.outcome)),
            select = c(exposure, outcome, id.exposure, id.outcome))
for (j in 1:nrow(d)) {
    x <- subset(harm_radial, exposure == d$exposure[j] & outcome == d$outcome[j])
    message(x$exposure[1], " - ", x$outcome[1])
    radial_data <-format_radial(x$beta.exposure,x$beta.outcome,
                                                 x$se.exposure,x$se.outcome,
                                                 x$SNP)
    ivw_model <- ivw_radial(radial_data,0.05/nrow(radial_data),3,0.0001)
#    message(typeof(ivw_model$outliers))
    if (is.character(ivw_model$outliers)) {
        keep <- subset(harm_radial, exposure == d$exposure[j] & outcome == d$outcome[j])
    } else {
        keep <- subset(harm_radial, exposure == d$exposure[j] & outcome == d$outcome[j] & !SNP %in% ivw_model$outliers$SNP)
        radial_outliers <- rbind(radial_outliers, subset(harm_radial, exposure == d$exposure[j] & outcome == d$outcome[j] & SNP %in% ivw_model$outliers$SNP))
    }
    radial_result <- rbind(radial_result, keep)
}

results <- mr(radial_result) %>%
    filter(method=="Wald ratio" | method=="Inverse variance weighted") %>%
    mutate(method="IVW, -RadialMR outliers")

radial_outliers1 <- radial_outliers %>%
    select(SNP, exposure, outcome) %>%
    group_by(exposure, outcome) %>%
    summarise(snp.outliers=n()) %>%
    ungroup()

radial_outliers2 <- radial_outliers %>%
    select(SNP, exposure, outcome) %>%
    group_by(exposure, outcome) %>%
    summarise(snp.list = paste0(SNP, collapse = " ")) %>%
    ungroup()

radial_outliers <- left_join(radial_outliers1, radial_outliers2)


#HETREOGENEITY
het<-mr_heterogeneity(radial_result)
het <- het %>%
    mutate(I2=100*(Q-Q_df)/Q) %>%
    mutate(I2=if_else(I2<0, 0, I2)) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/
het <- subset(het, method=="Inverse variance weighted") %>%
    mutate(method="IVW, -RadialMR outliers")


#Merge results
mrradial <- merge(results, het, by=c("id.exposure", "id.outcome", "outcome", "exposure", "method"), all.x = TRUE, all.y = T) #check for column names
names(mrradial)
mrradial_all <- merge(mrradial, radial_outliers, all.x=TRUE)

mrradial_all$snp.outliers <-as.numeric(mrradial_all$snp.outliers)

mrradial <- mrradial_all %>%
    mutate(snp.outliers = if_else(is.na(snp.outliers), 0, snp.outliers))

save(mrradial, file = "data/mrradial_thyroid_CBC.RData")
