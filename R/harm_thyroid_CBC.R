library(readr)
library(readxl)
library(openxlsx)
library(TwoSampleMR)
library(dplyr)
library(tidyverse)

#################################################
#EXPOSURE: TSH
#
#################################################
TSH_TOC <- read_excel("data-raw/MR_overview_instruments_March2021_RS.xlsx",
                                                   sheet = "TSH_ThyroidOmicsConsortium",
                                                   skip = 2)


TSH_TOC <- TSH_TOC %>%
  mutate(var=paste(Chr, Pos, sep=":")) %>%
  mutate(id.exposure=1) %>%
  mutate(exposure="TSH_TOC")   %>%
  mutate(F.stat=(Effect_TSH*Effect_TSH)/(StdErr_TSH*StdErr_TSH)) %>% #F: 30-536
  filter(F.stat>9.99) %>%
  filter(RSID!="rs12390237") #X-chromosome variant, ABO and BCAS3 removed already
names(TSH_TOC)

dat_TSH_TOC <- data.frame(SNP = TSH_TOC$RSID, beta.exposure= TSH_TOC$Effect_TSH, se.exposure= TSH_TOC$StdErr_TSH, id.exposure=TSH_TOC$id.exposure, exposure=TSH_TOC$exposure,
                              units.exposure="SD", pval.exposure=TSH_TOC$Pvalue_TSH, samplesize.exposure=TSH_TOC$Samplesize_TSH, var=TSH_TOC$var, gene=TSH_TOC$GeneName,
                              effect_allele.exposure=toupper(TSH_TOC$Effect_allele), other_allele.exposure=toupper(TSH_TOC$Other_allele),  eaf.exposure= TSH_TOC$EAF_TSH)


TSH_HUNT <- read_excel("data-raw/MR_overview_instruments_March2021_RS.xlsx",
                                                   sheet = "TSH_HUNT-MGI-TOC", skip = 2)

TSH_HUNT <- TSH_HUNT %>%
  mutate(var=paste(Chr, Pos, sep=":")) %>%
  mutate(id.exposure=3) %>%
  mutate(exposure="TSH_HUNT") %>%
  mutate(gene=`GeneName*`) %>%
  mutate(F.stat=(EffectAlternate_normalrangeTSH_Hunt*EffectAlternate_normalrangeTSH_Hunt)/(SE_normalrangeTSH_Hunt*SE_normalrangeTSH_Hunt)) %>%
  filter(F.stat!='NA') %>%
  filter(RSID!="rs546738875" & RSID!="rs145153320") #MAF<0.005 and table states 28, where 30 are given

TSH_HUNT <- TSH_HUNT %>%
  filter(F.stat>9.99) #5/28 SNPs have F<10

names(TSH_HUNT)
dat_TSH_HUNT <- data.frame(SNP = TSH_HUNT$RSID, beta.exposure= TSH_HUNT$EffectAlternate_normalrangeTSH_Hunt, se.exposure= TSH_HUNT$SE_normalrangeTSH_Hunt, id.exposure=TSH_HUNT$id.exposure, exposure=TSH_HUNT$exposure,
                          units.exposure="SD", pval.exposure=TSH_HUNT$`Pvalue_normalrangeTSH_Hunt - survived FDR<5%`, samplesize.exposure=TSH_HUNT$Samplesize_normalrangeTSH_HUNT, var=TSH_HUNT$var, gene=TSH_HUNT$gene,
                          effect_allele.exposure= TSH_HUNT$Effect_allele, other_allele.exposure= TSH_HUNT$Other_allele,  eaf.exposure= TSH_HUNT$EAF_normalrangeTSH_Hunt)

#binding TSH_TOC and TSH_HUNT
TSH <- bind_rows(dat_TSH_TOC, dat_TSH_HUNT) %>%
  mutate(exposure="TSH_all", id.exposure=1)


#FT4
FT4<- read_excel("data-raw/MR_overview_instruments_March2021_RS.xlsx",
                                                   sheet = "fT4_ThyroidOmicsConsortium",
                                                   skip = 2)

FT4_TOC <- FT4 %>%
  mutate(var=paste(Chr, Pos, sep=":")) %>%
  mutate(id.exposure=2) %>%
  mutate(exposure="FT4_TOC") %>%
  filter(RSID!='NA')
names(FT4_TOC)
FT4_TOC <- data.frame(SNP = FT4_TOC$RSID, beta.exposure= FT4_TOC$Effect_TSH, se.exposure= FT4_TOC$StdErr_TSH, id.exposure=FT4_TOC$id.exposure, exposure=FT4_TOC$exposure,
                          units.exposure="SD", pval.exposure=FT4_TOC$Pvalue_TSH, samplesize.exposure=FT4_TOC$Samplesize_fT4, var=FT4_TOC$var, gene=FT4_TOC$GeneName,
                          effect_allele.exposure=toupper(FT4_TOC$Effect_allele), other_allele.exposure=toupper(FT4_TOC$Other_allele),  eaf.exposure= FT4_TOC$EAF_TSH)


FT4_DIO <- FT4_TOC[FT4_TOC$gene %in% c("DIO1","DIO2"), ] %>%
  mutate(id.exposure=3, exposure="FT4_DIO1_2")

FT4_DIO3OS <- FT4_TOC[FT4_TOC$gene %in% c("DIO3OS"), ] %>%
  mutate(id.exposure=4, exposure="FT4_DIO3OS")

AITD <- read_excel("data-raw/MR_overview_instruments_March2021_RS.xlsx",
                                                   sheet = "Hypothyroidism_Decode-UKBB",
                                                   skip = 2)

AITD <-AITD %>%
  mutate(var=paste(Chr, `Position(hg19)`, sep=":")) %>%
  mutate(id.exposure=5) %>%
  mutate(exposure="AITD")  %>%
  filter(!is.na(Effect_allele)) %>%
  filter(is.na(fjern))


names(AITD)
dat_AITD <- data.frame(SNP = AITD$RSID, beta.exposure= AITD$Effect_hypothyroidism, se.exposure= AITD$`StdErr_hypothyroidism*`, id.exposure=AITD$id.exposure, exposure=AITD$exposure,
                          units.exposure="log odds", pval.exposure=AITD$Pvalue_hypothyroidism, ncase.exposure=AITD$Samplesize_UKBB_Decode, ncontrol.exposure=725172, prevalence.exposure=0.038, var=AITD$var, gene=AITD$GeneName,
                          effect_allele.exposure= AITD$Effect_allele, other_allele.exposure= AITD$Other_allele,  eaf.exposure= AITD$`EAF_hypothyroidism*`)

dat_AITD$beta.exposure <-as.numeric(dat_AITD$beta.exposure)
dat_AITD$se.exposure <-as.numeric(dat_AITD$se.exposure)
dat_AITD$pval.exposure <-as.numeric(dat_AITD$pval.exposure)

dat_AITD <- dat_AITD %>%
  mutate(F.stat=(beta.exposure*beta.exposure)/(se.exposure*se.exposure)) %>%
  filter(F.stat>9.99) %>%
  filter(SNP!="rs146750254" & SNP!="rs11052754" & SNP!="rs76226393")  #3 SNPs to be excluded according to updated file from Rosalie
#check rs146750254 and rs11052754 betas and SEs roughly same size, F<2

redundant_hypo <- read_excel("data-raw/MR_overview_instruments_March2021_RS.xlsx",
                                                   sheet = "redundant_AITD") #SNPs with same chromosomal location, but different rs IDs
`%notin%` <- Negate(`%in%`)
dat_AITD <- dat_AITD %>%
  filter(SNP %notin% redundant_hypo$redundant)
dat_AITD <- dat_AITD[!duplicated(dat_AITD$var), ]


AITD <- dat_AITD %>%
  filter(nchar(other_allele.exposure)==1) #OA has to have exactly 1 character, i.e., be a SNP

#"Old" thyroid exposures: overt hypo from 23andMe, and subclin hypo and hyper from Teumer
dat_old <- read_excel("data-raw/old_thyroid_exposures.xlsx")


#Subclinical hypothyroidism, Teumer id.exposure=6
dat_SCH <- dat_old %>%
  filter(Exposure=="hypothyroidism, subclinical")
SCH <- data.frame(gene=dat_SCH$`Nearest gene`, SNP=dat_SCH$`SNP rs ID`, var=dat_SCH$var, beta.exposure=dat_SCH$Beta, id.exposure=6, effect_allele.exposure=dat_SCH$`Effect allele`, se.exposure=dat_SCH$`Standard error`,
                    units.exposure="log odds", eaf.exposure=dat_SCH$eaf, other_allele.exposure=dat_SCH$`Other allele`,
                  ncase.exposure=3440, ncontrol.exposure=49983, prevalence.exposure=0.0644, pval.exposure=dat_SCH$p, exposure="SCH")

#Overt hypothyroidism, 23andMe id.exposure=7
dat_ohypo <- dat_old %>%
  filter(Exposure=="hypothyroidism, overt")
ohypo <- data.frame(gene=dat_ohypo$`Nearest gene`, SNP=dat_ohypo$`SNP rs ID`, var=dat_ohypo$var, beta.exposure=dat_ohypo$Beta, id.exposure=7, effect_allele.exposure=dat_ohypo$`Effect allele`, se.exposure=dat_ohypo$`Standard error`,
                  units.exposure="log odds", eaf.exposure=dat_ohypo$eaf, other_allele.exposure=dat_ohypo$`Other allele`, ncase.exposure=17558, ncontrol.exposure=117083, prevalence.exposure=0.13,  pval.exposure=dat_ohypo$p, exposure="hypothyroidism, overt")


# SUbclinical hyperthyroidism, id.exposure=8
dat_hyper <- dat_old %>%
  filter(Exposure=="hyperthyroidism, subclinical")

hyper <- data.frame(gene=dat_hyper$`Nearest gene`, SNP=dat_hyper$`SNP rs ID`, var=dat_hyper$var, beta.exposure=dat_hyper$Beta, id.exposure=8, effect_allele.exposure=dat_hyper$`Effect allele`, se.exposure=dat_hyper$`Standard error`,
                    units.exposure="log odds", eaf.exposure=dat_hyper$eaf, other_allele.exposure=dat_hyper$`Other allele`, ncase.exposure=1840, ncontrol.exposure=49983, prevalence.exposure=0.0355, pval.exposure=dat_hyper$p, exposure="hyperthyroidism")

hyper <- hyper %>%
  filter(SNP!="rs925488") #Remove FOXE1 variant associated with both hypo-and hyperthyroidism


#################################################
#Combining EXPOSURES, EA
#################################################
SCH$eaf.exposure <-as.numeric(SCH$eaf.exposure)
ohypo$eaf.exposure <-as.numeric(ohypo$eaf.exposure)
hyper$eaf.exposure <-as.numeric(hyper$eaf.exposure)

dat_exposure <- bind_rows(TSH, FT4_TOC, FT4_DIO, FT4_DIO3OS, AITD, SCH, ohypo, hyper) %>%
  mutate(F.stat=round((beta.exposure*beta.exposure)/(se.exposure*se.exposure)))


#################################################
#################################################
#OUTCOME: Neale lab's reticulocytes and bilirubin
#################################################
CBC <- read_delim("data-raw/CBC_Neale-20210909.tsv",
                   "\t", escape_double = FALSE, trim_ws = TRUE)

spec(CBC)


dat_CBC <- data.frame(SNP = CBC$rs, beta.outcome= CBC$beta, se.outcome= CBC$seoutcome,  outcome=CBC$outcome_desc, #id.exposure=CBC$id.exposure,
                       units.outcome="SD", pval.outcome=CBC$pval, samplesize.outcome=CBC$n_complete_samples, var=CBC$varid, #gene.exposure="Gene", mr
                       effect_allele.outcome=CBC$alt, other_allele.outcome=CBC$ref,   eaf.outcome= CBC$AF)
dat_CBC <- dat_CBC %>% relocate(var)
dat_CBC$var <- gsub("[_, A, C, T, G]","", dat_CBC$var)

dat_CBC <- dat_CBC %>%
   filter(SNP %in% dat_exposure$SNP)
dat_CBC <- dat_CBC %>%
  group_by(outcome) %>%
  dplyr::mutate(id.outcome = cur_group_id()) %>%
  ungroup()

dat_CBC <- dat_CBC %>%
    mutate(id.outcome=id.outcome+20)
dat_CBC  <- dat_CBC  %>%
  filter(var %in% dat_exposure$var) %>% distinct()
#################################################
#OUTCOME:BCX, outcomes 1-15
#################################################
BCX <- read.delim("data-raw/BCX2_20210909.tsv")
names(BCX)

BCX <- data.frame(var = BCX$rs_number, beta.outcome= BCX$beta, se.outcome= BCX$se, outcome=BCX$outcome,
                  units.outcome="SD", pval.outcome=BCX$pvalue, samplesize.outcome=BCX$n_samples,  #gene.exposure="Gene", mr
                  effect_allele.outcome=BCX$reference_allele, other_allele.outcome=BCX$other_allele,   eaf.outcome= BCX$af)
#gsub(pattern, replacement, x)
BCX$var <- gsub('[_, A, C, T, G]','', BCX$var)


BCX <- BCX %>%
  filter(var %in% dat_exposure$var)

dat <- dat_exposure %>%
  select(var, SNP)
dat_BCX <- left_join(BCX, dat)

dat_BCX <- dat_BCX %>%
  group_by(outcome) %>%
  dplyr::mutate(id.outcome = cur_group_id()) %>% #assigned id.outcome: 1-15
  ungroup()

#Somehow some of the SNPs in the outcomes are duplicates (e.g. rs3775291 C/T, C/G?? for MPV), but different summary stats, so keeping the ones with highest sample size, as harmonizing data does not solve the problem
#3668 ->3645
dat_BCX <- dat_BCX %>%
  filter(var %in% dat_exposure$var) %>%
  group_by(SNP, outcome) %>%
  filter(samplesize.outcome==max(samplesize.outcome)) %>%
  ungroup()

#################################################
#OUTCOME:Pernicious anemia
#################################################
pernicious_anemia <- read_delim("data-raw/b12_meta_ukbb_estbb_finngen_061020.out",
                                                 "\t", escape_double = FALSE, trim_ws = TRUE)
names(pernicious_anemia)

dat_pa <- data.frame(SNP = pernicious_anemia$rs_number, beta.outcome= pernicious_anemia$beta, se.outcome= pernicious_anemia$se, outcome="pernicious_anemia", id.outcome=17,
                  units.outcome="log odds", prevalence.outcome=0.0033, ncontrol.outcome=pernicious_anemia$n_samples,  ncase.outcome=round(0.0033*pernicious_anemia$n_samples),
                  pval.outcome=pernicious_anemia$'p-value', #gene.exposure="Gene", mr
                  effect_allele.outcome=pernicious_anemia$reference_allele, other_allele.outcome=pernicious_anemia$other_allele,   eaf.outcome= pernicious_anemia$eaf)

dat_pa <- dat_pa %>%
  filter(SNP %in% dat_exposure$SNP)
#################################################
#Combining OUTCOMES
#################################################
dat_outcome <- bind_rows(dat_BCX, dat_pa, dat_CBC) %>% distinct()


#################################################
#HARMONISING DATA
#################################################
harm_dat <- harmonise_data(dat_exposure, dat_outcome) %>%
  filter(mr_keep=="TRUE") #Removing SNPs for being palindromic with intermediate allele frequencies
list <- harm_dat %>% select(id.outcome, outcome) %>% distinct()

harm_dat <- harm_dat %>%
  mutate(pval.outcome=if_else(pval.outcome==0, 1e-300, pval.outcome))   %>%
  filter(id.outcome!=1 & id.outcome!=2 & id.outcome!=5 & id.outcome!=9 & id.outcome!=10 & id.outcome!=11 & id.outcome!=12  & id.outcome!=15)
list <- harm_dat %>% select(id.outcome, outcome) %>% distinct()

save(harm_dat, file = "data/harm_dat.RData")
#################################################
#Steiger filtering
#################################################
harm_dat.steiger<-steiger_filtering(harm_dat) %>%
  filter(steiger_dir=="TRUE")
save(harm_dat.steiger, file = "data/harm_dat_sf.RData")
