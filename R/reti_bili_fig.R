library(TwoSampleMR)
library(dplyr)
library(forestplot)


load("data/harm_dat.RData")
load("data/harm_dat_sf.RData")
##############################
#MR analyses
##############################
harm_dat <- harm_dat %>%
    mutate(pval.outcome=if_else(pval.outcome==0, 1e-300, pval.outcome))   %>%
    filter(id.exposure==3 |id.exposure==5 | id.exposure==7) %>%
    filter(id.outcome==30 | id.outcome==33)


res <- mr(harm_dat, method_list=c( "mr_ivw", "mr_weighted_median", "mr_egger_regression"))
ivw <- subset_on_method(res,
                        single_snp_method = "Wald ratio",
                        multi_snp_method = "Inverse variance weighted"
)

fig.1 <- ivw %>%
    filter(id.exposure==3 |id.exposure==5 | id.exposure==7) %>%
    filter(id.outcome==30 | id.outcome==33)



fig.1 <-fig.1 %>%
    mutate(beta=sprintf("%5.3f", b)) %>%
    mutate(lower=b-1.96*se) %>%
    mutate(upper=b+1.96*se) %>%
    mutate(lci=sprintf("%5.3f", lower)) %>%
    mutate(uci=sprintf("%5.3f", upper)) %>%
    mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA,NA,5,7,3,NA,5,7,3)
id.outcome <- c(NA,NA, 30,30,30,NA, 33,33,33)
id <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
panel.fig1 <- data.frame(id.exposure, id.outcome, id)

fig.1 <-merge(panel.fig1, fig.1, all.x = TRUE, all.y = TRUE)
fig.1 <-fig.1[order(fig.1$id),]


text.fig1<-cbind(
    c("Thyroid exposures:", "Outcome: Reticulocyte count","AITD", "Hypothyroidism, overt", "FT4 regulated by DIO1/DIO2",   "Outcome: Total bilirubin",  "AITD", "Hypothyroidism, overt", "FT4 regulated by DIO1/DIO2"),
    c(fig.1$est.size))
text.fig1[1, 2] <- "\u03b2 (95% CI)"


#creating plot
pdf(file = 'reti_bili.pdf', onefile=F)
reti_bili<-forestplot(text.fig1,
               fn.ci_norm = fpDrawCircleCI,
               txt_gp=fpTxtGp(ticks=gpar(cex=0.8),
                              xlab=gpar(cex=0.8),
                              label=gpar(cex=0.8)),
               title="Thyroid function, reticulocytes and bilirubin",
               graph.pos = 2,
               mean=cbind(fig.1$b),
               upper =cbind(fig.1$upper),
               lower= cbind(fig.1$lower),
               lwd.ci = 1,
               vertices=TRUE,
               boxsize = 0.2,
               #clip=c(-0.075, 1.75),
               is.summary=c(TRUE,TRUE, rep(FALSE,3),TRUE, rep(FALSE,3)),
               col = fpColors(box = "black",
                              line = "black"),
               xticks = c(-0.1, -0.05, 0, 0.05, 0.1),
              # zero=0,
              xlab=expression(paste(beta, " (95% CI)")))
dev.off()
