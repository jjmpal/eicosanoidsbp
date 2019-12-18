library("turkumetabolites")
library(dplyr)
#library(ggplot2)
#library(broom)
#library(tidyr)
#library(gridExtra)
require(survminer)

left <- turkumetabolites::checkcoxassumption(
  rset = combined.risk.set, 
  forcedcovariates = "AGE + sex + tot_chol + hdl_chol + curr_smk + curr_diab", 
  singlecovariates = loopedcovariates,
  events = list(model_1 = "cvdth", model_2 = "dth"))

pdf("checkproporationalhazard-left.pdf", paper="a4")
for(outer in names(left)) {
  for(inner in names(left[[outer]])) {
    survminer::ggcoxzph(left[[outer]][[inner]], 
                        newpage = FALSE, 
                        caption = paste(outer, inner),
                        ggtheme = theme(plot.margin = unit(c(1,1,1,1), "mm"),
                                        text = element_text(size=6))) %>% 
      print
  }
}
dev.off()

right <- coxmodels.for.paneltwo <- turkumetabolites::checkcoxassumption(
  rset = combined.risk.set,
  forcedcovariates = "AGE + sex + tot_chol + hdl_chol + curr_smk + curr_diab + HRX + MAP", 
  singlecovariates = loopedcovariates,
  events = list(model_1 = "cvdth", model_2 = "dth"))

pdf("checkproporationalhazard-right.pdf", paper="a4")
for(outer in names(right)) {
  for(inner in names(right[[outer]])) {
    survminer::ggcoxzph(right[[outer]][[inner]], 
                        newpage = FALSE, 
                        caption = paste(outer, inner),
                        ggtheme = theme(plot.margin = unit(c(1,1,1,1), "mm"),
                                        text = element_text(size=6))) %>% 
      print
  }
}
dev.off()
