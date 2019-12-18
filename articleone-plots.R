mycolors <- function(i) {
  colors <- list("red" = "#FF7575",
                 "blue" = "#4B81B6",
                 "gray" = "#bebebe",
                 "white" = "#FFFAFA",
                 "purewhite" = "#FFFFFF")
  get(i, colors)
}

mzrt2mz <- function(mzid) {
  mzid %>%
    gsub("mzid_", "", .) %>%
    gsub("_.*", "", .) %>%
    as.numeric()
}

plot.manhattanplot <- function(linres,
                               bonfp = NULL,
                               response = "SBP",
                               nlabels = 10,
                               xlab = "Eicosanoids Rank Ordered by Mass to Charge Ratio",
                               ylab = "Negative Log P Value",
                               color_scheme = c("positive" = mycolors("red"),
                                                "negative" = mycolors("blue"),
                                                "insignificant" = mycolors("gray"))) {
    dat <- linres %>%
    dplyr::filter(response == "SBP", grepl('^mzid', term)) %>%
    dplyr::arrange(term) %>%
    dplyr::mutate(neg_log10_pvalue = -1 * log10(p.value),
                  coloring = ifelse(p.value > bonfp, "insignificant",
                             ifelse(estimate < 0, "negative", "positive")),
                  mz = mzrt2mz(term))

  label_dat <- dat %>%
    dplyr::arrange(desc(neg_log10_pvalue)) %>%
    dplyr::slice(seq_len(nlabels)) %>%
    dplyr::select(term, neg_log10_pvalue, mz, mzrt, coloring)

  ggplot2::ggplot(dat, ggplot2::aes(x = mz, y = neg_log10_pvalue)) +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::geom_point(ggplot2::aes(color = coloring), size = 1.5) +
    ggrepel::geom_text_repel(data = label_dat,
                             ggplot2::aes(label = mzrt),
                             segment.size = 0.1,
                             size = 3) +
    ggplot2::geom_hline(yintercept = -1 * log10(bonfp), linetype = 2) +
    ggplot2::scale_x_continuous(breaks = seq(200, 650, 50)) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0)) +
    ggplot2::theme(
      text = element_text(size = 14),
      legend.position="none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.background = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = color_scheme)
}

spearmancorrelation  <- function(data, mzids) {
  dat <- data$data.full %>% dplyr::select(mzids)
  colnames(dat) <- colnames(dat) %>% pub.mzrt(., data$metabolite.names)
  cor <- cor(dat, method = 'spearman')
}

replicationforestplot <- function(forest.fr02, forest.fhs, file) {
    png(width = 1600, height = 700, res = 150, file = file)
    forestplot::forestplot(
                    labeltext = cbind(c("Eicosanoid", "", forest.fhs$name),
                                      c("FINRISK\n β (95% CI)", "", forest.fr02$mean_ci),
                                      c("FHS\n β (95% CI)", "", forest.fhs$mean_ci)),
                    mean = cbind(c(NA, NA, forest.fr02$estimate),
                                 c(NA, NA, forest.fhs$estimate)),
                    lower = cbind(c(NA, NA, forest.fr02$conf.low),
                                  c(NA, NA, forest.fhs$conf.low)),
                    upper = cbind(c(NA, NA, forest.fr02$conf.high),
                                  c(NA, NA, forest.fhs$conf.high)),
                    legend = c("FINRISK", "FHS"),
                    align = c("l", "l", "l"),
                    graph.pos = 4,
                    title = "",
                    xlog = FALSE,
                    xlab = " β (95%-CI)",
                    txt_gp = fpTxtGp(label = gpar(cex = 1.25),
                                     ticks = gpar(cex = 1.1),
                                     xlab = gpar(cex = 1.2),
                                     title = gpar(cex = 1.2)),
                    xticks = seq(-1, 4),
                    clip =c(-1, 4),
                    col = fpColors(box=c("blue", "darkred")),
                    zero = 0, 
                    lineheight = unit(14, "mm"), 
                    boxsize = 0.2, 
                    colgap = unit(8, "mm"),
                    lwd.ci = 1)
    dev.off()
}

plot.riskmodel <- function(riskmodel.fhs,
                           riskmodel.fr02,
                           variable,
                           file,
                           ylab = "Odds for hypertension (95%-CI)") {
    riskmodel <- rbind(riskmodel.fhs %>% mutate(cohort = "FHS"),
                       riskmodel.fr02 %>% mutate(cohort = "FINRISK") %>% select(-fo, -loop)) %>%
        filter(response == variable) %>%
        add_row(term = "fwdbonf_riskclass1", adjustment = c("adjusted", "unadjusted"),
                cohort = "FHS", estimate = 1, conf.low = 1, conf.high = 1) %>%
        add_row(term = "fwdbonf_riskclass1", adjustment = c("adjusted", "unadjusted"),
                cohort = "FINRISK", estimate = 1, conf.low = 1, conf.high = 1) %>%
        dplyr::mutate(adjustment = relevel(as.factor(adjustment), ref = "unadjusted"))
    
    g.riskforest <- ggplot(riskmodel, aes(x = term,
                                          y = estimate,
                                          ymin = conf.low,
                                          ymax = conf.high,
                                          col = cohort)) +
        facet_wrap(~adjustment,
                   strip.position = "top", 
                   labeller = as_labeller(c("unadjusted" = "Unadjusted",
                                            "adjusted" = "Multivariable-adjusted"))) +
        ylab(ylab) +
        xlab("Eicosanoid risk score") +
        geom_pointrange(position = position_dodge(width = 0.4), shape = 20) +
        geom_hline(yintercept = 1, linetype="dotted") +
        scale_colour_manual(name = "Cohorts",
                            labels = c("FINRISK", "FHS"),
                            breaks = c("FINRISK", "FHS"),
                            values = c("blue", "red")) +
        scale_x_discrete(labels=c("fwdbonf_riskpersd" = "per\n1-SD",
                                  "fwdbonf_riskclass1" = "Q1",
                                  "fwdbonf_riskclass2" = "Q2",
                                  "fwdbonf_riskclass3" = "Q3",
                                  "fwdbonf_riskclass4" = "Q4")) +
        scale_y_continuous(breaks = pretty_breaks()) +
        theme_classic() +
        theme(strip.background = element_blank(),
              strip.placement = "outside")
    
    ggsave(file = file, plot = g.riskforest, height = 3, width = 6, dpi = 300)
}
