compare.names <- function() {
    tribble(~term, ~name,
            "mzid_339.179900_1.7205", "11-dehydro-2,3-dinor-TXB2",
            "mzid_279.196600_3.7247", "12-HHTrE",
            "mzid_295.227900_4.8902", "295.2279/4.89 (putative eicosanoid)",
            "mzid_319.228000_5.6733", "319.2280/5.67 (unknown)",
            "mzid_331.267900_6.5274", "Adrenic acid",
            "mzid_265.181000_3.5705", "265.1810/3.57 (putative eicosanoid)")
}

compare.fhs <- function() {
    tribble(
        ~term, ~estimate, ~std.error, ~statistic, ~p.value, ~conf.low, ~conf.high,
        "mzid_319.228000_5.6733", 1.420476819, 0.305985562, 4.642300149, 3.60E-06, 0.820756138, 2.0201975,
        "mzid_279.196600_3.7247", 1.375545359, 0.307302607, 4.476191641, 7.90E-06, 0.773243318, 1.9778474,
        "mzid_331.267900_6.5274", 1.258063711, 0.319231395, 3.940914746, 8.31E-05, 0.632381674, 1.883745747,
        "mzid_265.181000_3.5705", 1.885710794, 0.306688267, 6.148623854, 8.91E-10, 1.284612836, 2.486808752) %>%
        full_join(., compare.names(), by = c("term" = "term")) %>%
        mutate(mean_ci = ifelse(is.na(estimate),
                                NA,
                                sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)),
               mzrt = pub.mzrt.naive(term)) %>%
        arrange(name)
}

riskmodel.results.fhs <- function() {
    tribble(
        ~adjustment, ~response, ~term, ~estimate, ~std.error, ~statistic, ~p.value, ~conf.low, ~conf.high,
        "adjusted", "HT", "fwdbonf_riskclass2", 1.267338843, 0.117173977, 2.021944707, 0.043182064, 1.007414265, 1.59493293,
        "adjusted", "HT", "fwdbonf_riskclass3", 1.451320224, 0.11830646, 3.148379577, 0.001641784, 1.151270487, 1.830801589,
        "adjusted", "HT", "fwdbonf_riskclass4", 1.987662079, 0.121086397, 5.673297166, 1.40E-08, 1.568774721, 2.522097434,
        "unadjusted", "HT", "fwdbonf_riskclass2", 1.409127152, 0.106304062, 3.226315771, 0.001253949, 1.144402773, 1.736172228,
        "unadjusted", "HT", "fwdbonf_riskclass3", 1.656403115, 0.107070629, 4.713229547, 2.44E-06, 1.343417696, 2.044239652,
        "unadjusted", "HT", "fwdbonf_riskclass4", 2.315762473, 0.109887253, 7.64182351, 2.14E-14, 1.868582799, 2.874988278,
        "adjusted", "HT", "fwdbonf_riskpersd", 1.321190381, 0.044740907, 6.225469066, 4.80E-10, 1.211146839, 1.443464549,
        "unadjusted", "HT", "fwdbonf_riskpersd", 1.379976023, 0.041950032, 7.677374895, 1.62E-14, 1.272187535, 1.499628884,
        "adjusted", "SBP", "fwdbonf_riskclass2", 2.770114235, 0.855433249, 3.238258785, 0.001216418, 1.093495875, 4.446732595,
        "adjusted", "SBP", "fwdbonf_riskclass3", 3.669461577, 0.859204682, 4.270765341, 2.01E-05, 1.985451345, 5.353471809,
        "adjusted", "SBP", "fwdbonf_riskclass4", 6.811279502, 0.866077722, 7.864513001, 5.21E-15, 5.11379836, 8.508760645,
        "unadjusted", "SBP", "fwdbonf_riskclass2", 3.766433566, 0.897118806, 4.198366529, 2.77E-05, 2.008113018, 5.524754115,
        "unadjusted", "SBP", "fwdbonf_riskclass3", 4.903496503, 0.897118806, 5.465827349, 5.01E-08, 3.145175955, 6.661817052,
        "unadjusted", "SBP", "fwdbonf_riskclass4", 8.447070576, 0.897432868, 9.412481845, 9.65E-21, 6.688134477, 10.20600668,
        "adjusted", "SBP", "fwdbonf_riskpersd", 2.215848154, 0.305864072, 7.244551944, 5.56E-13, 1.61636559, 2.815330719,
        "unadjusted", "SBP", "fwdbonf_riskpersd", 2.721582913, 0.318103303, 8.555657498, 1.88E-17, 2.098111895, 3.345053931)
}


compare.fr02 <- function(lmrank, eicosanoids) {
    lmrank %>%
        filter(term %in% eicosanoids, response == "SBP") %>%
        full_join(., compare.names(), by = c("term" = "term")) %>%
        mutate(mean_ci = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
               mzrt = pub.mzrt.naive(term)) %>%
        arrange(name)
}
