remapRegularisationMzids <- function(regularisation, mapping) {
    regularisation %>%
        mutate(oterm = term) %>%
        group_by(oterm) %>%
        mutate(term = ifelse(oterm %in% names(mapping),
                             mapping[[unique(oterm)]],
                             oterm)) %>%
        ungroup %>%
        dplyr::filter(grepl("mzid_", term)) %>%
        select(term, estimate) %>%
        group_by(term) %>%
        summarise(estimate = sum(estimate)) %>%
        ungroup
}

c2l <- function(...) {
    l <- as.list(c(...))
    names(l) <- c(...)
    l
}

mygrep <- function(..., word, ignorecase = TRUE, complement = FALSE) {
    c(...)[xor(grepl(word, c(...), ignore.case = ignorecase), (complement == TRUE))]
}


pub.mzrt.naive <- function(mzid, formatter = "%.4f / %.2f") {
    sprintf(formatter,
            as.numeric(gsub("_.*", "", gsub("mzid_", "", mzid))),
            as.numeric(gsub(".*_", "", mzid)))
}

pub.mzrt <- function(mzid,
                     metabolitenames,
                     mark = "",
                     formatter = "%.4f / %.2f",
                     maxlength = 20) {
  matchfound <- mzid %in% metabolitenames$mzid
  names(mzid) <- mzid
  
  mzid[matchfound] <-
      substr(trimws(lapply(mzid[matchfound], function(id)
          metabolitenames$name[metabolitenames$mzid == id])), 1, maxlength)
  
  mzid[!matchfound] <-
    lapply(mzid[!matchfound], function(id)
      sprintf(formatter,
              as.numeric(gsub("_.*", "", gsub("mzid_", "", id))),
              as.numeric(gsub(".*_", "", id))))

  return(mzid)
}

calculateglm <- function(dset,
                         loops,
                         responses = list(),
                         binomials = c("HT", "HT8"),
                         covariates = c(),
                         maxcores = 20,
                         filter = "mzid") {
  parallel::mclapply(c2l(loops), function(loop) {
      lapply(responses, function(response) {
          fo <- sprintf("%s ~ %s", response, paste(c(loop, covariates), collapse = "+"))
          stats::glm(formula = as.formula(fo),
                     family = ifelse(response %in% binomials, stats::binomial, stats::gaussian),
                     data = dset) %>%
              broom::tidy() %>%
              dplyr::filter(grepl(filter, term)) %>%
              dplyr::mutate(conf.low = estimate - qnorm(1- 0.05/2) * std.error,
                            conf.high = estimate + qnorm(1- 0.05/2) * std.error,
                            response = response,
                            fo = fo) %>%
              { if (response %in% binomials) dplyr::mutate(., estimate = exp(estimate),
                                                           conf.low = exp(conf.low),
                                                           conf.high = exp(conf.high)) else . }
      }) %>%
      purrr::map_df(~as.data.frame(.x), .id="response")
  }, mc.cores = min(length(responses), maxcores)) %>%
      purrr::map_df(~as.data.frame(.x), .id="loop")
}

glm.binomial <- function(dset, term, covariates = c(), modelstr = "HT ~ %s", grepterm = "risk") {
  fo <- sprintf(modelstr, paste(c(term, covariates), collapse = " + "))
  glm(as.formula(fo), family = binomial(link=logit), data=dset) %>%
    broom::tidy(conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE) %>%
    filter(grepl(grepterm, term)) %>%
    mutate(mean_ci = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high), 
           p=pub.p(p.value))
}

glm.gaussian <- function(dset, term, covariates = c(), modelstr = "SBP ~ %s", grepterm = "risk") {
  fo <- sprintf(modelstr, paste(c(term, covariates), collapse = " + "))
  glm(as.formula(fo), family = stats::gaussian, data=dset) %>%
    broom::tidy(conf.int = TRUE, conf.level = 0.95) %>%
    filter(grepl(grepterm, term)) %>%
    mutate(mean_ci = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high), 
           p=pub.p(p.value)) %>%
    select(mean_ci, p)
}

#' Regularisation models
regularisationmodel <- function(dset,
                                mzids,
                                bonfp = NULL,
                                model_string = "SBP ~ %s",
                                forced = c("AGE", "sex", "BMI", "curr_smk",
                                           "curr_diab", "HRX", "plate"),
                                lambda = "lambda.1se") {
    regularisation = list()
  
    lm.lower <- stats::as.formula(sprintf(model_string,
                                          paste(c(forced), collapse='+')))
    lm.upper <- stats::as.formula(sprintf(model_string,
                                          paste(c(mzids, forced), collapse='+')))
  
                                        # Forward selection by Bonferroni p

    regularisation$fwdbonf <- stepwise.bonfp(full.model = lm.upper,
                                             initial.model = lm.lower,
                                             alpha.to.enter = bonfp,
                                             alpha.to.leave = bonfp,
                                             forced = forced,
                                             data=dset) %>%
        broom::tidy() %>%
        dplyr::mutate(conf.low = estimate - qnorm(1- 0.05/2) * std.error,
                      conf.high = estimate + qnorm(1- 0.05/2) * std.error)
  
                                        # Forward selection AIC
    regularisation$fwdaic <- step(lm(lm.lower, data=dset),
                                  trace = 1, direction = "forward",
                                  scope = list(lower = lm.lower, upper = lm.upper),
                                  steps = 1000000) %>%
        broom::tidy()
  
                                        # Lasso
    regularisation$lasso <- regularisationmodel.lasso(dset, mzids, lambda)
  
    return(regularisation)
}

#' Regularisation model LASSO
regularisationmodel.lasso <- function(dset, mzids, lambda = "lambda.1se") {
    lasso.continuous <- dset %>%
        dplyr::select(AGE, BMI, mzids)
    lasso.formula <- stats::as.formula("SBP ~ sex + curr_smk + curr_diab + HRX + plate")
    lasso.factors <- stats::model.matrix(lasso.formula, data=dset)[, -1]
    lasso.x <- cbind(as.matrix(lasso.continuous), as.matrix(lasso.factors))
    lasso.y <- dset$SBP
    lasso.penalty <- ifelse(grepl("mzid", colnames(lasso.x)), 1, 0)
    
    glmnet::cv.glmnet(lasso.x, y=lasso.y, alpha=1, penalty.factor = lasso.penalty) %>%
        stats::coef(s = lambda) %>%
        broom::tidy() %>%
        dplyr::rename(term = row, estimate = value)
}

#' Perform a stepwise linear regression using F tests of significance.
#' Based on the stepwise-function by Paul A. Rubin.
#' https://msu.edu/~rubin/code/stepwise_demo.nb.html
stepwise.bonfp <- function(full.model,
                           initial.model,
                           alpha.to.enter,
                           alpha.to.leave,
                           forced = c(),
                           data = NULL) {
    
    if (alpha.to.enter > alpha.to.leave) {
      warning("Your alpha-to-enter is greater than your alpha-to-leave, which could throw the function into an infinite loop.\n")
      return(NA)
    }
    if (is.null(data)) {
      data <- parent.frame()
    }
    
    if (is.character(full.model)) {
      fm <- as.formula(full.model)
    } else {
      fm <- as.formula(capture.output(print(full.model, showEnv = F)))
    }
    if (is.character(initial.model)) {
      im <- as.formula(initial.model)
    } else {
      im <- as.formula(capture.output(print(initial.model, showEnv = F)))
    }
    
    # Fit the full model.
    full <- lm(fm, data);
    # Sanity check: do not allow an overspecified full model.
    if (full$df.residual < 1) {
      warning("Your full model does not have enough observations to properly estimate it.\n")
      return(NA)
    }
    
    msef <- (summary(full)$sigma)^2;  # MSE of full model
    n <- length(full$residuals);  # sample size
    
    current <- lm(im, data);
    counter <- 0
    while (TRUE) {
      counter <- counter + 1
      temp <- summary(current);
      #print(temp$coefficients);
      p <- dim(temp$coefficients)[1]; # size
      mse <- (temp$sigma)^2; # MSE
      cp <- (n - p)*mse/msef - (n - 2*p);  # Mallow's cp
      fit <- sprintf("step=%i, S = %f, R-sq = %f, R-sq(adj) = %f, C-p = %f",
                     counter, temp$sigma, temp$r.squared, temp$adj.r.squared, cp);
      write(fit, file = "");
      # Try to drop a term (but only if more than one is left).
      if (p > 1) {
        d <- drop1(current, test = "F")
        d <- d[!(rownames(d) %in% forced), ]
        
        pmax <- suppressWarnings(max(d[, 6], na.rm = TRUE));
        
        if (pmax > alpha.to.leave) {
          var <- rownames(d)[d[,6] == pmax];
          
          # If an intercept is present, it will be the first name in the list.
          # There also could be ties for worst p-value.
          # Taking the second entry if there is more than one is a safe solution to both issues.
          if (length(var) > 1) {
            var <- var[2];
          }
          # Print out the variable to be dropped.
          write(paste("--- Dropping", var), file = "");
          # Modify the formulat to drop the chosen variable (by subtracting it from the current formula).
          f <- formula(current);
          f <- as.formula(paste(f[2], "~", paste(f[3], var, sep = " - ")), env = environment(f));
          # Fit the modified model and loop.
          current <- lm(f, data);
          next;
        }
      }
      # If we get here, we failed to drop a term; try adding one.
      # Note: add1 throws an error if nothing can be added (current == full), which we trap with tryCatch.
      a <- tryCatch(
        add1(current, fm, test = "F"),
        error = function(e) NULL
      );
      if (is.null(a)) {
        # There are no unused variables (or something went splat), so we bail out.
        break;
      }
      # Find the minimum p-value of any term (skipping the terms with no p-value). In case none of the remaining terms have a p-value (true of the intercept and any linearly dependent predictors), suppress warnings about an empty list. The test for a suitable candidate to drop will fail since pmin will be set to infinity.
      pmin <- suppressWarnings(min(a[, 6], na.rm = TRUE));
      if (pmin < alpha.to.enter) {
        # We have a candidate for addition to the model. Get the variable's name.
        var <- rownames(a)[a[,6] == pmin];
        # We have the same issue with ties and the presence of an intercept term, and the same solution, as above.
        if (length(var) > 1) {
          var <- var[2];
        }
        # Print the variable being added.
        write(paste("+++ Adding", var), file = "");
        # Add it to the current formula.
        f <- formula(current);
        f <- as.formula(paste(f[2], "~", paste(f[3], var, sep = " + ")), env = environment(f));
        # Fit the modified model and loop.
        current <- lm(f, data = data);
        next;
      }
      # If we get here, we failed to make any changes to the model; time to declare victory and exit.
      break;
    }
    current
  }


#' Formats p values
pub.p <- function(p) {
  p <- as.numeric(p)
  ifelse(p < 0.01, ifelse(p<0.001, "p<0.001", sprintf("%.3f", p)), sprintf("%.2f", p))
}

ctolist <- function(c) {
  l <- as.list(c)
  names(l) <- l
  l
}
  
  
bonf.adjust <- function(..., n = 1) {
    mod <- c(...) * length(c(...)) / n
    ifelse(mod < 0, 0, ifelse(mod > 1, 1, mod))
}

pull.signf <- function(df, pullterm = "term", filterterm = NULL, plimit = 0.05) {
    df %>%
        {if (!is.null(filterterm)) dplyr::filter(., response %in% filterterm) else .} %>% 
        dplyr::filter(qval < plimit) %>%
        dplyr::pull(var = pullterm) %>% 
        sort %>%
        unique
}

filter.model <- function(df, filterterm = NULL, plimit = NULL) {
    df %>%
        {if (!is.null(filterterm)) dplyr::filter(., response %in% filterterm) else .} %>%
        {if (!is.null(plimit)) dplyr::filter(., qval < plimit) else .}
}

ret.n <- function(df, model = "SBP") {
    tests <- filter.model(df, model)
    tests.signf <- pull.signf(tests)
    return(list("all" = nrow(tests),
                "sig" =  length(tests.signf),
                "insig" = nrow(tests) - length(tests.signf)))
}

comparecohorts <- function(fr02, fhs, mapping, eicosanoids) {
    compare.pairs <- mapping[eicosanoids]
    
    compare.fhs <- fhs %>%
        mutate(oterm = term) %>%
        group_by(oterm) %>%
        mutate(term = ifelse(oterm %in% names(compare.pairs),
                             compare.pairs[[unique(oterm)]],
                             oterm)) %>%
        ungroup %>%
        filter(term %in% compare.pairs, response == "SBP8") %>%
        select(term, estimate, conf.low, conf.high, p.value)  %>%
        mutate(beta_ci = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high))
    
    compare.fr02 <- fr02 %>%
        filter(term %in% names(compare.pairs), response == "SBP") %>%
        select(term, estimate, conf.low, conf.high, p.value) %>%
        mutate(beta_ci = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high))
        
    list(FR02 = compare.fr02, FHS = compare.fhs)
}


pub.lmrank <- function(...,
                       by="term",
                       mark = "",
                       models = c2l("MAP", "SBP", "DBP", "PP", "HT"),
                       arrangebyp = FALSE) {
  linrunall <- dplyr::bind_rows(...) %>%
    dplyr::filter(grepl("mzid", term), qval < 0.05) %>%
    dplyr::mutate("betase" = sprintf("%.2f±%.2f", estimate, std.error),
                  p = pub.p(p.value)) %>%
    dplyr::select(response, term, betase, p)
  
  linrunarray <- lapply(models,
                        function(model) {linrunall %>%
                            dplyr::filter(response == model) %>%
                            dplyr::select(-response)})

  linrunarray %>%
    Reduce(function(dtf1, dtf2) dplyr::full_join(dtf1, dtf2, by=by, suffix=c(".1", ".2")), .) %>%
    dplyr::arrange(term)
}
