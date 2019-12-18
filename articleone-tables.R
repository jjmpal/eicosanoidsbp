characteristics <- function(dset, type = "fr02") {
  if (type == "fr02") {
    factors <- c("HT", "HRX", "sex", "curr_smk", "curr_diab", "plate")
    names <- list("Age, y (SD)" = "AGE",
                  "Female, N (%)" = "sex",
                  "BMI, kg/m² (SD)" = "BMI",
                  "Systolic blood pressure, mmHg (SD)" = "SBP",
                  "Diastolic blood pressure, mmHg (SD)" = "DBP",
                  "Pulse pressure, mmHg (SD)" = "PP",
                  "Mean arterial pressure, mmHg (SD)" = "MAP",
                  "Hypertension, N (%)" = "HT",
                  "Antihypertensive medication, N (%)" = "HRX",
                  "Current smoker, N (%)" = "curr_smk",
                  "Diabetes mellitus, N (%)" = "curr_diab")
  } else {
    factors <- c("HT8", "HRX8", "sex", "CURRSMK8", "curr_diab8", "plate")
    names <- list("Age, y (SD)" = "AGE8",
                  "Female, N (%)" = "sex",
                  "BMI, kg/m² (SD)" = "BMI8",
                  "Systolic blood pressure, mmHg (SD)" = "SBP8",
                  "Diastolic blood pressure, mmHg (SD)" = "DBP8",
                  "Pulse pressure, mmHg (SD)" = "PP8",
                  "Mean arterial pressure, mmHg (SD)" = "MAP8",
                  "Hypertension, N (%)" = "HT8",
                  "Antihypertensive medication, N (%)" = "HRX8",
                  "Current smoker, N (%)" = "CURRSMK8",
                  "Diabetes mellitus, N (%)" = "curr_diab8")
  }
  tableobject <- tableone::CreateTableOne(data = dset, vars = unlist(names), factorVars = factors)
  tablecsv <- print(tableobject,
                    exact = "stage",
                    quote = FALSE,
                    noSpaces = TRUE,
                    printToggle = FALSE,
                    digits = 1,
                    pDigits = 3,
                    contDigits=1)
  
  title <- "Characteristics"
  overall <- paste0("Cases, n=", dim(dset)[1])
  
  tablecsv %>%
    as.data.frame %>%
    dplyr::filter(row_number() > 1) %>%
    dplyr::mutate(!!title := names(names), !!overall := Overall) %>%
    dplyr::select(-Overall)
}

typologyformatter <- function(data, font = 12, typology) {
  flex <- flextable(data = data) %>%
    flextable::theme_booktabs()
  
  if (!missing(typology)) {
    flex <- flex %>%
      set_header_df(mapping = typology, key = "col_keys") %>%
      merge_h(part = "header")
  }
  
  flex %>%
    flextable::bold(bold = FALSE, part = "header") %>%
    flextable::bold(bold = FALSE, part = "body") %>%
    flextable::fontsize(size = font, part = "header") %>%
    flextable::fontsize(size = font, part = "body") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(align = "left", part = "header", j = 1) %>%
    flextable::align(align = "left", part = "body", j = 1)
}
