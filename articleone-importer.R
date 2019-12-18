importdata <- function(logtransform = FALSE,
                       normalize = FALSE,
                       replacenawithmin = FALSE) {
    dset <- importfinriskidata()
    
    if (replacenawithmin == TRUE) {
        dset$data.full <- dset$data.full %>% 
            mutate_at(vars(starts_with('mzid_')), list(replacewithmin))
    }

    if (logtransform == TRUE) {
        dset$data.full <- dset$data.full %>%
            mutate_at(vars(starts_with('mzid_')), list(log))
    }
    
    if (normalize == TRUE) {
        dset$data.full <- dset$data.full %>%
            mutate_at(vars(starts_with('mzid_')), list(scale))
    }
    return(dset)
}

importfinriskidata <- function(
                               fmeta = "eicdata/metabolites/fr02-biodatacoreImport-meta.rds",
                               fabund = "eicdata/metabolites/fr02-eicosanoids-MAD.rds",
                               by.y = "key",
                               by.x = "sample_run_id") {
    abund <- readRDS(fabund)
    meta <- readRDS(fmeta)

    dat <- merge(meta, abund$data, by.x = by.x, by.y = by.y, all = TRUE)

  filtered_dset <- dat %>%
      dplyr::filter(!is.na(Sample_ID),
                    !is.na(SYSTM),
                    !is.na(DIASM),
                  !is.na(BP_TREAT),
                  !is.na(sex),
                  !is.na(age_at_baseline),
                  !is.na(plate),
                  !is.na(BMI),
                  !is.na(CURR_SMOKE),
                  !is.na(DIAB),
                  !is.na(plate),
                  !is.na(KOL),
                  !is.na(HDL),
                  !is.na(CHD_AGE),
                  !is.na(PREVAL_CHD),
                  !is.na(INCIDENT_CHD),
                  !is.na(HFAIL),
                  !is.na(PREVAL_HFAIL),
                  !is.na(INCIDENT_HFAIL),
                  !is.na(STR),
                  !is.na(PREVAL_STR),
                  !is.na(INCIDENT_STR),
                  !is.na(CVD),
                  !is.na(PREVAL_CVD),
                  !is.na(INCIDENT_CVD),
                  !is.na(DEATH),
                  plate_well != "0010_0087") %>%
    dplyr::mutate(AGE = age_at_baseline,
                  SBP = SYSTM,
                  DBP = as.numeric(DIASM),
                  HRX = factor(BP_TREAT),
                  sex = factor(sex),
                  curr_smk = factor(CURR_SMOKE),
                  BMI = BMI,
                  curr_diab = factor(PREVAL_DIAB),
                  tot_chol = KOL,
                  hdl_chol = HDL,
                  chdtime = CHD_AGEDIFF,
                  chd = INCIDENT_CHD,
                  hxchd = PREVAL_CHD,
                  chftime = HFAIL_AGEDIFF,
                  chf = INCIDENT_HFAIL,
                  hxchf = PREVAL_HFAIL,
                  cvatime = STR_AGEDIFF,
                  cva = INCIDENT_STR,
                  hxcva = PREVAL_STR,
                  hardchdtime = CVD_AGEDIFF,
                  hardchd = INCIDENT_CVD,
                  hxhardchd = PREVAL_CVD,
                  dthtime = DEATH_AGEDIFF,
                  dth = DEATH,
                  hxdth = as.integer(0)) %>%
    dplyr::mutate(PP = SBP - DBP,
                  HT = factor(ifelse(SBP >= 140 | DBP >= 90 | HRX == 1, 1, 0)),
                  sex = factor(ifelse(sex == "Female", 1, 0)),
                  MAP = 2./3.*DBP + 1./3.*SBP)

  return(list(data.full = filtered_dset, 
              metabolite.names = abund$mapping,
              metabolite.mzids = abund$mapping %>% pull(mzid),
              original.abund.dim = dim(abund$data),
              original.meta.dim = dim(meta)))
}


replacewithmin <- function(...) {
    list <- c(...)
    min.value <- min(list, na.rm = TRUE)
    ifelse(is.na(list), min.value, list)
}
