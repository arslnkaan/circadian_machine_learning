#!/usr/bin/env Rscript

# ============================================================
# XR-seq PUP vs ADULT
#
# INTEGRATE limorhyde2 + dryR
#
# ============================================================
#
# PURPOSE
# ------------------------------------------------------------
#
# Merge:
#
#   limorhyde2:
#       continuous developmental rhythm changes
#
#   dryR:
#       categorical rhythmic models + BIC weights
#
#
# Then:
#
#   1. Compare agreement between methods
#   2. Summarize limorhyde2 effects across dryR Models 1-5
#   3. Generate conservative consensus classifications
#   4. Subclassify dryR Model 5 using limorhyde2
#   5. Compare liver / kidney / brain
#   6. Identify shared and tissue-specific remodeling
#
#
# INPUT
# ------------------------------------------------------------
#
# Argument 1:
#   limorhyde2 result directory OR master TSV
#
# Argument 2:
#   dryR result directory OR master TSV
#
# Argument 3:
#   output directory (optional)
#
#
# Example:
#
# Rscript 03_XRseq_integrate_limorhyde2_dryR.R \
#   /path/limorhyde2_ts_pup_vs_adult \
#   /path/dryR_ts_pup_vs_adult \
#   /path/XRseq_TS_consensus
#
#
# EXPECTED FILES
# ------------------------------------------------------------
#
# limorhyde2:
#
#   ALL_TISSUES_gene_tissue_limorhyde2_master.tsv
#
# dryR:
#
#   ALL_TISSUES_gene_tissue_dryR_master.tsv
#
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

#!/usr/bin/env Rscript

# ============================================================
# XR-seq PUP vs ADULT
# INTEGRATE limorhyde2 + dryR
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  
  library(data.table)
  library(ggplot2)
  
})


# ============================================================
# 1. USER SETTINGS — RSTUDIO
# ============================================================

LIMOR_DIR <- "limorhyde2_ts_pup_vs_adult"

DRYR_DIR <- "dryR_ts_pup_vs_adult"

OUTDIR <- "limorhyde2_dryR_consensus"


# ============================================================
# 2. ANALYSIS SETTINGS
# ============================================================

BICW_THRESHOLD <- 0.60

EMPIRICAL_QUANTILE <- 0.95

MIN_MODEL4_REFERENCE <- 50L

MIN_MODEL4_FALLBACK <- 20L

TISSUE_ORDER <- c(
  "liver",
  "kidney",
  "brain"
)

MODEL_ORDER <- c(
  "Non_rhythmic",
  "Loss_in_adult",
  "Gain_in_adult",
  "Shared_rhythm",
  "Altered_rhythm"
)


# ============================================================
# 3. INPUT FILES
# ============================================================

LIMOR_FILE <- file.path(
  LIMOR_DIR,
  "ALL_TISSUES_gene_tissue_limorhyde2_master.tsv"
)

DRYR_FILE <- file.path(
  DRYR_DIR,
  "ALL_TISSUES_gene_tissue_dryR_master.tsv"
)


if (!file.exists(LIMOR_FILE)) {
  stop("Cannot find:\n", LIMOR_FILE)
}

if (!file.exists(DRYR_FILE)) {
  stop("Cannot find:\n", DRYR_FILE)
}


dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

PLOTDIR <- file.path(
  OUTDIR,
  "plots"
)

dir.create(
  PLOTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 4. READ TABLES
# ============================================================

limor <- fread(
  LIMOR_FILE
)

dryr <- fread(
  DRYR_FILE
)

# ============================================================
# 6. BASIC COLUMN QC
# ============================================================

required_keys <- c(
  "gene",
  "tissue"
)


for (
  x in required_keys
) {
  
  if (
    !x %in%
    names(limor)
  ) {
    
    stop(
      "limorhyde2 table is missing column: ",
      x
    )
    
  }
  
  
  if (
    !x %in%
    names(dryr)
  ) {
    
    stop(
      "dryR table is missing column: ",
      x
    )
    
  }
  
}


limor[
  ,
  gene := as.character(gene)
]


dryr[
  ,
  gene := as.character(gene)
]


limor[
  ,
  tissue := tolower(
    as.character(tissue)
  )
]


dryr[
  ,
  tissue := tolower(
    as.character(tissue)
  )
]


# ============================================================
# 7. CHECK DUPLICATE GENE x TISSUE ROWS
# ============================================================

limor_duplicates <- limor[
  ,
  .N,
  by = .(
    gene,
    tissue
  )
][
  N > 1
]


dryr_duplicates <- dryr[
  ,
  .N,
  by = .(
    gene,
    tissue
  )
][
  N > 1
]


if (
  nrow(limor_duplicates) > 0
) {
  
  fwrite(
    
    limor_duplicates,
    
    file.path(
      OUTDIR,
      "ERROR_limorhyde2_duplicate_gene_tissue.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  stop(
    "Duplicate gene+tissue rows detected in limorhyde2."
  )
  
}


if (
  nrow(dryr_duplicates) > 0
) {
  
  fwrite(
    
    dryr_duplicates,
    
    file.path(
      OUTDIR,
      "ERROR_dryR_duplicate_gene_tissue.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  stop(
    "Duplicate gene+tissue rows detected in dryR."
  )
  
}


# ============================================================
# 8. VERIFY REQUIRED ANALYSIS COLUMNS
# ============================================================

required_limor <- c(
  
  "diff_peak_trough_amp",
  
  "diff_peak_phase",
  
  "diff_mesor",
  
  "diff_rhy_dist",
  
  "rms_diff_rhy"
  
)


missing_limor <- setdiff(
  required_limor,
  names(limor)
)


if (
  length(missing_limor) > 0
) {
  
  stop(
    
    "\nMissing limorhyde2 columns:\n",
    
    paste(
      missing_limor,
      collapse = ", "
    ),
    
    "\n\nAvailable columns:\n",
    
    paste(
      names(limor),
      collapse = ", "
    )
    
  )
  
}


required_dryr <- c(
  
  "chosen_model",
  
  "chosen_model_BICW"
  
)


missing_dryr <- setdiff(
  required_dryr,
  names(dryr)
)


if (
  length(missing_dryr) > 0
) {
  
  stop(
    
    "\nMissing dryR columns:\n",
    
    paste(
      missing_dryr,
      collapse = ", "
    ),
    
    "\n\nAvailable columns:\n",
    
    paste(
      names(dryr),
      collapse = ", "
    )
    
  )
  
}


# ============================================================
# 9. PREFIX SOURCE-SPECIFIC COLUMNS
# ============================================================

limor_nonkey <- setdiff(
  names(limor),
  required_keys
)


dryr_nonkey <- setdiff(
  names(dryr),
  required_keys
)


setnames(
  
  limor,
  
  limor_nonkey,
  
  paste0(
    "limor__",
    limor_nonkey
  )
  
)


setnames(
  
  dryr,
  
  dryr_nonkey,
  
  paste0(
    "dryR__",
    dryr_nonkey
  )
  
)


# ============================================================
# 10. MERGE METHODS
# ============================================================

message("")
message(
  "Merging limorhyde2 and dryR..."
)


merged <- merge(
  
  limor,
  
  dryr,
  
  by = c(
    "gene",
    "tissue"
  ),
  
  all = TRUE,
  
  sort = FALSE
  
)


message(
  "Merged gene+tissue rows: ",
  nrow(merged)
)


# ============================================================
# 11. STANDARDIZED CORE VARIABLES
# ============================================================

merged[
  ,
  chosen_model :=
    as.integer(
      dryR__chosen_model
    )
]


merged[
  ,
  chosen_model_BICW :=
    as.numeric(
      dryR__chosen_model_BICW
    )
]


merged[
  ,
  limor_diff_amp :=
    as.numeric(
      limor__diff_peak_trough_amp
    )
]


merged[
  ,
  limor_diff_phase :=
    as.numeric(
      limor__diff_peak_phase
    )
]


merged[
  ,
  limor_abs_diff_phase :=
    abs(
      limor_diff_phase
    )
]


merged[
  ,
  limor_diff_mesor :=
    as.numeric(
      limor__diff_mesor
    )
]


merged[
  ,
  limor_diff_rhy_dist :=
    as.numeric(
      limor__diff_rhy_dist
    )
]


merged[
  ,
  limor_rms_diff_rhy :=
    as.numeric(
      limor__rms_diff_rhy
    )
]


# ============================================================
# 12. METHOD AVAILABILITY
# ============================================================

merged[
  ,
  has_limorhyde2 :=
    !is.na(
      limor_diff_amp
    ) |
    !is.na(
      limor_rms_diff_rhy
    )
]


merged[
  ,
  has_dryR :=
    !is.na(
      chosen_model
    )
]


merged[
  ,
  both_methods :=
    has_limorhyde2 &
    has_dryR
]


# ============================================================
# 13. STANDARD dryR MODEL LABEL
# ============================================================

merged[
  ,
  dryR_model := fcase(
    
    chosen_model == 1,
    "Non_rhythmic",
    
    chosen_model == 2,
    "Loss_in_adult",
    
    chosen_model == 3,
    "Gain_in_adult",
    
    chosen_model == 4,
    "Shared_rhythm",
    
    chosen_model == 5,
    "Altered_rhythm",
    
    default =
      NA_character_
    
  )
]


merged[
  ,
  dryR_model := factor(
    dryR_model,
    levels = MODEL_ORDER
  )
]


merged[
  ,
  dryR_high_confidence :=
    !is.na(
      chosen_model_BICW
    ) &
    chosen_model_BICW >=
    BICW_THRESHOLD
]


# ============================================================
# 14. METHOD COVERAGE SUMMARY
# ============================================================

coverage_summary <- merged[
  ,
  .(
    
    n_rows =
      .N,
    
    n_limorhyde2 =
      sum(
        has_limorhyde2
      ),
    
    n_dryR =
      sum(
        has_dryR
      ),
    
    n_both =
      sum(
        both_methods
      ),
    
    n_limor_only =
      sum(
        has_limorhyde2 &
          !has_dryR
      ),
    
    n_dryR_only =
      sum(
        has_dryR &
          !has_limorhyde2
      )
    
  ),
  by = tissue
]


fwrite(
  
  coverage_summary,
  
  file.path(
    OUTDIR,
    "01_method_coverage_summary.tsv"
  ),
  
  sep = "\t"
  
)


print(
  coverage_summary
)


# ============================================================
# 15. DRYR MODEL COUNTS
# ============================================================

dryr_model_counts <- merged[
  has_dryR == TRUE,
  .(
    
    n_genes =
      .N,
    
    n_high_confidence =
      sum(
        dryR_high_confidence
      )
    
  ),
  by = .(
    tissue,
    dryR_model
  )
]


dryr_model_counts[
  ,
  percent :=
    100 *
    n_genes /
    sum(n_genes),
  by = tissue
]


fwrite(
  
  dryr_model_counts,
  
  file.path(
    OUTDIR,
    "02_dryR_model_counts.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 16. EMPIRICAL MODEL-4 REFERENCE THRESHOLDS
# ============================================================

safe_q <- function(
    x,
    p
) {
  
  x <- x[
    is.finite(x)
  ]
  
  
  if (
    length(x) == 0
  ) {
    
    return(
      NA_real_
    )
    
  }
  
  
  as.numeric(
    
    quantile(
      x,
      probs = p,
      na.rm = TRUE,
      names = FALSE
    )
    
  )
  
}


get_model4_thresholds <- function(
    tissue_name
) {
  
  x <- merged[
    tissue == tissue_name &
      chosen_model == 4 &
      both_methods == TRUE
  ]
  
  
  x_high <- x[
    dryR_high_confidence == TRUE
  ]
  
  
  if (
    nrow(x_high) >=
    MIN_MODEL4_REFERENCE
  ) {
    
    ref <- x_high
    
    reference_source <-
      "high_confidence_model4"
    
  } else if (
    nrow(x) >=
    MIN_MODEL4_FALLBACK
  ) {
    
    ref <- x
    
    reference_source <-
      "all_model4_fallback"
    
  } else {
    
    return(
      
      data.table(
        
        tissue =
          tissue_name,
        
        reference_source =
          "insufficient_model4_reference",
        
        n_reference =
          nrow(x),
        
        amp_threshold =
          NA_real_,
        
        phase_threshold =
          NA_real_,
        
        rms_threshold =
          NA_real_,
        
        rhy_dist_threshold =
          NA_real_
        
      )
      
    )
    
  }
  
  
  data.table(
    
    tissue =
      tissue_name,
    
    reference_source =
      reference_source,
    
    n_reference =
      nrow(ref),
    
    amp_threshold =
      safe_q(
        abs(
          ref$limor_diff_amp
        ),
        EMPIRICAL_QUANTILE
      ),
    
    phase_threshold =
      safe_q(
        ref$limor_abs_diff_phase,
        EMPIRICAL_QUANTILE
      ),
    
    rms_threshold =
      safe_q(
        ref$limor_rms_diff_rhy,
        EMPIRICAL_QUANTILE
      ),
    
    rhy_dist_threshold =
      safe_q(
        ref$limor_diff_rhy_dist,
        EMPIRICAL_QUANTILE
      )
    
  )
  
}


thresholds <- rbindlist(
  
  lapply(
    
    intersect(
      TISSUE_ORDER,
      unique(
        merged$tissue
      )
    ),
    
    get_model4_thresholds
    
  )
  
)


fwrite(
  
  thresholds,
  
  file.path(
    OUTDIR,
    "03_empirical_model4_thresholds.tsv"
  ),
  
  sep = "\t"
  
)


print(
  thresholds
)


# ============================================================
# 17. ADD THRESHOLDS TO MASTER TABLE
# ============================================================

merged <- merge(
  
  merged,
  
  thresholds,
  
  by =
    "tissue",
  
  all.x =
    TRUE,
  
  sort =
    FALSE
  
)


# ============================================================
# 18. LIMORHYDE2 EFFECT FLAGS
# ============================================================

merged[
  ,
  limor_large_amp :=
    !is.na(
      amp_threshold
    ) &
    !is.na(
      limor_diff_amp
    ) &
    abs(
      limor_diff_amp
    ) >
    amp_threshold
]


merged[
  ,
  limor_large_phase :=
    !is.na(
      phase_threshold
    ) &
    !is.na(
      limor_abs_diff_phase
    ) &
    limor_abs_diff_phase >
    phase_threshold
]


merged[
  ,
  limor_large_shape :=
    !is.na(
      rms_threshold
    ) &
    !is.na(
      limor_rms_diff_rhy
    ) &
    limor_rms_diff_rhy >
    rms_threshold
]


merged[
  ,
  limor_large_rhy_dist :=
    !is.na(
      rhy_dist_threshold
    ) &
    !is.na(
      limor_diff_rhy_dist
    ) &
    limor_diff_rhy_dist >
    rhy_dist_threshold
]


# ============================================================
# 19. DIRECTIONAL AGREEMENT
#
# This applies specifically to dryR:
#
# Model 2 = Loss
# Model 3 = Gain
#
# For those classes, limorhyde2 should independently show:
#
# Loss -> Adult-Pup amplitude < 0
# Gain -> Adult-Pup amplitude > 0
# ============================================================

merged[
  ,
  directional_concordance :=
    fcase(
      
      chosen_model == 2 &
        !is.na(
          limor_diff_amp
        ) &
        limor_diff_amp < 0,
      "Concordant",
      
      chosen_model == 2 &
        !is.na(
          limor_diff_amp
        ) &
        limor_diff_amp >= 0,
      "Discordant",
      
      chosen_model == 3 &
        !is.na(
          limor_diff_amp
        ) &
        limor_diff_amp > 0,
      "Concordant",
      
      chosen_model == 3 &
        !is.na(
          limor_diff_amp
        ) &
        limor_diff_amp <= 0,
      "Discordant",
      
      default =
        NA_character_
      
    )
]


# ============================================================
# 20. CONSENSUS BROAD CLASS
# ============================================================

merged[
  ,
  consensus_broad :=
    fcase(
      
      !both_methods,
      "Incomplete_method_data",
      
      !dryR_high_confidence,
      "Ambiguous_low_BICW",
      
      chosen_model == 1,
      "Stable_nonrhythmic",
      
      chosen_model == 2 &
        limor_diff_amp < 0,
      "Loss_in_adult",
      
      chosen_model == 2 &
        limor_diff_amp >= 0,
      "Discordant_loss_direction",
      
      chosen_model == 3 &
        limor_diff_amp > 0,
      "Gain_in_adult",
      
      chosen_model == 3 &
        limor_diff_amp <= 0,
      "Discordant_gain_direction",
      
      chosen_model == 4,
      "Stable_rhythmic",
      
      chosen_model == 5,
      "Altered_rhythm",
      
      default =
        "Other"
      
    )
]


# ============================================================
# 21. SUBTYPE ALTERED RHYTHM
#
# dryR Model 5 says the rhythm changed.
#
# limorhyde2 tells us HOW it changed.
# ============================================================

merged[
  ,
  consensus_subtype :=
    consensus_broad
]


merged[
  consensus_broad ==
    "Altered_rhythm" &
    is.na(
      amp_threshold
    ),
  consensus_subtype :=
    "Altered_rhythm_unsubtyped"
]


merged[
  consensus_broad ==
    "Altered_rhythm" &
    !is.na(
      amp_threshold
    ) &
    limor_large_amp &
    limor_large_phase,
  consensus_subtype :=
    "Altered_amplitude_and_phase"
]


merged[
  consensus_broad ==
    "Altered_rhythm" &
    !is.na(
      amp_threshold
    ) &
    limor_large_amp &
    !limor_large_phase,
  consensus_subtype :=
    "Altered_amplitude"
]


merged[
  consensus_broad ==
    "Altered_rhythm" &
    !is.na(
      phase_threshold
    ) &
    !limor_large_amp &
    limor_large_phase,
  consensus_subtype :=
    "Altered_phase"
]


merged[
  consensus_broad ==
    "Altered_rhythm" &
    !limor_large_amp &
    !limor_large_phase &
    limor_large_shape,
  consensus_subtype :=
    "Altered_shape"
]


merged[
  consensus_broad ==
    "Altered_rhythm" &
    !limor_large_amp &
    !limor_large_phase &
    !limor_large_shape,
  consensus_subtype :=
    "Altered_dryR_only"
]


# ============================================================
# 22. CONSENSUS CONFIDENCE LABEL
# ============================================================

merged[
  ,
  consensus_confidence :=
    fcase(
      
      !both_methods,
      "Incomplete",
      
      !dryR_high_confidence,
      "Low_BICW",
      
      grepl(
        "^Discordant",
        consensus_broad
      ),
      "High_confidence_discordant",
      
      dryR_high_confidence,
      "High_confidence_concordant",
      
      default =
        "Other"
      
    )
]


# ============================================================
# 23. MODEL-4 LIMORHYDE2 OUTLIER QC
#
# These are dryR shared-rhythm genes whose limorhyde2
# developmental effect lies beyond our empirical reference.
# ============================================================

merged[
  ,
  model4_limorhyde2_outlier :=
    chosen_model == 4 &
    (
      limor_large_amp |
        limor_large_phase |
        limor_large_shape
    )
]


# ============================================================
# 24. WRITE LONG CONSENSUS MASTER
# ============================================================

setorder(
  merged,
  tissue,
  gene
)


fwrite(
  
  merged,
  
  file.path(
    OUTDIR,
    "XRseq_consensus_master.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 25. DRYR / LIMORHYDE2 DIRECTIONAL CONCORDANCE
# ============================================================

direction_summary <- merged[
  chosen_model %in%
    c(
      2,
      3
    ),
  .(
    
    n =
      .N,
    
    n_concordant =
      sum(
        directional_concordance ==
          "Concordant",
        na.rm = TRUE
      ),
    
    n_discordant =
      sum(
        directional_concordance ==
          "Discordant",
        na.rm = TRUE
      ),
    
    percent_concordant =
      100 *
      mean(
        directional_concordance ==
          "Concordant",
        na.rm = TRUE
      )
    
  ),
  by = .(
    tissue,
    dryR_model,
    dryR_high_confidence
  )
]


fwrite(
  
  direction_summary,
  
  file.path(
    OUTDIR,
    "04_dryR_limorhyde2_directional_concordance.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 26. LONG TABLE OF LIMORHYDE2 METRICS
# ============================================================

metric_long <- melt(
  
  merged[
    both_methods == TRUE
  ],
  
  id.vars = c(
    "gene",
    "tissue",
    "chosen_model",
    "dryR_model",
    "chosen_model_BICW",
    "dryR_high_confidence"
  ),
  
  measure.vars = c(
    
    "limor_diff_amp",
    
    "limor_abs_diff_phase",
    
    "limor_diff_mesor",
    
    "limor_rms_diff_rhy",
    
    "limor_diff_rhy_dist"
    
  ),
  
  variable.name =
    "metric",
  
  value.name =
    "value"
  
)


# ============================================================
# 27. LIMORHYDE2 EFFECTS BY dryR MODEL
# ============================================================

safe_median <- function(
    x
) {
  
  if (
    all(
      is.na(x)
    )
  ) {
    
    return(
      NA_real_
    )
    
  }
  
  median(
    x,
    na.rm = TRUE
  )
  
}


safe_quantile <- function(
    x,
    p
) {
  
  x <- x[
    !is.na(x)
  ]
  
  
  if (
    length(x) == 0
  ) {
    
    return(
      NA_real_
    )
    
  }
  
  
  as.numeric(
    quantile(
      x,
      p,
      names = FALSE
    )
  )
  
}


model_metric_summary <- metric_long[
  ,
  .(
    
    n =
      .N,
    
    n_nonmissing =
      sum(
        !is.na(value)
      ),
    
    mean =
      mean(
        value,
        na.rm = TRUE
      ),
    
    median =
      safe_median(
        value
      ),
    
    q25 =
      safe_quantile(
        value,
        0.25
      ),
    
    q75 =
      safe_quantile(
        value,
        0.75
      )
    
  ),
  by = .(
    tissue,
    dryR_model,
    metric
  )
]


fwrite(
  
  model_metric_summary,
  
  file.path(
    OUTDIR,
    "05_limorhyde2_metrics_by_dryR_model.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 28. KRUSKAL-WALLIS TESTS
#
# Question:
#
# Do limorhyde2 effect distributions differ among dryR
# Models 1-5?
# ============================================================

run_kruskal <- function(
    tissue_name,
    metric_name
) {
  
  x <- metric_long[
    tissue ==
      tissue_name &
      metric ==
      metric_name &
      !is.na(value) &
      !is.na(dryR_model)
  ]
  
  
  group_counts <- x[
    ,
    .N,
    by = dryR_model
  ][
    N >= 2
  ]
  
  
  x <- x[
    dryR_model %in%
      group_counts$dryR_model
  ]
  
  
  if (
    uniqueN(
      x$dryR_model
    ) < 2
  ) {
    
    return(
      
      data.table(
        
        tissue =
          tissue_name,
        
        metric =
          metric_name,
        
        statistic =
          NA_real_,
        
        df =
          NA_real_,
        
        p_value =
          NA_real_
        
      )
      
    )
    
  }
  
  
  fit <- kruskal.test(
    value ~ dryR_model,
    data = x
  )
  
  
  data.table(
    
    tissue =
      tissue_name,
    
    metric =
      metric_name,
    
    statistic =
      as.numeric(
        fit$statistic
      ),
    
    df =
      as.numeric(
        fit$parameter
      ),
    
    p_value =
      fit$p.value
    
  )
  
}


kruskal_results <- rbindlist(
  
  lapply(
    
    unique(
      metric_long$tissue
    ),
    
    function(tt) {
      
      rbindlist(
        
        lapply(
          
          unique(
            metric_long$metric
          ),
          
          function(mm) {
            
            run_kruskal(
              tt,
              mm
            )
            
          }
          
        )
        
      )
      
    }
    
  )
  
)


kruskal_results[
  ,
  p_adj_BH :=
    p.adjust(
      p_value,
      method = "BH"
    )
]


fwrite(
  
  kruskal_results,
  
  file.path(
    OUTDIR,
    "06_kruskal_limorhyde2_metrics_across_dryR_models.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 29. PAIRWISE WILCOXON TESTS
# ============================================================

run_pairwise <- function(
    tissue_name,
    metric_name
) {
  
  x <- metric_long[
    tissue ==
      tissue_name &
      metric ==
      metric_name &
      !is.na(value) &
      !is.na(dryR_model)
  ]
  
  
  groups <- unique(
    as.character(
      x$dryR_model
    )
  )
  
  
  groups <- groups[
    !is.na(groups)
  ]
  
  
  if (
    length(groups) < 2
  ) {
    
    return(
      NULL
    )
    
  }
  
  
  pairs <- combn(
    groups,
    2,
    simplify = FALSE
  )
  
  
  rbindlist(
    
    lapply(
      
      pairs,
      
      function(pair) {
        
        a <- x[
          as.character(
            dryR_model
          ) ==
            pair[1],
          value
        ]
        
        
        b <- x[
          as.character(
            dryR_model
          ) ==
            pair[2],
          value
        ]
        
        
        if (
          length(a) < 3 ||
          length(b) < 3
        ) {
          
          return(
            NULL
          )
          
        }
        
        
        test <- suppressWarnings(
          
          wilcox.test(
            a,
            b,
            exact = FALSE
          )
          
        )
        
        
        data.table(
          
          tissue =
            tissue_name,
          
          metric =
            metric_name,
          
          model_A =
            pair[1],
          
          model_B =
            pair[2],
          
          n_A =
            length(a),
          
          n_B =
            length(b),
          
          median_A =
            median(
              a,
              na.rm = TRUE
            ),
          
          median_B =
            median(
              b,
              na.rm = TRUE
            ),
          
          median_difference_A_minus_B =
            median(
              a,
              na.rm = TRUE
            ) -
            median(
              b,
              na.rm = TRUE
            ),
          
          p_value =
            test$p.value
          
        )
        
      }
      
    ),
    
    fill = TRUE
    
  )
  
}


pairwise_results <- rbindlist(
  
  lapply(
    
    unique(
      metric_long$tissue
    ),
    
    function(tt) {
      
      rbindlist(
        
        lapply(
          
          unique(
            metric_long$metric
          ),
          
          function(mm) {
            
            run_pairwise(
              tt,
              mm
            )
            
          }
          
        ),
        
        fill = TRUE
        
      )
      
    }
    
  ),
  
  fill = TRUE
  
)


if (
  nrow(pairwise_results) > 0
) {
  
  pairwise_results[
    ,
    p_adj_BH :=
      p.adjust(
        p_value,
        method = "BH"
      ),
    by = .(
      tissue,
      metric
    )
  ]
  
  
  fwrite(
    
    pairwise_results,
    
    file.path(
      OUTDIR,
      "07_pairwise_wilcoxon_dryR_models.tsv"
    ),
    
    sep = "\t"
    
  )
  
}


# ============================================================
# 30. CONSENSUS COUNTS
# ============================================================

consensus_counts <- merged[
  ,
  .N,
  by = .(
    tissue,
    consensus_broad
  )
]


setnames(
  consensus_counts,
  "N",
  "n_genes"
)


consensus_counts[
  ,
  percent :=
    100 *
    n_genes /
    sum(n_genes),
  by = tissue
]


fwrite(
  
  consensus_counts,
  
  file.path(
    OUTDIR,
    "08_consensus_class_counts.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 31. CONSENSUS SUBTYPE COUNTS
# ============================================================

subtype_counts <- merged[
  ,
  .N,
  by = .(
    tissue,
    consensus_subtype
  )
]


setnames(
  subtype_counts,
  "N",
  "n_genes"
)


subtype_counts[
  ,
  percent :=
    100 *
    n_genes /
    sum(n_genes),
  by = tissue
]


fwrite(
  
  subtype_counts,
  
  file.path(
    OUTDIR,
    "09_consensus_subtype_counts.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 32. CROSS-TISSUE LONG TABLE
# ============================================================

cross_long <- merged[
  ,
  .(
    
    gene,
    
    tissue,
    
    consensus_broad,
    
    consensus_subtype,
    
    chosen_model,
    
    chosen_model_BICW,
    
    dryR_high_confidence,
    
    limor_diff_amp,
    
    limor_diff_phase,
    
    limor_rms_diff_rhy,
    
    limor_diff_rhy_dist
    
  )
]


# ============================================================
# 33. CROSS-TISSUE STATE TABLE
# ============================================================

state_wide <- dcast(
  
  cross_long,
  
  gene ~ tissue,
  
  value.var =
    "consensus_broad"
  
)


subtype_wide <- dcast(
  
  cross_long,
  
  gene ~ tissue,
  
  value.var =
    "consensus_subtype"
  
)


bicw_wide <- dcast(
  
  cross_long,
  
  gene ~ tissue,
  
  value.var =
    "chosen_model_BICW"
  
)


# ============================================================
# 34. RENAME WIDE COLUMNS
# ============================================================

for (
  tt in TISSUE_ORDER
) {
  
  if (
    tt %in%
    names(state_wide)
  ) {
    
    setnames(
      
      state_wide,
      
      tt,
      
      paste0(
        "state_",
        tt
      )
      
    )
    
  }
  
  
  if (
    tt %in%
    names(subtype_wide)
  ) {
    
    setnames(
      
      subtype_wide,
      
      tt,
      
      paste0(
        "subtype_",
        tt
      )
      
    )
    
  }
  
  
  if (
    tt %in%
    names(bicw_wide)
  ) {
    
    setnames(
      
      bicw_wide,
      
      tt,
      
      paste0(
        "BICW_",
        tt
      )
      
    )
    
  }
  
}


# ============================================================
# 35. MERGE CROSS-TISSUE TABLES
# ============================================================

cross <- merge(
  
  state_wide,
  
  subtype_wide,
  
  by =
    "gene",
  
  all =
    TRUE
  
)


cross <- merge(
  
  cross,
  
  bicw_wide,
  
  by =
    "gene",
  
  all =
    TRUE
  
)


# ============================================================
# 36. ENSURE ALL THREE TISSUE COLUMNS EXIST
# ============================================================

for (
  tt in TISSUE_ORDER
) {
  
  state_col <- paste0(
    "state_",
    tt
  )
  
  
  subtype_col <- paste0(
    "subtype_",
    tt
  )
  
  
  bicw_col <- paste0(
    "BICW_",
    tt
  )
  
  
  if (
    !state_col %in%
    names(cross)
  ) {
    
    cross[
      ,
      (state_col) :=
        NA_character_
    ]
    
  }
  
  
  if (
    !subtype_col %in%
    names(cross)
  ) {
    
    cross[
      ,
      (subtype_col) :=
        NA_character_
    ]
    
  }
  
  
  if (
    !bicw_col %in%
    names(cross)
  ) {
    
    cross[
      ,
      (bicw_col) :=
        NA_real_
    ]
    
  }
  
}


# ============================================================
# 37. CROSS-TISSUE CLASSIFICATION FUNCTION
# ============================================================

classify_cross_tissue <- function(
    liver,
    kidney,
    brain
) {
  
  states <- c(
    
    liver =
      liver,
    
    kidney =
      kidney,
    
    brain =
      brain
    
  )
  
  
  if (
    any(
      is.na(states)
    )
  ) {
    
    return(
      "Incomplete_tissue_coverage"
    )
    
  }
  
  
  problematic <- grepl(
    
    "^(Ambiguous|Discordant|Incomplete|Other)",
    
    states
    
  )
  
  
  if (
    any(problematic)
  ) {
    
    return(
      "Contains_ambiguous_or_discordant"
    )
    
  }
  
  
  gains <-
    states ==
    "Gain_in_adult"
  
  
  losses <-
    states ==
    "Loss_in_adult"
  
  
  altered <-
    states ==
    "Altered_rhythm"
  
  
  stable_rhythmic <-
    states ==
    "Stable_rhythmic"
  
  
  stable_nonrhythmic <-
    states ==
    "Stable_nonrhythmic"
  
  
  remodeled <-
    gains |
    losses |
    altered
  
  
  n_remodeled <-
    sum(
      remodeled
    )
  
  
  # ----------------------------------------------------------
  # Opposite developmental directions
  # ----------------------------------------------------------
  
  if (
    any(gains) &&
    any(losses)
  ) {
    
    return(
      "Opposite_gain_loss"
    )
    
  }
  
  
  # ----------------------------------------------------------
  # All stable
  # ----------------------------------------------------------
  
  if (
    n_remodeled == 0
  ) {
    
    if (
      all(
        stable_rhythmic
      )
    ) {
      
      return(
        "Stable_rhythmic_all_3"
      )
      
    }
    
    
    if (
      all(
        stable_nonrhythmic
      )
    ) {
      
      return(
        "Stable_nonrhythmic_all_3"
      )
      
    }
    
    
    return(
      "Stable_all_3_mixed_state"
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Shared remodeling in all tissues
  # ----------------------------------------------------------
  
  if (
    all(gains)
  ) {
    
    return(
      "Shared_gain_all_3"
    )
    
  }
  
  
  if (
    all(losses)
  ) {
    
    return(
      "Shared_loss_all_3"
    )
    
  }
  
  
  if (
    all(altered)
  ) {
    
    return(
      "Shared_altered_all_3"
    )
    
  }
  
  
  if (
    n_remodeled == 3
  ) {
    
    return(
      "Remodeled_all_3_mixed"
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Remodeling in two tissues
  # ----------------------------------------------------------
  
  if (
    n_remodeled == 2
  ) {
    
    return(
      "Remodeled_two_tissues"
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Tissue-specific remodeling
  # ----------------------------------------------------------
  
  if (
    n_remodeled == 1
  ) {
    
    tissue_name <- names(
      remodeled
    )[
      which(
        remodeled
      )
    ]
    
    
    return(
      
      paste0(
        tissue_name,
        "_specific_remodeling"
      )
      
    )
    
  }
  
  
  "Other"
  
}


# ============================================================
# 38. APPLY CROSS-TISSUE CLASSIFICATION
# ============================================================

cross[
  ,
  cross_tissue_class :=
    mapply(
      
      classify_cross_tissue,
      
      state_liver,
      
      state_kidney,
      
      state_brain,
      
      USE.NAMES = FALSE
      
    )
]


# ============================================================
# 39. NUMBER OF OBSERVED / REMODELED TISSUES
# ============================================================

cross[
  ,
  n_tissues_observed :=
    rowSums(
      
      !is.na(
        
        cbind(
          state_liver,
          state_kidney,
          state_brain
        )
        
      )
      
    )
]


remodel_states <- c(
  "Gain_in_adult",
  "Loss_in_adult",
  "Altered_rhythm"
)


cross[
  ,
  n_remodeled_tissues :=
    rowSums(
      
      cbind(
        
        state_liver %in%
          remodel_states,
        
        state_kidney %in%
          remodel_states,
        
        state_brain %in%
          remodel_states
        
      )
      
    )
]


cross[
  ,
  cross_tissue_pattern :=
    paste(
      
      state_liver,
      
      state_kidney,
      
      state_brain,
      
      sep = " | "
      
    )
]


# ============================================================
# 40. CROSS-TISSUE BICW SUMMARY
# ============================================================

cross[
  ,
  min_BICW_across_observed :=
    apply(
      
      cbind(
        BICW_liver,
        BICW_kidney,
        BICW_brain
      ),
      
      1,
      
      function(x) {
        
        if (
          all(
            is.na(x)
          )
        ) {
          
          return(
            NA_real_
          )
          
        }
        
        
        min(
          x,
          na.rm = TRUE
        )
        
      }
      
    )
]


cross[
  ,
  mean_BICW_across_observed :=
    rowMeans(
      
      cbind(
        BICW_liver,
        BICW_kidney,
        BICW_brain
      ),
      
      na.rm = TRUE
      
    )
]


# ============================================================
# 41. WRITE CROSS-TISSUE MASTER
# ============================================================

setorder(
  
  cross,
  
  -n_remodeled_tissues,
  
  -min_BICW_across_observed
  
)


fwrite(
  
  cross,
  
  file.path(
    OUTDIR,
    "XRseq_cross_tissue_consensus.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 42. CROSS-TISSUE CLASS COUNTS
# ============================================================

cross_counts <- cross[
  ,
  .N,
  by =
    cross_tissue_class
]


setnames(
  cross_counts,
  "N",
  "n_genes"
)


cross_counts[
  ,
  percent :=
    100 *
    n_genes /
    sum(
      n_genes
    )
]


setorder(
  cross_counts,
  -n_genes
)


fwrite(
  
  cross_counts,
  
  file.path(
    OUTDIR,
    "10_cross_tissue_class_counts.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 43. HIGH-CONFIDENCE SHARED REMODELING
# ============================================================

shared_remodeling <- cross[
  cross_tissue_class %in%
    c(
      
      "Shared_gain_all_3",
      
      "Shared_loss_all_3",
      
      "Shared_altered_all_3",
      
      "Remodeled_all_3_mixed"
      
    )
]


setorder(
  
  shared_remodeling,
  
  -min_BICW_across_observed
  
)


fwrite(
  
  shared_remodeling,
  
  file.path(
    OUTDIR,
    "11_shared_remodeling_genes.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 44. TISSUE-SPECIFIC REMODELING
# ============================================================

tissue_specific <- cross[
  cross_tissue_class %in%
    c(
      
      "liver_specific_remodeling",
      
      "kidney_specific_remodeling",
      
      "brain_specific_remodeling"
      
    )
]


setorder(
  
  tissue_specific,
  
  cross_tissue_class,
  
  -min_BICW_across_observed
  
)


fwrite(
  
  tissue_specific,
  
  file.path(
    OUTDIR,
    "12_tissue_specific_remodeling_genes.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 45. OPPOSITE-DIRECTION GENES
# ============================================================

opposite_direction <- cross[
  cross_tissue_class ==
    "Opposite_gain_loss"
]


setorder(
  
  opposite_direction,
  
  -min_BICW_across_observed
  
)


fwrite(
  
  opposite_direction,
  
  file.path(
    OUTDIR,
    "13_opposite_gain_loss_genes.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 46. PLOT:
# limorhyde2 DELTA AMPLITUDE BY dryR MODEL
# ============================================================

plot_dt <- copy(
  metric_long
)


plot_dt[
  ,
  dryR_model :=
    factor(
      dryR_model,
      levels = MODEL_ORDER
    )
]


p_amp <- ggplot(
  
  plot_dt[
    metric ==
      "limor_diff_amp"
  ],
  
  aes(
    x =
      dryR_model,
    y =
      value
  )
  
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_boxplot(
    outlier.shape = NA
  ) +
  
  facet_wrap(
    ~ tissue,
    scales = "free_y"
  ) +
  
  labs(
    
    title =
      "limorhyde2 amplitude change across dryR models",
    
    x =
      "dryR model",
    
    y =
      "Delta peak-to-trough amplitude (Adult - Pup)"
    
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
    
  )


ggsave(
  
  file.path(
    PLOTDIR,
    "01_limor_delta_amplitude_by_dryR_model.pdf"
  ),
  
  p_amp,
  
  width = 11,
  
  height = 5
  
)


# ============================================================
# 47. PLOT:
# ABSOLUTE PHASE CHANGE BY dryR MODEL
# ============================================================

p_phase <- ggplot(
  
  plot_dt[
    metric ==
      "limor_abs_diff_phase"
  ],
  
  aes(
    x =
      dryR_model,
    y =
      value
  )
  
) +
  
  geom_boxplot(
    outlier.shape = NA
  ) +
  
  facet_wrap(
    ~ tissue
  ) +
  
  labs(
    
    title =
      "limorhyde2 phase change across dryR models",
    
    x =
      "dryR model",
    
    y =
      "Absolute phase difference (hours)"
    
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
    
  )


ggsave(
  
  file.path(
    PLOTDIR,
    "02_limor_absolute_phase_change_by_dryR_model.pdf"
  ),
  
  p_phase,
  
  width = 11,
  
  height = 5
  
)


# ============================================================
# 48. PLOT:
# RMS RHYTHMIC DIFFERENCE BY dryR MODEL
# ============================================================

p_rms <- ggplot(
  
  plot_dt[
    metric ==
      "limor_rms_diff_rhy"
  ],
  
  aes(
    x =
      dryR_model,
    y =
      value
  )
  
) +
  
  geom_boxplot(
    outlier.shape = NA
  ) +
  
  facet_wrap(
    ~ tissue,
    scales = "free_y"
  ) +
  
  labs(
    
    title =
      "limorhyde2 RMS rhythmic difference across dryR models",
    
    x =
      "dryR model",
    
    y =
      "RMS difference in mean-centered rhythmic curves"
    
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
    
  )


ggsave(
  
  file.path(
    PLOTDIR,
    "03_limor_RMS_difference_by_dryR_model.pdf"
  ),
  
  p_rms,
  
  width = 11,
  
  height = 5
  
)


# ============================================================
# 49. PLOT:
# dryR BICW vs limorhyde2 RMS DIFFERENCE
# ============================================================

p_bicw_rms <- ggplot(
  
  merged[
    both_methods == TRUE
  ],
  
  aes(
    
    x =
      chosen_model_BICW,
    
    y =
      limor_rms_diff_rhy
    
  )
  
) +
  
  geom_point(
    alpha = 0.15,
    size = 0.6
  ) +
  
  facet_grid(
    tissue ~ dryR_model,
    scales = "free_y"
  ) +
  
  geom_vline(
    
    xintercept =
      BICW_THRESHOLD,
    
    linetype =
      "dashed"
    
  ) +
  
  labs(
    
    title =
      "dryR model confidence vs limorhyde2 rhythmic difference",
    
    x =
      "dryR chosen-model BIC weight",
    
    y =
      "limorhyde2 RMS rhythmic difference"
    
  ) +
  
  theme_classic(
    base_size = 11
  )


ggsave(
  
  file.path(
    PLOTDIR,
    "04_dryR_BICW_vs_limor_RMS_difference.pdf"
  ),
  
  p_bicw_rms,
  
  width = 12,
  
  height = 8
  
)


# ============================================================
# 50. PLOT:
# CONSENSUS CLASS DISTRIBUTION
# ============================================================

p_consensus <- ggplot(
  
  consensus_counts,
  
  aes(
    
    x =
      tissue,
    
    y =
      percent,
    
    fill =
      consensus_broad
    
  )
  
) +
  
  geom_col() +
  
  labs(
    
    title =
      "Consensus developmental XR-seq rhythmicity",
    
    x =
      NULL,
    
    y =
      "Genes (%)",
    
    fill =
      "Consensus class"
    
  ) +
  
  theme_classic(
    base_size = 12
  )


ggsave(
  
  file.path(
    PLOTDIR,
    "05_consensus_class_distribution.pdf"
  ),
  
  p_consensus,
  
  width = 9,
  
  height = 6
  
)


# ============================================================
# 51. PLOT:
# CROSS-TISSUE CLASS COUNTS
# ============================================================

p_cross <- ggplot(
  
  cross_counts,
  
  aes(
    
    x =
      reorder(
        cross_tissue_class,
        n_genes
      ),
    
    y =
      n_genes
    
  )
  
) +
  
  geom_col() +
  
  coord_flip() +
  
  labs(
    
    title =
      "Cross-tissue developmental XR-seq patterns",
    
    x =
      NULL,
    
    y =
      "Number of genes"
    
  ) +
  
  theme_classic(
    base_size = 12
  )


ggsave(
  
  file.path(
    PLOTDIR,
    "06_cross_tissue_pattern_counts.pdf"
  ),
  
  p_cross,
  
  width = 8,
  
  height = 7
  
)


# ============================================================
# 52. MODEL 5 SUBTYPE PLOT
# ============================================================

model5_counts <- merged[
  consensus_broad ==
    "Altered_rhythm",
  .N,
  by = .(
    tissue,
    consensus_subtype
  )
]


if (
  nrow(
    model5_counts
  ) > 0
) {
  
  setnames(
    model5_counts,
    "N",
    "n_genes"
  )
  
  
  model5_counts[
    ,
    percent :=
      100 *
      n_genes /
      sum(
        n_genes
      ),
    by = tissue
  ]
  
  
  fwrite(
    
    model5_counts,
    
    file.path(
      OUTDIR,
      "14_model5_subtype_counts.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  p_model5 <- ggplot(
    
    model5_counts,
    
    aes(
      
      x =
        tissue,
      
      y =
        percent,
      
      fill =
        consensus_subtype
      
    )
    
  ) +
    
    geom_col() +
    
    labs(
      
      title =
        "Subtype of high-confidence dryR Model 5 genes",
      
      x =
        NULL,
      
      y =
        "Model 5 genes (%)",
      
      fill =
        "Subtype"
      
    ) +
    
    theme_classic(
      base_size = 12
    )
  
  
  ggsave(
    
    file.path(
      PLOTDIR,
      "07_model5_subtypes.pdf"
    ),
    
    p_model5,
    
    width = 9,
    
    height = 6
    
  )
  
}


# ============================================================
# 53. SAVE ANALYSIS SETTINGS
# ============================================================

settings <- data.table(
  
  setting = c(
    
    "BICW_THRESHOLD",
    
    "EMPIRICAL_QUANTILE",
    
    "MIN_MODEL4_REFERENCE",
    
    "MIN_MODEL4_FALLBACK"
    
  ),
  
  value = c(
    
    as.character(
      BICW_THRESHOLD
    ),
    
    as.character(
      EMPIRICAL_QUANTILE
    ),
    
    as.character(
      MIN_MODEL4_REFERENCE
    ),
    
    as.character(
      MIN_MODEL4_FALLBACK
    )
    
  )
  
)


fwrite(
  
  settings,
  
  file.path(
    OUTDIR,
    "analysis_settings.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 54. SESSION INFO
# ============================================================

capture.output(
  
  sessionInfo(),
  
  file =
    file.path(
      OUTDIR,
      "sessionInfo.txt"
    )
  
)


# ============================================================
# 55. FINAL SUMMARY
# ============================================================

message("")
message(
  "===================================================="
)

message(
  "INTEGRATION COMPLETE"
)

message(
  "===================================================="
)


message("")
message(
  "Main gene x tissue consensus:"
)

message(
  
  file.path(
    OUTDIR,
    "XRseq_consensus_master.tsv"
  )
  
)


message("")
message(
  "Cross-tissue consensus:"
)

message(
  
  file.path(
    OUTDIR,
    "XRseq_cross_tissue_consensus.tsv"
  )
  
)


message("")
message(
  "Shared remodeling genes:"
)

message(
  
  file.path(
    OUTDIR,
    "11_shared_remodeling_genes.tsv"
  )
  
)


message("")
message(
  "Tissue-specific remodeling genes:"
)

message(
  
  file.path(
    OUTDIR,
    "12_tissue_specific_remodeling_genes.tsv"
  )
  
)


message("")
message(
  "Opposite gain/loss genes:"
)

message(
  
  file.path(
    OUTDIR,
    "13_opposite_gain_loss_genes.tsv"
  )
  
)


message("")
message(
  "Plots:"
)

message(
  PLOTDIR
)

message("")