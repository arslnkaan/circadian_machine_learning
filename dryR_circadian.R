#!/usr/bin/env Rscript

# ============================================================
# XR-seq PUP vs ADULT DIFFERENTIAL RHYTHMICITY
#
# dryR / drylm
#
# + replicate-mean XR-seq atlas
# + Pup and Adult f_24 rhythmicity tests
# + BIC model classification
# + robust dryR result conversion
# + resume/reuse existing dryR fits
#
# ============================================================
#
# INPUT
# ------------------------------------------------------------
#
# One file per sample:
#
#   adult_brain_ZT00_R1_ts_rpkm.txt
#   adult_brain_ZT00_R2_ts_rpkm.txt
#   adult_brain_ZT04_R1_ts_rpkm.txt
#   ...
#   pup_liver_ZT20_R2_ts_rpkm.txt
#
# Each file:
#
#   gene    RPKM
#   Xpa     1.234
#   Per2    2.345
#
#
# FILE NAME FORMAT
# ------------------------------------------------------------
#
# age_tissue_ZTxx_Rx_class_rpkm.txt
#
# Example:
#
# adult_brain_ZT00_R1_ts_rpkm.txt
#
#
# dryR TWO-CONDITION MODELS
# ------------------------------------------------------------
#
# Group order is explicitly:
#
#   group 1 = Pup
#   group 2 = Adult
#
# Therefore:
#
#   Model 1 = Non-rhythmic in Pup and Adult
#
#   Model 2 = Rhythmic in Pup only
#             LOSS_IN_ADULT
#
#   Model 3 = Rhythmic in Adult only
#             GAIN_IN_ADULT
#
#   Model 4 = Rhythmic in both with same rhythm
#             SHARED_RHYTHM
#
#   Model 5 = Rhythmic in both but amplitude and/or
#             phase differs
#             ALTERED_RHYTHM
#
#
# NOTE
# ------------------------------------------------------------
#
# chosen_model_BICW is NOT a p-value.
#
# It is the BIC/Schwarz weight supporting the selected
# dryR model relative to the alternative models.
#
# f_24 is additionally run separately on Pup and Adult.
# That analysis DOES provide:
#
#   pval
#   padj (BH FDR)
#
#
# MEAN ATLAS
# ------------------------------------------------------------
#
# Mean atlas:
#
#   R1 and R2 are averaged descriptively.
#
# dryR:
#
#   R1 and R2 remain separate.
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  
  library(data.table)
  library(dryR)
  library(ggplot2)
  
})


# ============================================================
# 1. COMMAND LINE
# ============================================================

args <- commandArgs(
  trailingOnly = TRUE
)


# ------------------------------------------------------------
# Usage:
#
# Rscript 02_XRseq_dryR_pup_vs_adult.R \
#     /path/to/rpkm/files \
#     ts
#
# Argument 1:
# directory containing RPKM files
#
# Argument 2:
# ts / nts / total
# ------------------------------------------------------------


INPUT_DIR <- ifelse(
  length(args) >= 1,
  args[1],
  "."
)


REPAIR_CLASS <- ifelse(
  length(args) >= 2,
  tolower(args[2]),
  "ts"
)


if (
  !REPAIR_CLASS %in%
  c(
    "ts",
    "nts",
    "total"
  )
) {
  
  stop(
    "REPAIR_CLASS must be ts, nts, or total."
  )
  
}


OUTDIR <- file.path(
  
  INPUT_DIR,
  
  paste0(
    "dryR_",
    REPAIR_CLASS,
    "_pup_vs_adult"
  )
  
)


dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. ANALYSIS SETTINGS
# ============================================================

PERIOD <- 24


# ------------------------------------------------------------
# LOW-SIGNAL FILTER
#
# Keep a gene if:
#
#   RPKM >= MIN_RPKM
#
# in at least MIN_SAMPLES_PER_AGE samples
# in EITHER Pup or Adult.
# ------------------------------------------------------------

MIN_RPKM <- 0.1

MIN_SAMPLES_PER_AGE <- 3L


# ------------------------------------------------------------
# Gaussian dryR input
#
# RPKM is transformed:
#
#   log2(RPKM + PSEUDOCOUNT)
# ------------------------------------------------------------

PSEUDOCOUNT <- 0.1


# ------------------------------------------------------------
# dryR BIC weight threshold
#
# This is a model-support threshold,
# NOT a statistical p-value threshold.
# ------------------------------------------------------------

BICW_THRESHOLD <- 0.60


# ------------------------------------------------------------
# Single-age f_24 significance
# ------------------------------------------------------------

F24_FDR_THRESHOLD <- 0.05


# ------------------------------------------------------------
# Reuse dryR_fit.rds if it already exists AND contains
# exactly the same genes as the current filtered dataset.
#
# This should allow your already-completed liver fit
# to be reused after the previous conversion error.
# ------------------------------------------------------------

REUSE_EXISTING_FIT <- TRUE


# ------------------------------------------------------------
# dryR global summary PDF can be huge for ~20k genes.
# ------------------------------------------------------------

MAKE_DRYR_SUMMARY_PDF <- FALSE


# ============================================================
# 3. CPU SETTINGS
# ============================================================

n_cores <- suppressWarnings(
  
  as.integer(
    
    Sys.getenv(
      "SLURM_CPUS_PER_TASK",
      ""
    )
    
  )
  
)


if (
  is.na(n_cores) ||
  n_cores < 1
) {
  
  detected_cores <-
    parallel::detectCores()
  
  if (
    is.na(detected_cores)
  ) {
    
    detected_cores <- 1L
    
  }
  
  
  n_cores <- min(
    4L,
    detected_cores
  )
  
}


message(
  "Using ",
  n_cores,
  " cores."
)


# ============================================================
# 4. FIND INPUT FILES
# ============================================================

message("")
message(
  "Searching for RPKM files..."
)

message(
  "Input directory: ",
  INPUT_DIR
)

message(
  "Repair class: ",
  REPAIR_CLASS
)


all_files <- list.files(
  
  INPUT_DIR,
  
  pattern =
    "_rpkm\\.txt$",
  
  full.names = TRUE
  
)


if (
  length(all_files) == 0
) {
  
  stop(
    "No *_rpkm.txt files found in: ",
    INPUT_DIR
  )
  
}


# ============================================================
# 5. PARSE METADATA FROM FILENAMES
# ============================================================

# Example:
#
# adult_brain_ZT00_R1_ts_rpkm.txt


filename_pattern <- paste0(
  
  "(?i)^",
  
  "(pup|adult)_",
  
  "(liver|kidney|brain)_",
  
  "ZT([0-9]{1,2})_",
  
  "R([0-9]+)_",
  
  "(",
  REPAIR_CLASS,
  ")",
  
  "_rpkm\\.txt$"
  
)


keep_file <- grepl(
  
  filename_pattern,
  
  basename(all_files),
  
  perl = TRUE
  
)


files <- all_files[
  keep_file
]


if (
  length(files) == 0
) {
  
  stop(
    
    "No files matched expected format.\n\n",
    
    "Example:\n",
    
    "adult_brain_ZT00_R1_",
    REPAIR_CLASS,
    "_rpkm.txt"
    
  )
  
}


proto <- data.frame(
  
  age = character(),
  
  tissue = character(),
  
  time = character(),
  
  replicate = character(),
  
  repair_class = character(),
  
  stringsAsFactors = FALSE
  
)


parsed <- strcapture(
  
  filename_pattern,
  
  basename(files),
  
  proto = proto,
  
  perl = TRUE
  
)


metadata <- as.data.table(
  parsed
)


metadata[
  ,
  file := files
]


metadata[
  ,
  sample := sub(
    "\\.txt$",
    "",
    basename(file)
  )
]


metadata[
  ,
  age := tolower(age)
]


metadata[
  ,
  tissue := tolower(tissue)
]


metadata[
  ,
  repair_class :=
    tolower(repair_class)
]


metadata[
  ,
  time := as.numeric(time)
]


metadata[
  ,
  replicate := as.integer(replicate)
]


# ------------------------------------------------------------
# Explicit Pup -> Adult ordering.
#
# IMPORTANT for dryR models 2 and 3.
# ------------------------------------------------------------

metadata[
  ,
  age_factor := factor(
    age,
    levels = c(
      "pup",
      "adult"
    )
  )
]


setorder(
  
  metadata,
  
  tissue,
  
  age_factor,
  
  time,
  
  replicate
  
)


# ============================================================
# 6. SAVE DETECTED SAMPLE TABLE
# ============================================================

fwrite(
  
  metadata,
  
  file.path(
    OUTDIR,
    "detected_samples.tsv"
  ),
  
  sep = "\t"
  
)


message("")
message(
  "Detected ",
  nrow(metadata),
  " sample files."
)


print(
  
  metadata[
    ,
    .(
      sample,
      tissue,
      age,
      time,
      replicate,
      repair_class
    )
  ]
  
)


# ============================================================
# 7. EXPERIMENTAL DESIGN QC
# ============================================================

message("")
message(
  "Checking experimental design..."
)


if (
  anyDuplicated(
    metadata$sample
  )
) {
  
  stop(
    "Duplicate sample names detected."
  )
  
}


# ============================================================
# 7A. REPLICATE QC
# ============================================================

replicate_qc <- metadata[
  ,
  .N,
  by = .(
    tissue,
    age,
    time
  )
]


fwrite(
  
  replicate_qc,
  
  file.path(
    OUTDIR,
    "replicate_QC.tsv"
  ),
  
  sep = "\t"
  
)


if (
  any(
    replicate_qc$N != 2
  )
) {
  
  print(
    replicate_qc[
      N != 2
    ]
  )
  
  
  stop(
    "Expected exactly 2 replicates ",
    "for every tissue / age / ZT."
  )
  
}


# ============================================================
# 7B. TIMEPOINT QC
# ============================================================

timepoint_qc <- metadata[
  ,
  .(
    
    n_timepoints =
      uniqueN(time),
    
    timepoints =
      paste(
        sort(
          unique(time)
        ),
        collapse = ","
      )
    
  ),
  by = .(
    tissue,
    age
  )
]


fwrite(
  
  timepoint_qc,
  
  file.path(
    OUTDIR,
    "timepoint_QC.tsv"
  ),
  
  sep = "\t"
  
)


print(
  timepoint_qc
)


if (
  any(
    timepoint_qc$n_timepoints != 6
  )
) {
  
  stop(
    "Expected exactly 6 timepoints ",
    "for every tissue / age."
  )
  
}


# ============================================================
# 7C. PUP/ADULT TIMEPOINT MATCHING
# ============================================================

for (
  tissue_name in unique(
    metadata$tissue
  )
) {
  
  pup_times <- sort(
    
    unique(
      
      metadata[
        tissue == tissue_name &
          age == "pup",
        time
      ]
      
    )
    
  )
  
  
  adult_times <- sort(
    
    unique(
      
      metadata[
        tissue == tissue_name &
          age == "adult",
        time
      ]
      
    )
    
  )
  
  
  if (
    !identical(
      pup_times,
      adult_times
    )
  ) {
    
    stop(
      "Pup and Adult timepoints differ in ",
      tissue_name,
      "."
    )
    
  }
  
}


# ============================================================
# 8. READ ONE RPKM FILE
# ============================================================

read_rpkm_file <- function(
    file,
    sample_name
) {
  
  dat <- fread(
    file
  )
  
  
  if (
    !all(
      c(
        "gene",
        "RPKM"
      ) %in%
      names(dat)
    )
  ) {
    
    stop(
      
      "\nFile does not contain gene and RPKM columns:\n",
      
      file,
      
      "\nColumns found:\n",
      
      paste(
        names(dat),
        collapse = ", "
      )
      
    )
    
  }
  
  
  dat <- dat[
    ,
    .(
      
      gene =
        as.character(gene),
      
      RPKM =
        as.numeric(RPKM)
      
    )
  ]
  
  
  if (
    anyDuplicated(
      dat$gene
    )
  ) {
    
    duplicate_genes <- unique(
      
      dat$gene[
        duplicated(
          dat$gene
        )
      ]
      
    )
    
    
    stop(
      
      "Duplicate genes found in:\n",
      
      file,
      
      "\nExamples: ",
      
      paste(
        head(
          duplicate_genes,
          10
        ),
        collapse = ", "
      )
      
    )
    
  }
  
  
  if (
    anyNA(
      dat$gene
    ) ||
    any(
      dat$gene == ""
    )
  ) {
    
    stop(
      "Missing/blank gene identifiers in: ",
      file
    )
    
  }
  
  
  if (
    anyNA(
      dat$RPKM
    )
  ) {
    
    stop(
      "NA RPKM values in: ",
      file
    )
    
  }
  
  
  if (
    any(
      !is.finite(
        dat$RPKM
      )
    )
  ) {
    
    stop(
      "Non-finite RPKM values in: ",
      file
    )
    
  }
  
  
  if (
    any(
      dat$RPKM < 0
    )
  ) {
    
    stop(
      "Negative RPKM values in: ",
      file
    )
    
  }
  
  
  dat[
    ,
    sample := sample_name
  ]
  
  
  dat
  
}


# ============================================================
# 9. READ ALL SAMPLE FILES
# ============================================================

message("")
message(
  "Reading RPKM files..."
)


long_rpkm <- rbindlist(
  
  lapply(
    
    seq_len(
      nrow(metadata)
    ),
    
    function(i) {
      
      message(
        
        "[",
        i,
        "/",
        nrow(metadata),
        "] ",
        
        basename(
          metadata$file[i]
        )
        
      )
      
      
      read_rpkm_file(
        
        file =
          metadata$file[i],
        
        sample_name =
          metadata$sample[i]
        
      )
      
    }
    
  ),
  
  use.names = TRUE
  
)


message("")
message(
  "Total gene/sample measurements: ",
  nrow(long_rpkm)
)


# ============================================================
# 10. BUILD COMPLETE GENE x SAMPLE RPKM MATRIX
# ============================================================

message("")
message(
  "Building complete RPKM matrix..."
)


rpkm_wide <- dcast(
  
  long_rpkm,
  
  gene ~ sample,
  
  value.var = "RPKM"
  
)


if (
  anyNA(
    rpkm_wide
  )
) {
  
  stop(
    
    "\nDifferent sample files contain different gene sets.\n",
    
    "Missing genes will NOT be interpreted as RPKM = 0."
    
  )
  
}


message(
  "Genes present in every sample: ",
  nrow(rpkm_wide)
)


fwrite(
  
  rpkm_wide,
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_",
      REPAIR_CLASS,
      "_RPKM_matrix.tsv"
    )
    
  ),
  
  sep = "\t"
  
)


# ============================================================
# 11. GENERATE REPLICATE-MEAN XR-seq ATLAS
# ============================================================

message("")
message(
  "===================================================="
)

message(
  "GENERATING MEAN XR-seq ATLAS"
)

message(
  "===================================================="
)


atlas_long <- merge(
  
  long_rpkm,
  
  metadata[
    ,
    .(
      sample,
      age,
      tissue,
      time,
      replicate,
      repair_class
    )
  ],
  
  by = "sample",
  
  all.x = TRUE,
  
  sort = FALSE
  
)


if (
  anyNA(
    atlas_long$age
  ) ||
  anyNA(
    atlas_long$tissue
  ) ||
  anyNA(
    atlas_long$time
  )
) {
  
  stop(
    "Metadata could not be assigned to all RPKM measurements."
  )
  
}


mean_atlas_long <- atlas_long[
  ,
  .(
    
    n_replicates =
      .N,
    
    mean_RPKM =
      mean(
        RPKM
      ),
    
    sd_RPKM =
      if (.N > 1) {
        sd(
          RPKM
        )
      } else {
        NA_real_
      },
    
    sem_RPKM =
      if (.N > 1) {
        
        sd(
          RPKM
        ) /
          sqrt(.N)
        
      } else {
        
        NA_real_
        
      }
    
  ),
  by = .(
    gene,
    tissue,
    age,
    time,
    repair_class
  )
]


if (
  any(
    mean_atlas_long$n_replicates != 2L
  )
) {
  
  stop(
    "Some atlas values were not calculated from exactly 2 replicates."
  )
  
}


setorder(
  
  mean_atlas_long,
  
  tissue,
  
  age,
  
  gene,
  
  time
  
)


fwrite(
  
  mean_atlas_long,
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_TISSUES_",
      REPAIR_CLASS,
      "_mean_atlas_long.tsv"
    )
    
  ),
  
  sep = "\t"
  
)


# ============================================================
# 11A. WIDE MEAN ATLAS
# ============================================================

mean_atlas_cast <- copy(
  mean_atlas_long
)


mean_atlas_cast[
  ,
  atlas_column := sprintf(
    "%s_%s_ZT%02d",
    age,
    tissue,
    as.integer(time)
  )
]


mean_atlas_wide <- dcast(
  
  mean_atlas_cast,
  
  gene ~ atlas_column,
  
  value.var =
    "mean_RPKM"
  
)


fwrite(
  
  mean_atlas_wide,
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_TISSUES_",
      REPAIR_CLASS,
      "_mean_RPKM_atlas.tsv"
    )
    
  ),
  
  sep = "\t"
  
)


# ============================================================
# 11B. WIDE MEAN + SD + SEM + N ATLAS
# ============================================================

mean_atlas_stats_wide <- dcast(
  
  mean_atlas_cast,
  
  gene ~ atlas_column,
  
  value.var = c(
    "mean_RPKM",
    "sd_RPKM",
    "sem_RPKM",
    "n_replicates"
  )
  
)


fwrite(
  
  mean_atlas_stats_wide,
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_TISSUES_",
      REPAIR_CLASS,
      "_mean_RPKM_atlas_with_variability.tsv"
    )
    
  ),
  
  sep = "\t"
  
)


# ============================================================
# 12. ROBUST dryR OBJECT -> data.table CONVERSION
#
# THIS FIXES THE ERROR:
#
# Error in setnames(ans, "rn", keep.rownames[1L])
#
# We DO NOT use:
#
#   keep.rownames = "gene"
#
# Instead gene IDs are inserted explicitly.
# ============================================================

to_gene_dt <- function(
    x,
    gene_ids = NULL
) {
  
  if (
    is.null(x)
  ) {
    
    stop(
      "Attempted to convert a NULL dryR object."
    )
    
  }
  
  
  if (
    is.null(gene_ids)
  ) {
    
    gene_ids <- rownames(x)
    
  }
  
  
  if (
    is.null(gene_ids)
  ) {
    
    stop(
      "Could not determine gene identifiers."
    )
    
  }
  
  
  if (
    length(gene_ids) != nrow(x)
  ) {
    
    stop(
      
      "Gene identifier length mismatch: ",
      
      length(gene_ids),
      
      " gene IDs vs ",
      
      nrow(x),
      
      " rows."
      
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Convert through ordinary data.frame first.
  #
  # This avoids the data.table keep.rownames bug.
  # ----------------------------------------------------------
  
  tmp <- as.data.frame(
    
    x,
    
    stringsAsFactors = FALSE
    
  )
  
  
  dt <- as.data.table(
    tmp
  )
  
  
  # ----------------------------------------------------------
  # Remove an existing gene column if somehow present.
  # ----------------------------------------------------------
  
  if (
    "gene" %in%
    names(dt)
  ) {
    
    dt[
      ,
      gene := NULL
    ]
    
  }
  
  
  dt[
    ,
    gene :=
      as.character(
        gene_ids
      )
  ]
  
  
  setcolorder(
    
    dt,
    
    c(
      "gene",
      setdiff(
        names(dt),
        "gene"
      )
    )
    
  )
  
  
  dt
  
}


# ============================================================
# 13. f_24 RESULT -> TABLE
# ============================================================

f24_to_table <- function(
    f24_object,
    condition_name,
    gene_ids
) {
  
  if (
    is.null(
      f24_object$parameters
    )
  ) {
    
    stop(
      "f_24 output does not contain $parameters."
    )
    
  }
  
  
  dt <- to_gene_dt(
    
    f24_object$parameters,
    
    gene_ids =
      gene_ids
    
  )
  
  
  old_names <- setdiff(
    names(dt),
    "gene"
  )
  
  
  setnames(
    
    dt,
    
    old_names,
    
    paste0(
      old_names,
      "_",
      condition_name
    )
    
  )
  
  
  dt
  
}


# ============================================================
# 14. BIC-WEIGHT MATRIX -> TABLE
# ============================================================

bicw_to_table <- function(
    x,
    prefix,
    gene_ids
) {
  
  dt <- to_gene_dt(
    
    x,
    
    gene_ids =
      gene_ids
    
  )
  
  
  model_columns <- setdiff(
    names(dt),
    "gene"
  )
  
  
  setnames(
    
    dt,
    
    model_columns,
    
    paste0(
      prefix,
      seq_along(
        model_columns
      )
    )
    
  )
  
  
  dt
  
}


# ============================================================
# 15. CIRCULAR PHASE DIFFERENCE
# ============================================================

circular_phase_difference <- function(
    adult,
    pup,
    period = 24
) {
  
  (
    (
      adult -
        pup +
        period / 2
    ) %%
      period
  ) -
    period / 2
  
}


# ============================================================
# 16. CHECK EXISTING dryR FIT
# ============================================================

existing_fit_is_compatible <- function(
    fit,
    gene_ids
) {
  
  if (
    is.null(fit$parameters)
  ) {
    
    return(FALSE)
    
  }
  
  
  fit_gene_ids <-
    rownames(
      fit$parameters
    )
  
  
  if (
    is.null(
      fit_gene_ids
    )
  ) {
    
    return(FALSE)
    
  }
  
  
  identical(
    as.character(
      fit_gene_ids
    ),
    as.character(
      gene_ids
    )
  )
  
}


# ============================================================
# 17. RUN dryR FOR ONE TISSUE
# ============================================================

run_dryr_tissue <- function(
    tissue_name
) {
  
  message("")
  message(
    "===================================================="
  )
  
  message(
    "TISSUE: ",
    toupper(
      tissue_name
    )
  )
  
  message(
    "===================================================="
  )
  
  
  # ========================================================
  # Metadata
  # ========================================================
  
  meta <- copy(
    
    metadata[
      tissue == tissue_name
    ]
    
  )
  
  
  if (
    nrow(meta) == 0
  ) {
    
    stop(
      "No samples found for ",
      tissue_name
    )
    
  }
  
  
  # --------------------------------------------------------
  # Factor ordering is deliberate:
  #
  # Pup = group 1
  # Adult = group 2
  # --------------------------------------------------------
  
  meta[
    ,
    condition := factor(
      age,
      levels = c(
        "pup",
        "adult"
      )
    )
  ]
  
  
  setorder(
    
    meta,
    
    condition,
    
    time,
    
    replicate
    
  )
  
  
  # ========================================================
  # Extract RPKM matrix
  # ========================================================
  
  sample_cols <-
    meta$sample
  
  
  y_raw <- as.matrix(
    
    rpkm_wide[
      ,
      ..sample_cols
    ]
    
  )
  
  
  rownames(
    y_raw
  ) <- rpkm_wide$gene
  
  
  storage.mode(
    y_raw
  ) <- "numeric"
  
  
  stopifnot(
    
    identical(
      colnames(y_raw),
      meta$sample
    )
    
  )
  
  
  # ========================================================
  # SIGNAL FILTER
  # ========================================================
  
  pup_samples <- meta[
    age == "pup",
    sample
  ]
  
  
  adult_samples <- meta[
    age == "adult",
    sample
  ]
  
  
  pup_signal <- rowSums(
    
    y_raw[
      ,
      pup_samples,
      drop = FALSE
    ] >= MIN_RPKM
    
  ) >= MIN_SAMPLES_PER_AGE
  
  
  adult_signal <- rowSums(
    
    y_raw[
      ,
      adult_samples,
      drop = FALSE
    ] >= MIN_RPKM
    
  ) >= MIN_SAMPLES_PER_AGE
  
  
  # --------------------------------------------------------
  # Signal required in either Pup OR Adult.
  # --------------------------------------------------------
  
  signal_keep <-
    pup_signal |
    adult_signal
  
  
  # --------------------------------------------------------
  # Remove completely invariant profiles.
  # --------------------------------------------------------
  
  row_max <- apply(
    y_raw,
    1,
    max
  )
  
  
  row_min <- apply(
    y_raw,
    1,
    min
  )
  
  
  variable_keep <-
    row_max >
    row_min
  
  
  keep <-
    signal_keep &
    variable_keep
  
  
  # ========================================================
  # Tissue output directory
  # ========================================================
  
  tissue_dir <- file.path(
    OUTDIR,
    tissue_name
  )
  
  
  dir.create(
    
    tissue_dir,
    
    recursive = TRUE,
    
    showWarnings = FALSE
    
  )
  
  
  # ========================================================
  # FILTER QC
  # ========================================================
  
  filter_qc <- data.table(
    
    gene =
      rownames(
        y_raw
      ),
    
    pup_signal =
      pup_signal,
    
    adult_signal =
      adult_signal,
    
    variable =
      variable_keep,
    
    retained =
      keep
    
  )
  
  
  fwrite(
    
    filter_qc,
    
    file.path(
      tissue_dir,
      "gene_filter_QC.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  message(
    "Genes before filtering: ",
    nrow(y_raw)
  )
  
  
  message(
    "Genes with Pup signal: ",
    sum(pup_signal)
  )
  
  
  message(
    "Genes with Adult signal: ",
    sum(adult_signal)
  )
  
  
  message(
    "Genes retained: ",
    sum(keep)
  )
  
  
  message(
    "Genes removed: ",
    sum(!keep)
  )
  
  
  y_raw <- y_raw[
    keep,
    ,
    drop = FALSE
  ]
  
  
  gene_ids <-
    rownames(
      y_raw
    )
  
  
  # ========================================================
  # SAVE TISSUE MEAN ATLAS FOR RETAINED GENES
  # ========================================================
  
  tissue_atlas <- mean_atlas_long[
    
    tissue ==
      tissue_name &
      
      gene %in%
      gene_ids
    
  ]
  
  
  setorder(
    tissue_atlas,
    gene,
    age,
    time
  )
  
  
  fwrite(
    
    tissue_atlas,
    
    file.path(
      tissue_dir,
      "mean_RPKM_atlas_genes_used_in_dryR.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # LOG2 TRANSFORMATION
  # ========================================================
  
  y <- log2(
    y_raw +
      PSEUDOCOUNT
  )
  
  
  fwrite(
    
    data.table(
      gene = rownames(y),
      y
    ),
    
    file.path(
      tissue_dir,
      "log2_RPKM_used_for_dryR.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # SAVE DRYR METADATA
  # ========================================================
  
  dryr_metadata <- meta[
    ,
    .(
      sample,
      age,
      time,
      replicate
    )
  ]
  
  
  fwrite(
    
    dryr_metadata,
    
    file.path(
      tissue_dir,
      "dryR_metadata.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # GROUP VECTOR
  #
  # IMPORTANT:
  # Pup is first factor level.
  # Adult is second factor level.
  # ========================================================
  
  group_vector <- factor(
    
    meta$age,
    
    levels = c(
      "pup",
      "adult"
    )
    
  )
  
  
  # ========================================================
  # RUN OR REUSE drylm
  # ========================================================
  
  fit_file <- file.path(
    tissue_dir,
    "dryR_fit.rds"
  )
  
  
  use_existing <- FALSE
  
  
  if (
    REUSE_EXISTING_FIT &&
    file.exists(
      fit_file
    )
  ) {
    
    message("")
    message(
      "Existing dryR fit found."
    )
    
    
    candidate_fit <- tryCatch(
      
      readRDS(
        fit_file
      ),
      
      error =
        function(e) {
          NULL
        }
      
    )
    
    
    if (
      !is.null(
        candidate_fit
      ) &&
      existing_fit_is_compatible(
        candidate_fit,
        gene_ids
      )
    ) {
      
      dry_fit <-
        candidate_fit
      
      use_existing <-
        TRUE
      
      
      message(
        "Existing fit matches current retained genes."
      )
      
      message(
        "Reusing: ",
        fit_file
      )
      
    } else {
      
      message(
        "Existing fit does not match current gene set."
      )
      
      message(
        "drylm will be rerun."
      )
      
    }
    
  }
  
  
  if (
    !use_existing
  ) {
    
    message("")
    message(
      "Running drylm..."
    )
    
    
    dry_fit <- dryR::drylm(
      
      data =
        y,
      
      group =
        group_vector,
      
      time =
        meta$time,
      
      period =
        PERIOD,
      
      sample_name =
        meta$sample,
      
      n.cores =
        n_cores
      
    )
    
    
    saveRDS(
      
      dry_fit,
      
      fit_file
      
    )
    
  }
  
  
  message("")
  message(
    "drylm fit available."
  )
  
  
  # ========================================================
  # VERIFY dryR OUTPUT STRUCTURE
  # ========================================================
  
  required_dryr_objects <- c(
    
    "parameters",
    
    "BICW_rhythm",
    
    "BICW_mean"
    
  )
  
  
  missing_dryr_objects <- setdiff(
    
    required_dryr_objects,
    
    names(
      dry_fit
    )
    
  )
  
  
  if (
    length(
      missing_dryr_objects
    ) > 0
  ) {
    
    stop(
      
      "dryR output is missing: ",
      
      paste(
        missing_dryr_objects,
        collapse = ", "
      )
      
    )
    
  }
  
  
  # ========================================================
  # DRYR PARAMETER TABLE
  #
  # ROBUST CONVERSION -- NO keep.rownames
  # ========================================================
  
  dry_params <- to_gene_dt(
    
    dry_fit$parameters,
    
    gene_ids =
      gene_ids
    
  )
  
  
  message("")
  message(
    "dryR parameter columns:"
  )
  
  message(
    paste(
      names(dry_params),
      collapse = ", "
    )
  )
  
  
  # ========================================================
  # VERIFY EXPECTED PUP/ADULT PARAMETER NAMES
  # ========================================================
  
  expected_columns <- c(
    
    "mean_pup",
    
    "amp_pup",
    
    "phase_pup",
    
    "mean_adult",
    
    "amp_adult",
    
    "phase_adult",
    
    "chosen_model",
    
    "chosen_model_BICW",
    
    "chosen_model_mean",
    
    "chosen_model_mean_BICW"
    
  )
  
  
  missing_expected <- setdiff(
    
    expected_columns,
    
    names(
      dry_params
    )
    
  )
  
  
  if (
    length(
      missing_expected
    ) > 0
  ) {
    
    stop(
      
      "\nUnexpected dryR parameter output.\n",
      
      "Missing expected columns:\n",
      
      paste(
        missing_expected,
        collapse = ", "
      ),
      
      "\n\nAvailable columns:\n",
      
      paste(
        names(
          dry_params
        ),
        collapse = ", "
      )
      
    )
  }
  
  
  
  # ========================================================
  # MODEL LABELS
  # ========================================================
  
  dry_params[
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
        paste0(
          "Model_",
          chosen_model
        )
      
    )
  ]
  
  
  # ========================================================
  # BIC CONFIDENCE
  # ========================================================
  
  dry_params[
    ,
    high_confidence_BICW :=
      chosen_model_BICW >=
      BICW_THRESHOLD
  ]
  
  
  dry_params[
    ,
    dryR_model_high_confidence :=
      fifelse(
        
        high_confidence_BICW,
        
        dryR_model,
        
        "Ambiguous"
        
      )
  ]
  
  
  # ========================================================
  # DIFFERENTIAL RHYTHMICITY
  #
  # Models:
  #
  # 2 = Loss
  # 3 = Gain
  # 5 = Altered
  # ========================================================
  
  dry_params[
    ,
    differentially_rhythmic :=
      chosen_model %in%
      c(
        2,
        3,
        5
      )
  ]
  
  
  dry_params[
    ,
    differentially_rhythmic_high_confidence :=
      differentially_rhythmic &
      high_confidence_BICW
  ]
  
  
  # ========================================================
  # ADULT - PUP EFFECTS
  # ========================================================
  
  dry_params[
    ,
    delta_mean_adult_minus_pup :=
      mean_adult -
      mean_pup
  ]
  
  
  dry_params[
    ,
    delta_amp_adult_minus_pup :=
      amp_adult -
      amp_pup
  ]
  
  
  dry_params[
    ,
    delta_relamp_adult_minus_pup :=
      relamp_adult -
      relamp_pup
  ]
  
  
  # ========================================================
  # CIRCULAR PHASE DIFFERENCE
  # ========================================================
  
  dry_params[
    ,
    delta_phase_adult_minus_pup :=
      circular_phase_difference(
        
        adult =
          phase_adult,
        
        pup =
          phase_pup,
        
        period =
          PERIOD
        
      )
  ]
  
  
  # --------------------------------------------------------
  # Phase comparison is meaningful only when both
  # conditions have rhythmic components.
  #
  # Models 4 and 5.
  # --------------------------------------------------------
  
  dry_params[
    !chosen_model %in%
      c(
        4,
        5
      ),
    delta_phase_adult_minus_pup :=
      NA_real_
  ]
  
  
  dry_params[
    ,
    tissue :=
      tissue_name
  ]
  
  
  fwrite(
    
    dry_params,
    
    file.path(
      tissue_dir,
      "dryR_parameters_and_classification.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # ALL RHYTHMIC-MODEL BIC WEIGHTS
  # ========================================================
  
  rhythm_bicw <- bicw_to_table(
    
    dry_fit$BICW_rhythm,
    
    prefix =
      "BICW_rhythm_model_",
    
    gene_ids =
      gene_ids
    
  )
  
  
  rhythm_bicw[
    ,
    tissue :=
      tissue_name
  ]
  
  
  fwrite(
    
    rhythm_bicw,
    
    file.path(
      tissue_dir,
      "dryR_all_rhythm_model_BIC_weights.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # ALL MEAN-MODEL BIC WEIGHTS
  # ========================================================
  
  mean_bicw <- bicw_to_table(
    
    dry_fit$BICW_mean,
    
    prefix =
      "BICW_mean_model_",
    
    gene_ids =
      gene_ids
    
  )
  
  
  mean_bicw[
    ,
    tissue :=
      tissue_name
  ]
  
  
  fwrite(
    
    mean_bicw,
    
    file.path(
      tissue_dir,
      "dryR_all_mean_model_BIC_weights.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # SINGLE-CONDITION f_24:
  # PUP
  # ========================================================
  
  message("")
  message(
    "Running f_24 for Pup..."
  )
  
  
  pup_idx <- which(
    meta$age ==
      "pup"
  )
  
  
  f24_pup <- dryR::f_24(
    
    data =
      y[
        ,
        pup_idx,
        drop = FALSE
      ],
    
    time =
      meta$time[
        pup_idx
      ],
    
    period =
      PERIOD,
    
    sample_name =
      meta$sample[
        pup_idx
      ]
    
  )
  
  
  pup_stats <- f24_to_table(
    
    f24_pup,
    
    condition_name =
      "pup",
    
    gene_ids =
      gene_ids
    
  )
  
  
  # ========================================================
  # SINGLE-CONDITION f_24:
  # ADULT
  # ========================================================
  
  message(
    "Running f_24 for Adult..."
  )
  
  
  adult_idx <- which(
    meta$age ==
      "adult"
  )
  
  
  f24_adult <- dryR::f_24(
    
    data =
      y[
        ,
        adult_idx,
        drop = FALSE
      ],
    
    time =
      meta$time[
        adult_idx
      ],
    
    period =
      PERIOD,
    
    sample_name =
      meta$sample[
        adult_idx
      ]
    
  )
  
  
  adult_stats <- f24_to_table(
    
    f24_adult,
    
    condition_name =
      "adult",
    
    gene_ids =
      gene_ids
    
  )
  
  
  # ========================================================
  # WITHIN-AGE FDR RHYTHMICITY
  # ========================================================
  
  if (
    !"padj_pup" %in%
    names(
      pup_stats
    )
  ) {
    
    stop(
      "f_24 Pup results do not contain padj."
    )
    
  }
  
  
  if (
    !"padj_adult" %in%
    names(
      adult_stats
    )
  ) {
    
    stop(
      "f_24 Adult results do not contain padj."
    )
    
  }
  
  
  pup_stats[
    ,
    rhythmic_FDR_pup :=
      !is.na(
        padj_pup
      ) &
      padj_pup <
      F24_FDR_THRESHOLD
  ]
  
  
  adult_stats[
    ,
    rhythmic_FDR_adult :=
      !is.na(
        padj_adult
      ) &
      padj_adult <
      F24_FDR_THRESHOLD
  ]
  
  
  fwrite(
    
    pup_stats,
    
    file.path(
      tissue_dir,
      "f24_pup_rhythmicity.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  fwrite(
    
    adult_stats,
    
    file.path(
      tissue_dir,
      "f24_adult_rhythmicity.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # BUILD MASTER dryR TABLE
  # ========================================================
  
  master <- merge(
    
    dry_params,
    
    pup_stats,
    
    by =
      "gene",
    
    all.x =
      TRUE,
    
    sort =
      FALSE
    
  )
  
  
  master <- merge(
    
    master,
    
    adult_stats,
    
    by =
      "gene",
    
    all.x =
      TRUE,
    
    sort =
      FALSE
    
  )
  
  
  master <- merge(
    
    master,
    
    rhythm_bicw[
      ,
      !"tissue"
    ],
    
    by =
      "gene",
    
    all.x =
      TRUE,
    
    sort =
      FALSE
    
  )
  
  
  master <- merge(
    
    master,
    
    mean_bicw[
      ,
      !"tissue"
    ],
    
    by =
      "gene",
    
    all.x =
      TRUE,
    
    sort =
      FALSE
    
  )
  
  
  # ========================================================
  # f_24 DEVELOPMENTAL STATE
  # ========================================================
  
  master[
    ,
    f24_developmental_state := fcase(
      
      rhythmic_FDR_pup &
        rhythmic_FDR_adult,
      "Rhythmic_both",
      
      rhythmic_FDR_pup &
        !rhythmic_FDR_adult,
      "Loss_in_adult",
      
      !rhythmic_FDR_pup &
        rhythmic_FDR_adult,
      "Gain_in_adult",
      
      !rhythmic_FDR_pup &
        !rhythmic_FDR_adult,
      "Non_rhythmic",
      
      default =
        NA_character_
      
    )
  ]
  
  
  # ========================================================
  # dryR vs f_24 AGREEMENT
  #
  # Model 5 has no direct binary equivalent.
  #
  # If both ages pass f_24 FDR and dryR says Model 5,
  # we count this as compatible.
  # ========================================================
  
  master[
    ,
    dryR_f24_agreement := fcase(
      
      chosen_model == 1 &
        f24_developmental_state ==
        "Non_rhythmic",
      TRUE,
      
      chosen_model == 2 &
        f24_developmental_state ==
        "Loss_in_adult",
      TRUE,
      
      chosen_model == 3 &
        f24_developmental_state ==
        "Gain_in_adult",
      TRUE,
      
      chosen_model == 4 &
        f24_developmental_state ==
        "Rhythmic_both",
      TRUE,
      
      chosen_model == 5 &
        f24_developmental_state ==
        "Rhythmic_both",
      TRUE,
      
      default =
        FALSE
      
    )
  ]
  
  
  # ========================================================
  # ORDER MAIN COLUMNS
  # ========================================================
  
  important_columns <- c(
    
    "gene",
    
    "tissue",
    
    "chosen_model",
    
    "dryR_model",
    
    "chosen_model_BICW",
    
    "high_confidence_BICW",
    
    "dryR_model_high_confidence",
    
    "differentially_rhythmic",
    
    "differentially_rhythmic_high_confidence",
    
    "mean_pup",
    
    "mean_adult",
    
    "delta_mean_adult_minus_pup",
    
    "amp_pup",
    
    "amp_adult",
    
    "delta_amp_adult_minus_pup",
    
    "relamp_pup",
    
    "relamp_adult",
    
    "delta_relamp_adult_minus_pup",
    
    "phase_pup",
    
    "phase_adult",
    
    "delta_phase_adult_minus_pup",
    
    "pval_pup",
    
    "padj_pup",
    
    "rhythmic_FDR_pup",
    
    "pval_adult",
    
    "padj_adult",
    
    "rhythmic_FDR_adult",
    
    "f24_developmental_state",
    
    "dryR_f24_agreement"
    
  )
  
  
  important_columns <-
    important_columns[
      important_columns %in%
        names(master)
    ]
  
  
  setcolorder(
    
    master,
    
    c(
      
      important_columns,
      
      setdiff(
        names(master),
        important_columns
      )
      
    )
    
  )
  
  
  fwrite(
    
    master,
    
    file.path(
      tissue_dir,
      "gene_level_dryR_master.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # HIGH-CONFIDENCE DIFFERENTIAL GENES
  # ========================================================
  
  high_conf_diff <- master[
    
    differentially_rhythmic_high_confidence ==
      TRUE
    
  ]
  
  
  setorder(
    
    high_conf_diff,
    
    -chosen_model_BICW
    
  )
  
  
  fwrite(
    
    high_conf_diff,
    
    file.path(
      
      tissue_dir,
      
      paste0(
        "high_confidence_differential_rhythmicity_BICW_",
        BICW_THRESHOLD,
        ".tsv"
      )
      
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # MODEL COUNTS
  # ========================================================
  
  model_counts <- master[
    ,
    .(
      
      n_genes =
        .N,
      
      n_high_confidence =
        sum(
          high_confidence_BICW,
          na.rm = TRUE
        )
      
    ),
    by = .(
      chosen_model,
      dryR_model
    )
  ]
  
  
  setorder(
    model_counts,
    chosen_model
  )
  
  
  fwrite(
    
    model_counts,
    
    file.path(
      tissue_dir,
      "dryR_model_counts.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # HIGH-CONFIDENCE MODEL COUNTS
  # ========================================================
  
  high_confidence_counts <- master[
    high_confidence_BICW ==
      TRUE,
    .N,
    by = .(
      chosen_model,
      dryR_model
    )
  ]
  
  
  setnames(
    high_confidence_counts,
    "N",
    "n_genes"
  )
  
  
  setorder(
    high_confidence_counts,
    chosen_model
  )
  
  
  fwrite(
    
    high_confidence_counts,
    
    file.path(
      tissue_dir,
      "dryR_high_confidence_model_counts.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # DRYR vs f24 CONCORDANCE TABLE
  # ========================================================
  
  concordance <- master[
    ,
    .N,
    by = .(
      dryR_model,
      f24_developmental_state,
      dryR_f24_agreement
    )
  ]
  
  
  setnames(
    concordance,
    "N",
    "n_genes"
  )
  
  
  fwrite(
    
    concordance,
    
    file.path(
      tissue_dir,
      "dryR_vs_f24_concordance.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # PLOT 1:
  # BIC WEIGHT DISTRIBUTION
  # ========================================================
  
  p_bicw <- ggplot(
    
    master,
    
    aes(
      x =
        chosen_model_BICW
    )
    
  ) +
    
    geom_histogram(
      bins = 50
    ) +
    
    geom_vline(
      
      xintercept =
        BICW_THRESHOLD,
      
      linetype =
        "dashed",
      
      linewidth =
        0.6
      
    ) +
    
    labs(
      
      title =
        paste0(
          tools::toTitleCase(
            tissue_name
          ),
          " dryR model confidence"
        ),
      
      subtitle =
        paste0(
          toupper(
            REPAIR_CLASS
          ),
          " XR-seq"
        ),
      
      x =
        "Chosen model BIC weight",
      
      y =
        "Number of genes"
      
    ) +
    
    theme_classic(
      base_size = 13
    )
  
  
  ggsave(
    
    file.path(
      tissue_dir,
      "chosen_model_BICW_distribution.pdf"
    ),
    
    p_bicw,
    
    width = 7,
    
    height = 5
    
  )
  
  
  # ========================================================
  # PLOT 2:
  # MODEL COUNTS
  # ========================================================
  
  p_models <- ggplot(
    
    model_counts,
    
    aes(
      
      x =
        factor(
          chosen_model,
          levels =
            1:5
        ),
      
      y =
        n_genes
      
    )
    
  ) +
    
    geom_col() +
    
    scale_x_discrete(
      
      labels = c(
        
        "1" =
          "Non-rhythmic",
        
        "2" =
          "Loss",
        
        "3" =
          "Gain",
        
        "4" =
          "Shared",
        
        "5" =
          "Altered"
        
      )
      
    ) +
    
    labs(
      
      title =
        paste0(
          tools::toTitleCase(
            tissue_name
          ),
          " dryR rhythmicity models"
        ),
      
      subtitle =
        paste0(
          toupper(
            REPAIR_CLASS
          ),
          " XR-seq"
        ),
      
      x =
        "dryR model",
      
      y =
        "Number of genes"
      
    ) +
    
    theme_classic(
      base_size = 13
    )
  
  
  ggsave(
    
    file.path(
      tissue_dir,
      "dryR_model_counts.pdf"
    ),
    
    p_models,
    
    width = 8,
    
    height = 5
    
  )
  
  
  # ========================================================
  # PLOT 3:
  # HIGH-CONFIDENCE MODEL COUNTS
  # ========================================================
  
  if (
    nrow(
      high_confidence_counts
    ) > 0
  ) {
    
    p_models_high <- ggplot(
      
      high_confidence_counts,
      
      aes(
        
        x =
          factor(
            chosen_model,
            levels =
              1:5
          ),
        
        y =
          n_genes
        
      )
      
    ) +
      
      geom_col() +
      
      scale_x_discrete(
        
        labels = c(
          
          "1" =
            "Non-rhythmic",
          
          "2" =
            "Loss",
          
          "3" =
            "Gain",
          
          "4" =
            "Shared",
          
          "5" =
            "Altered"
          
        )
        
      ) +
      
      labs(
        
        title =
          paste0(
            tools::toTitleCase(
              tissue_name
            ),
            " high-confidence dryR models"
          ),
        
        subtitle =
          paste0(
            "BICW >= ",
            BICW_THRESHOLD
          ),
        
        x =
          "dryR model",
        
        y =
          "Number of genes"
        
      ) +
      
      theme_classic(
        base_size = 13
      )
    
    
    ggsave(
      
      file.path(
        tissue_dir,
        "dryR_high_confidence_model_counts.pdf"
      ),
      
      p_models_high,
      
      width = 8,
      
      height = 5
      
    )
    
  }
  
  
  # ========================================================
  # PLOT 4:
  # ADULT - PUP AMPLITUDE DISTRIBUTION
  # ========================================================
  
  p_delta_amp <- ggplot(
    
    master[
      chosen_model %in%
        c(
          2,
          3,
          4,
          5
        )
    ],
    
    aes(
      x =
        delta_amp_adult_minus_pup
    )
    
  ) +
    
    geom_histogram(
      bins = 60
    ) +
    
    geom_vline(
      
      xintercept = 0,
      
      linewidth = 0.4
      
    ) +
    
    labs(
      
      title =
        paste0(
          tools::toTitleCase(
            tissue_name
          ),
          ": Adult - Pup XR-seq amplitude"
        ),
      
      x =
        "Delta amplitude (Adult - Pup)",
      
      y =
        "Number of genes"
      
    ) +
    
    theme_classic(
      base_size = 13
    )
  
  
  ggsave(
    
    file.path(
      tissue_dir,
      "delta_amplitude_distribution.pdf"
    ),
    
    p_delta_amp,
    
    width = 7,
    
    height = 5
    
  )
  
  
  # ========================================================
  # OPTIONAL BUILT-IN dryR SUMMARY
  # ========================================================
  
  if (
    MAKE_DRYR_SUMMARY_PDF
  ) {
    
    message(
      "Generating built-in dryR summary PDF..."
    )
    
    
    try(
      
      dryR::plot_models_rhythm(
        
        dry_fit,
        
        paste0(
          tissue_dir,
          "/"
        ),
        
        period =
          PERIOD
        
      ),
      
      silent =
        TRUE
      
    )
    
  }
  
  
  # ========================================================
  # RETURN
  # ========================================================
  
  list(
    
    fit =
      dry_fit,
    
    master =
      master,
    
    model_counts =
      model_counts,
    
    high_conf_diff =
      high_conf_diff,
    
    rhythm_bicw =
      rhythm_bicw,
    
    mean_bicw =
      mean_bicw,
    
    f24_pup =
      pup_stats,
    
    f24_adult =
      adult_stats
    
  )
  
}


# ============================================================
# 18. RUN AVAILABLE TISSUES
# ============================================================

tissue_order <- c(
  "liver",
  "kidney",
  "brain"
)


available_tissues <- intersect(
  
  tissue_order,
  
  unique(
    metadata$tissue
  )
  
)


if (
  length(
    available_tissues
  ) == 0
) {
  
  stop(
    "No liver, kidney, or brain files detected."
  )
  
}


message("")
message(
  "Tissues to analyze: ",
  paste(
    available_tissues,
    collapse = ", "
  )
)


results <- lapply(
  
  available_tissues,
  
  run_dryr_tissue
  
)


names(
  results
) <- available_tissues


# ============================================================
# 19. COMBINE ALL TISSUES
# ============================================================

message("")
message(
  "Combining all tissue results..."
)


all_master <- rbindlist(
  
  lapply(
    results,
    `[[`,
    "master"
  ),
  
  fill = TRUE
  
)


all_high_conf <- rbindlist(
  
  lapply(
    results,
    `[[`,
    "high_conf_diff"
  ),
  
  fill = TRUE
  
)


all_rhythm_bicw <- rbindlist(
  
  lapply(
    results,
    `[[`,
    "rhythm_bicw"
  ),
  
  fill = TRUE
  
)


all_mean_bicw <- rbindlist(
  
  lapply(
    results,
    `[[`,
    "mean_bicw"
  ),
  
  fill = TRUE
  
)


all_f24_pup <- rbindlist(
  
  lapply(
    
    names(results),
    
    function(tissue_name) {
      
      x <- copy(results[[tissue_name]]$f24_pup)
      
      
      x[
        ,
        tissue :=
          tissue_name
      ]
      
      
      x
      
    }
    
  ),
  
  fill = TRUE
  
)


all_f24_adult <- rbindlist(
  
  lapply(
    
    names(results),
    
    function(tissue_name) {
      
      x <- copy(results[[tissue_name]]$f24_adult)
      
      
      x[
        ,
        tissue :=
          tissue_name
      ]
      
      
      x
      
    }
    
  ),
  
  fill = TRUE
  
)


# ============================================================
# 20. COMBINED MODEL COUNTS
# ============================================================

all_model_counts <- rbindlist(
  
  lapply(
    
    names(results),
    
    function(tissue_name) {
      
      x <- copy(results[[tissue_name]]$model_counts)
      
      
      x[
        ,
        tissue :=
          tissue_name
      ]
      
      
      x
      
    }
    
  ),
  
  fill = TRUE
  
)


# ============================================================
# 21. WRITE COMBINED OUTPUTS
# ============================================================

fwrite(
  
  all_master,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_gene_tissue_dryR_master.tsv"
  ),
  
  sep = "\t"
  
)


fwrite(
  
  all_model_counts,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_dryR_model_counts.tsv"
  ),
  
  sep = "\t"
  
)


fwrite(
  
  all_high_conf,
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_TISSUES_high_confidence_differential_rhythmicity_BICW_",
      BICW_THRESHOLD,
      ".tsv"
    )
    
  ),
  
  sep = "\t"
  
)


fwrite(
  
  all_rhythm_bicw,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_all_rhythm_model_BIC_weights.tsv"
  ),
  
  sep = "\t"
  
)


fwrite(
  
  all_mean_bicw,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_all_mean_model_BIC_weights.tsv"
  ),
  
  sep = "\t"
  
)


fwrite(
  
  all_f24_pup,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_f24_pup_rhythmicity.tsv"
  ),
  
  sep = "\t"
  
)


fwrite(
  
  all_f24_adult,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_f24_adult_rhythmicity.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 22. GLOBAL RANKING BY dryR BIC WEIGHT
# ============================================================

global_rank <- copy(
  all_master
)


setorder(
  
  global_rank,
  
  -chosen_model_BICW
  
)


fwrite(
  
  global_rank,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_ranked_by_dryR_BICW.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 23. HIGH-CONFIDENCE DIFFERENTIAL ONLY
#
# Model 2 = loss
# Model 3 = gain
# Model 5 = altered
# ============================================================

differential_rank <- all_master[
  
  differentially_rhythmic_high_confidence ==
    TRUE
  
]


setorder(
  
  differential_rank,
  
  tissue,
  
  -chosen_model_BICW
  
)


fwrite(
  
  differential_rank,
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_TISSUES_differential_models_2_3_5_BICW_",
      BICW_THRESHOLD,
      ".tsv"
    )
    
  ),
  
  sep = "\t"
  
)


# ============================================================
# 24. CROSS-TISSUE SUMMARY
# ============================================================

cross_tissue_summary <- all_master[
  ,
  .(
    
    n_gene_tissue_rows =
      .N,
    
    n_model1_nonrhythmic =
      sum(
        chosen_model == 1,
        na.rm = TRUE
      ),
    
    n_model2_loss =
      sum(
        chosen_model == 2,
        na.rm = TRUE
      ),
    
    n_model3_gain =
      sum(
        chosen_model == 3,
        na.rm = TRUE
      ),
    
    n_model4_shared =
      sum(
        chosen_model == 4,
        na.rm = TRUE
      ),
    
    n_model5_altered =
      sum(
        chosen_model == 5,
        na.rm = TRUE
      ),
    
    n_high_confidence =
      sum(
        high_confidence_BICW,
        na.rm = TRUE
      ),
    
    n_high_conf_differential =
      sum(
        differentially_rhythmic_high_confidence,
        na.rm = TRUE
      ),
    
    n_f24_pup_FDR =
      sum(
        rhythmic_FDR_pup,
        na.rm = TRUE
      ),
    
    n_f24_adult_FDR =
      sum(
        rhythmic_FDR_adult,
        na.rm = TRUE
      ),
    
    dryR_f24_agreement_percent =
      100 *
      mean(
        dryR_f24_agreement,
        na.rm = TRUE
      )
    
  ),
  by = tissue
]


fwrite(
  
  cross_tissue_summary,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_dryR_summary.tsv"
  ),
  
  sep = "\t"
  
)


print(
  cross_tissue_summary
)


# ============================================================
# 25. MODEL DISTRIBUTION BY TISSUE
# ============================================================

combined_model_distribution <- all_master[
  ,
  .N,
  by = .(
    tissue,
    chosen_model,
    dryR_model
  )
]


combined_model_distribution[
  ,
  percent :=
    100 *
    N /
    sum(N),
  by = tissue
]


setnames(
  
  combined_model_distribution,
  
  "N",
  
  "n_genes"
  
)


fwrite(
  
  combined_model_distribution,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_dryR_model_distribution.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 26. COMBINED MODEL DISTRIBUTION PLOT
# ============================================================

p_all_models <- ggplot(
  
  combined_model_distribution,
  
  aes(
    
    x =
      tissue,
    
    y =
      percent,
    
    fill =
      dryR_model
    
  )
  
) +
  
  geom_col() +
  
  labs(
    
    title =
      paste0(
        "Pup-to-adult ",
        toupper(
          REPAIR_CLASS
        ),
        " XR-seq rhythmicity remodeling"
      ),
    
    x =
      NULL,
    
    y =
      "Genes (%)",
    
    fill =
      "dryR model"
    
  ) +
  
  theme_classic(
    base_size = 13
  )


ggsave(
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_dryR_model_distribution.pdf"
  ),
  
  p_all_models,
  
  width = 8,
  
  height = 6
  
)


# ============================================================
# 27. ANALYSIS SETTINGS
# ============================================================

settings_table <- data.table(
  
  setting = c(
    
    "repair_class",
    
    "period",
    
    "min_rpkm",
    
    "min_samples_per_age",
    
    "pseudocount",
    
    "bicw_threshold",
    
    "f24_fdr_threshold",
    
    "n_cores",
    
    "reuse_existing_fit",
    
    "make_dryR_summary_pdf"
    
  ),
  
  value = c(
    
    REPAIR_CLASS,
    
    as.character(
      PERIOD
    ),
    
    as.character(
      MIN_RPKM
    ),
    
    as.character(
      MIN_SAMPLES_PER_AGE
    ),
    
    as.character(
      PSEUDOCOUNT
    ),
    
    as.character(
      BICW_THRESHOLD
    ),
    
    as.character(
      F24_FDR_THRESHOLD
    ),
    
    as.character(
      n_cores
    ),
    
    as.character(
      REUSE_EXISTING_FIT
    ),
    
    as.character(
      MAKE_DRYR_SUMMARY_PDF
    )
    
  )
  
)


fwrite(
  
  settings_table,
  
  file.path(
    OUTDIR,
    "analysis_settings.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 28. SESSION INFO
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
# DONE
# ============================================================

message("")
message(
  "===================================================="
)

message(
  "dryR ANALYSIS COMPLETE"
)

message(
  "===================================================="
)


message("")
message(
  "Raw RPKM matrix:"
)

message(
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_",
      REPAIR_CLASS,
      "_RPKM_matrix.tsv"
    )
    
  )
  
)


message("")
message(
  "Mean XR-seq atlas:"
)

message(
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_TISSUES_",
      REPAIR_CLASS,
      "_mean_RPKM_atlas.tsv"
    )
    
  )
  
)


message("")
message(
  "Mean atlas + SD/SEM:"
)

message(
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_TISSUES_",
      REPAIR_CLASS,
      "_mean_RPKM_atlas_with_variability.tsv"
    )
    
  )
  
)


message("")
message(
  "Main dryR gene x tissue table:"
)

message(
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_gene_tissue_dryR_master.tsv"
  )
  
)


message("")
message(
  "High-confidence differential rhythmicity:"
)

message(
  
  file.path(
    
    OUTDIR,
    
    paste0(
      "ALL_TISSUES_differential_models_2_3_5_BICW_",
      BICW_THRESHOLD,
      ".tsv"
    )
    
  )
  
)


message("")
message(
  "Pup f_24 rhythmicity:"
)

message(
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_f24_pup_rhythmicity.tsv"
  )
  
)


message("")
message(
  "Adult f_24 rhythmicity:"
)

message(
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_f24_adult_rhythmicity.tsv"
  )
  
)


message("")
message(
  "Analysis complete."
)
