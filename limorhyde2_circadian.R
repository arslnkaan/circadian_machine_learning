#!/usr/bin/env Rscript

# ============================================================
# XR-seq PUP vs ADULT DIFFERENTIAL RHYTHMICITY
# + REPLICATE-MEAN XR-seq ATLAS
#
# limorhyde2
#
# ============================================================
#
# INPUT FORMAT
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
# ANALYSIS
# ------------------------------------------------------------
#
# 1. Detect all sample files.
#
# 2. Parse:
#       age
#       tissue
#       ZT
#       replicate
#       repair class
#
# 3. Perform sample/design QC.
#
# 4. Construct complete gene x sample RPKM matrix.
#
# 5. Generate RAW-RPKM replicate-mean atlas:
#
#       mean RPKM
#       SD
#       SEM
#
#    for every:
#
#       gene x tissue x age x ZT
#
#    IMPORTANT:
#    R1/R2 are averaged ONLY for this descriptive atlas.
#
# 6. Run limorhyde2 separately for:
#
#       Liver:  Pup vs Adult
#       Kidney: Pup vs Adult
#       Brain:  Pup vs Adult
#
#    limorhyde2 uses R1 and R2 independently.
#
# 7. For every gene obtain:
#
#       Pup amplitude
#       Adult amplitude
#
#       Pup peak phase
#       Adult peak phase
#
#       Pup mesor
#       Adult mesor
#
#       Adult - Pup amplitude
#       Adult - Pup phase
#       Adult - Pup mesor
#
#       diff_rhy_dist
#       rms_diff_rhy
#
#
# NOTE
# ------------------------------------------------------------
#
# This script prepares XR-seq rhythmicity features.
#
# It does NOT yet train the final machine-learning model.
#
# The resulting ALL_TISSUES gene x tissue table will later
# be joined to:
#
#   - pup/adult RNA oscillating-gene states
#   - promoter motifs
#   - TF binding
#   - chromatin features
#   - other gene-level predictors
#
# for the final cross-tissue predictive model.
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  
  library(data.table)
  library(limorhyde2)
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
# Rscript 01_XRseq_limorhyde2_pup_vs_adult.R \
#     /path/to/rpkm/files \
#     ts
#
#
# Argument 1:
# directory containing RPKM files
#
# Argument 2:
# ts / nts / total
#
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


if (!REPAIR_CLASS %in% c(
  "ts",
  "nts",
  "total"
)) {
  
  stop(
    "REPAIR_CLASS must be ts, nts, or total."
  )
  
}


OUTDIR <- file.path(
  INPUT_DIR,
  paste0(
    "limorhyde2_",
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
# Six approximately 4-hourly timepoints:
#
# ZT00
# ZT04
# ZT08
# ZT12
# ZT16
# ZT20
# ------------------------------------------------------------

N_KNOTS <- 3L


# ------------------------------------------------------------
# LOW-SIGNAL FILTER
#
# A gene must have RPKM >= MIN_RPKM in at least
# MIN_SAMPLES_PER_AGE samples in EITHER pup OR adult.
#
# We deliberately do NOT require signal in both ages.
#
# Otherwise genes that gain or lose repair activity during
# development could be removed before differential analysis.
# ------------------------------------------------------------

MIN_RPKM <- 0.1

MIN_SAMPLES_PER_AGE <- 3L


# ------------------------------------------------------------
# RPKM transformation for limorhyde2.
#
# The mean atlas remains on RAW RPKM.
#
# limorhyde2 receives:
#
#   log2(RPKM + PSEUDOCOUNT)
# ------------------------------------------------------------

PSEUDOCOUNT <- 0.1


# ------------------------------------------------------------
# Calculate RMS rhythmicity statistics.
#
# Adds:
#
#   rms_amp
#   mean_rms_amp
#   diff_rms_amp
#   rms_diff_rhy
#
# rms_diff_rhy is particularly useful because it measures
# the RMS difference between the mean-centered fitted
# Pup and Adult rhythmic curves.
# ------------------------------------------------------------

CALCULATE_RMS <- TRUE


# ============================================================
# 3. PARALLELIZATION
# ============================================================

n_cores <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "SLURM_CPUS_PER_TASK",
      ""
    )
  )
)


if (is.na(n_cores) || n_cores < 1) {
  
  detected_cores <- parallel::detectCores()
  
  if (is.na(detected_cores)) {
    detected_cores <- 1L
  }
  
  n_cores <- min(
    4L,
    detected_cores
  )
  
}


if (
  requireNamespace(
    "doParallel",
    quietly = TRUE
  )
) {
  
  doParallel::registerDoParallel(
    cores = n_cores
  )
  
  message(
    "Using ",
    n_cores,
    " cores."
  )
  
} else {
  
  message(
    "doParallel not installed; ",
    "running without registered parallel backend."
  )
  
}


# ============================================================
# 4. FIND FILES
# ============================================================

message("")
message("Searching for RPKM files...")
message("Directory: ", INPUT_DIR)
message("Repair class: ", REPAIR_CLASS)


all_files <- list.files(
  INPUT_DIR,
  pattern = "_rpkm\\.txt$",
  full.names = TRUE
)


if (length(all_files) == 0) {
  
  stop(
    "No *_rpkm.txt files found in: ",
    INPUT_DIR
  )
  
}


# ============================================================
# 5. PARSE SAMPLE INFORMATION FROM FILENAMES
# ============================================================

# Expected:
#
# adult_brain_ZT00_R1_ts_rpkm.txt
#
# Captured:
#
# 1 = age
# 2 = tissue
# 3 = ZT
# 4 = replicate
# 5 = repair class


filename_pattern <- paste0(
  "(?i)^",
  "(pup|adult)_",
  "(liver|kidney|brain)_",
  "ZT([0-9]{1,2})_",
  "R([0-9]+)_",
  "(",
  REPAIR_CLASS,
  ")_rpkm\\.txt$"
)


keep_file <- grepl(
  filename_pattern,
  basename(all_files),
  perl = TRUE
)


files <- all_files[
  keep_file
]


if (length(files) == 0) {
  
  stop(
    "No files matched expected filename format for ",
    REPAIR_CLASS,
    ".\n\nExample expected filename:\n",
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
  repair_class := tolower(repair_class)
]


metadata[
  ,
  time := as.numeric(time)
]


metadata[
  ,
  replicate := as.integer(replicate)
]


setorder(
  metadata,
  tissue,
  age,
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
message("Checking experimental design...")


# ------------------------------------------------------------
# Check duplicate sample names
# ------------------------------------------------------------

if (anyDuplicated(metadata$sample)) {
  
  stop(
    "Duplicate sample names detected after filename parsing."
  )
  
}


# ------------------------------------------------------------
# Samples per timepoint
# ------------------------------------------------------------

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
  
  message("")
  message(
    "Age/tissue/ZT combinations not having exactly 2 samples:"
  )
  
  print(
    replicate_qc[
      N != 2
    ]
  )
  
  stop(
    "Expected exactly 2 replicates ",
    "for every age/tissue/time combination."
  )
  
}


# ------------------------------------------------------------
# Number of timepoints
# ------------------------------------------------------------

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
    "for every age/tissue."
  )
  
}


# ------------------------------------------------------------
# Check that Pup and Adult have matching timepoints
# within each tissue
# ------------------------------------------------------------

for (
  tissue_name in unique(metadata$tissue)
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
  
  if (!identical(
    pup_times,
    adult_times
  )) {
    
    stop(
      "Pup and Adult timepoints differ for tissue: ",
      tissue_name
    )
    
  }
  
}


# ============================================================
# 8. READ ONE SAMPLE FILE
# ============================================================

read_rpkm_file <- function(
    file,
    sample_name
) {
  
  dat <- fread(
    file
  )
  
  
  # --------------------------------------------------------
  # Require columns:
  #
  # gene
  # RPKM
  # --------------------------------------------------------
  
  if (
    !all(
      c(
        "gene",
        "RPKM"
      ) %in% names(dat)
    )
  ) {
    
    stop(
      "\nFile does not contain columns gene and RPKM:\n",
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
  
  
  # --------------------------------------------------------
  # Gene QC
  # --------------------------------------------------------
  
  if (
    anyDuplicated(
      dat$gene
    )
  ) {
    
    duplicated_genes <- unique(
      dat$gene[
        duplicated(
          dat$gene
        )
      ]
    )
    
    stop(
      "\nDuplicate genes in:\n",
      file,
      "\nExamples:\n",
      paste(
        head(
          duplicated_genes,
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
      "Missing/blank gene identifiers detected in: ",
      file
    )
    
  }
  
  
  if (
    anyNA(
      dat$RPKM
    )
  ) {
    
    stop(
      "NA RPKM values detected in: ",
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
      "Non-finite RPKM values detected in: ",
      file
    )
    
  }
  
  
  if (
    any(
      dat$RPKM < 0
    )
  ) {
    
    stop(
      "Negative RPKM values detected in: ",
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
message("Reading RPKM files...")


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
  "Total gene/sample measurements read: ",
  nrow(long_rpkm)
)


# ============================================================
# 10. BUILD COMPLETE GENE x SAMPLE RPKM MATRIX
# ============================================================

message("")
message(
  "Building complete gene x sample RPKM matrix..."
)


rpkm_wide <- dcast(
  
  long_rpkm,
  
  gene ~ sample,
  
  value.var = "RPKM"
  
)


# ------------------------------------------------------------
# Missing values mean gene sets differ between files.
#
# Do NOT silently interpret missing genes as RPKM = 0.
# ------------------------------------------------------------

if (
  anyNA(
    rpkm_wide
  )
) {
  
  stop(
    paste0(
      "\nDifferent files contain different gene sets.\n",
      "Missing gene/sample combinations were detected.\n",
      "The script is stopping rather than treating ",
      "missing genes as RPKM = 0."
    )
  )
  
}


message(
  "Genes present in every sample: ",
  nrow(rpkm_wide)
)


# ============================================================
# 11. SAVE COMPLETE RAW RPKM MATRIX
# ============================================================

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
# 12. GENERATE REPLICATE-MEAN XR-seq ATLAS
# ============================================================

message("")
message(
  "===================================================="
)
message(
  "GENERATING REPLICATE-MEAN XR-seq ATLAS"
)
message(
  "===================================================="
)


# ------------------------------------------------------------
# Attach sample metadata to every RPKM observation.
#
# long_rpkm:
#
# gene
# RPKM
# sample
#
# becomes:
#
# gene
# RPKM
# sample
# age
# tissue
# time
# replicate
# repair_class
# ------------------------------------------------------------

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
  ) ||
  anyNA(
    atlas_long$replicate
  )
) {
  
  stop(
    "Metadata could not be assigned to all RPKM measurements."
  )
  
}


# ------------------------------------------------------------
# Replicate-mean atlas.
#
# Mean/SD/SEM calculated from RAW RPKM,
# NOT log2-transformed RPKM.
# ------------------------------------------------------------

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


setorder(
  mean_atlas_long,
  tissue,
  age,
  gene,
  time
)


# ------------------------------------------------------------
# Verify that all atlas means were constructed from
# exactly two replicates.
# ------------------------------------------------------------

atlas_replicate_qc <- mean_atlas_long[
  ,
  .(
    min_n_replicates =
      min(
        n_replicates
      ),
    
    max_n_replicates =
      max(
        n_replicates
      ),
    
    n_gene_time_groups =
      .N
  ),
  by = .(
    tissue,
    age
  )
]


fwrite(
  
  atlas_replicate_qc,
  
  file.path(
    OUTDIR,
    "mean_atlas_replicate_QC.tsv"
  ),
  
  sep = "\t"
  
)


if (
  any(
    mean_atlas_long$n_replicates != 2L
  )
) {
  
  bad_atlas <- mean_atlas_long[
    n_replicates != 2L
  ]
  
  fwrite(
    
    bad_atlas,
    
    file.path(
      OUTDIR,
      "mean_atlas_bad_replicate_counts.tsv"
    ),
    
    sep = "\t"
    
  )
  
  stop(
    "Some gene/tissue/age/ZT atlas values were not based on exactly 2 replicates."
  )
  
}


# ============================================================
# 12A. LONG-FORM MEAN ATLAS
# ============================================================

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
# 12B. CREATE ATLAS COLUMN LABEL
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


# ------------------------------------------------------------
# Desired wide-column order:
#
# liver adult ZTs
# liver pup ZTs
# kidney adult ZTs
# kidney pup ZTs
# brain adult ZTs
# brain pup ZTs
#
# determined directly from detected metadata.
# ------------------------------------------------------------

atlas_order_dt <- unique(
  
  metadata[
    ,
    .(
      tissue,
      age,
      time
    )
  ]
  
)


atlas_order_dt[
  ,
  tissue_order := match(
    tissue,
    c(
      "liver",
      "kidney",
      "brain"
    )
  )
]


atlas_order_dt[
  ,
  age_order := match(
    age,
    c(
      "adult",
      "pup"
    )
  )
]


setorder(
  atlas_order_dt,
  tissue_order,
  age_order,
  time
)


atlas_order_dt[
  ,
  atlas_column := sprintf(
    "%s_%s_ZT%02d",
    age,
    tissue,
    as.integer(time)
  )
]


atlas_column_order <- atlas_order_dt$atlas_column


# ============================================================
# 12C. WIDE MEAN-RPKM ATLAS
# ============================================================

mean_atlas_wide <- dcast(
  
  mean_atlas_cast,
  
  gene ~ atlas_column,
  
  value.var = "mean_RPKM"
  
)


existing_atlas_order <- atlas_column_order[
  atlas_column_order %in%
    names(
      mean_atlas_wide
    )
]


setcolorder(
  
  mean_atlas_wide,
  
  c(
    "gene",
    existing_atlas_order
  )
  
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
# 12D. WIDE MEAN + SD + SEM + N ATLAS
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


# ------------------------------------------------------------
# Arrange each age/tissue/ZT together:
#
# mean
# SD
# SEM
# n
# ------------------------------------------------------------

stats_column_order <- unlist(
  
  lapply(
    
    atlas_column_order,
    
    function(x) {
      
      c(
        paste0(
          "mean_RPKM_",
          x
        ),
        paste0(
          "sd_RPKM_",
          x
        ),
        paste0(
          "sem_RPKM_",
          x
        ),
        paste0(
          "n_replicates_",
          x
        )
      )
      
    }
    
  ),
  
  use.names = FALSE
  
)


stats_column_order <- stats_column_order[
  stats_column_order %in%
    names(
      mean_atlas_stats_wide
    )
]


setcolorder(
  
  mean_atlas_stats_wide,
  
  c(
    "gene",
    stats_column_order
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


message("")
message(
  "Mean atlas complete."
)

message(
  "Genes in atlas: ",
  nrow(mean_atlas_wide)
)


# ============================================================
# 13. FUNCTION TO ANALYZE ONE TISSUE
# ============================================================

run_tissue_limorhyde <- function(
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
      "No samples found for tissue: ",
      tissue_name
    )
    
  }
  
  
  # --------------------------------------------------------
  # Require Pup and Adult.
  # --------------------------------------------------------
  
  if (
    !all(
      c(
        "pup",
        "adult"
      ) %in%
      unique(
        meta$age
      )
    )
  ) {
    
    stop(
      tissue_name,
      " does not contain both Pup and Adult samples."
    )
    
  }
  
  
  # --------------------------------------------------------
  # Explicit condition ordering.
  #
  # cond1 = pup
  # cond2 = adult
  #
  # Therefore:
  #
  # diff_* = ADULT - PUP
  # --------------------------------------------------------
  
  meta[
    ,
    cond := factor(
      age,
      levels = c(
        "pup",
        "adult"
      )
    )
  ]
  
  
  setorder(
    meta,
    cond,
    time,
    replicate
  )
  
  
  # ========================================================
  # Extract tissue RPKM matrix
  # ========================================================
  
  sample_cols <- meta$sample
  
  
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
  # LOW XR-seq SIGNAL FILTER
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
  # Keep signal in either age.
  # --------------------------------------------------------
  
  signal_keep <-
    pup_signal |
    adult_signal
  
  
  # --------------------------------------------------------
  # Remove completely invariant genes.
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
    row_max > row_min
  
  
  keep <-
    signal_keep &
    variable_keep
  
  
  qc_filter <- data.table(
    
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
  
  
  tissue_dir <- file.path(
    OUTDIR,
    tissue_name
  )
  
  
  dir.create(
    tissue_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  fwrite(
    
    qc_filter,
    
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
  
  
  # ========================================================
  # SAVE MEAN ATLAS FOR GENES ENTERING THIS MODEL
  # ========================================================
  
  model_genes <- rownames(
    y_raw
  )
  
  
  tissue_mean_atlas_model_genes <- mean_atlas_long[
    
    tissue == tissue_name &
      gene %in% model_genes
    
  ]
  
  
  setorder(
    tissue_mean_atlas_model_genes,
    gene,
    age,
    time
  )
  
  
  fwrite(
    
    tissue_mean_atlas_model_genes,
    
    file.path(
      tissue_dir,
      "mean_RPKM_atlas_genes_used_in_limorhyde2.tsv"
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
  
  
  # ========================================================
  # METADATA FOR LIMORHYDE2
  # ========================================================
  
  limorhyde_meta <- meta[
    ,
    .(
      sample,
      time,
      cond
    )
  ]
  
  
  # --------------------------------------------------------
  # Metadata rows must correspond exactly to matrix columns.
  # --------------------------------------------------------
  
  stopifnot(
    
    identical(
      limorhyde_meta$sample,
      colnames(y)
    )
    
  )
  
  
  fwrite(
    
    limorhyde_meta,
    
    file.path(
      tissue_dir,
      "limorhyde2_metadata.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # SAVE LOG2 RPKM MATRIX USED BY LIMORHYDE2
  # ========================================================
  
  y_output <- data.table(
    
    gene =
      rownames(y),
    
    y
    
  )
  
  
  fwrite(
    
    y_output,
    
    file.path(
      tissue_dir,
      "log2_RPKM_used_for_limorhyde2.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # FIT LIMORHYDE2
  # ========================================================
  
  message("")
  message(
    "Fitting periodic spline model..."
  )
  
  
  fit <- limorhyde2::getModelFit(
    
    y = y,
    
    metadata =
      as.data.frame(
        limorhyde_meta
      ),
    
    period = PERIOD,
    
    nKnots = N_KNOTS,
    
    timeColname = "time",
    
    condColname = "cond",
    
    sampleColname = "sample",
    
    method = "trend"
    
  )
  
  
  # ========================================================
  # POSTERIOR SHRINKAGE
  # ========================================================
  
  message(
    "Calculating posterior estimates..."
  )
  
  
  fit <- limorhyde2::getPosteriorFit(
    fit
  )
  
  
  saveRDS(
    
    fit,
    
    file.path(
      tissue_dir,
      "limorhyde2_fit.rds"
    )
    
  )
  
  
  # ========================================================
  # RHYTHM STATISTICS
  #
  # IMPORTANT:
  #
  # Keep this original object untouched until
  # getDiffRhythmStats() has been run.
  # ========================================================
  
  message(
    "Calculating Pup and Adult rhythm parameters..."
  )
  
  
  rhythm_stats_core <- as.data.table(
    
    limorhyde2::getRhythmStats(
      
      fit,
      
      rms =
        CALCULATE_RMS
      
    )
    
  )
  
  
  # ========================================================
  # DIFFERENTIAL RHYTHMICITY
  #
  # cond1 = pup
  # cond2 = adult
  #
  # Therefore:
  #
  # diff_* = ADULT - PUP
  # ========================================================
  
  message(
    "Calculating Adult - Pup differential rhythmicity..."
  )
  
  
  diff_stats_core <- as.data.table(
    
    limorhyde2::getDiffRhythmStats(
      
      fit,
      
      rhythm_stats_core,
      
      conds = c(
        "pup",
        "adult"
      ),
      
      dopar = TRUE
      
    )
    
  )
  
  
  # ========================================================
  # ADD TISSUE ANNOTATION AFTER LIMORHYDE2 CALCULATIONS
  # ========================================================
  
  rhythm_stats <- copy(
    rhythm_stats_core
  )
  
  
  rhythm_stats[
    ,
    tissue := tissue_name
  ]
  
  
  diff_stats <- copy(
    diff_stats_core
  )
  
  
  diff_stats[
    ,
    tissue := tissue_name
  ]
  
  
  fwrite(
    
    rhythm_stats,
    
    file.path(
      tissue_dir,
      "rhythm_stats_pup_and_adult.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  fwrite(
    
    diff_stats,
    
    file.path(
      tissue_dir,
      "differential_rhythm_Adult_minus_Pup.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # WIDE PUP / ADULT RHYTHM TABLE
  # ========================================================
  
  value_columns <- c(
    
    "peak_phase",
    
    "peak_value",
    
    "trough_phase",
    
    "trough_value",
    
    "peak_trough_amp",
    
    "mesor"
    
  )
  
  
  if (
    "rms_amp" %in%
    names(
      rhythm_stats_core
    )
  ) {
    
    value_columns <- c(
      value_columns,
      "rms_amp"
    )
    
  }
  
  
  rhythm_wide <- dcast(
    
    rhythm_stats_core,
    
    feature ~ cond,
    
    value.var =
      value_columns
    
  )
  
  
  # ========================================================
  # MASTER GENE TABLE
  # ========================================================
  
  master <- merge(
    
    rhythm_wide,
    
    diff_stats_core,
    
    by = "feature",
    
    all = TRUE
    
  )
  
  
  master[
    ,
    tissue := tissue_name
  ]
  
  
  setnames(
    master,
    "feature",
    "gene"
  )
  
  
  setcolorder(
    
    master,
    
    c(
      "gene",
      "tissue",
      setdiff(
        names(master),
        c(
          "gene",
          "tissue"
        )
      )
    )
    
  )
  
  
  fwrite(
    
    master,
    
    file.path(
      tissue_dir,
      "gene_level_limorhyde2_master.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # RANK GENES BY DEVELOPMENTAL RHYTHM DIFFERENCE
  # ========================================================
  
  ranked <- copy(
    diff_stats
  )
  
  
  if (
    "rms_diff_rhy" %in%
    names(ranked)
  ) {
    
    setorder(
      ranked,
      -rms_diff_rhy
    )
    
  } else if (
    "diff_rhy_dist" %in%
    names(ranked)
  ) {
    
    setorder(
      ranked,
      -diff_rhy_dist
    )
    
  }
  
  
  fwrite(
    
    ranked,
    
    file.path(
      tissue_dir,
      "genes_ranked_by_developmental_rhythm_change.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  # ========================================================
  # PLOT:
  #
  # Adult-Pup amplitude change
  # versus
  # Adult-Pup mesor change
  # ========================================================
  
  p_amp_mesor <- ggplot(
    
    diff_stats,
    
    aes(
      x = diff_mesor,
      y = diff_peak_trough_amp
    )
    
  ) +
    
    geom_hline(
      yintercept = 0,
      linewidth = 0.35
    ) +
    
    geom_vline(
      xintercept = 0,
      linewidth = 0.35
    ) +
    
    geom_point(
      alpha = 0.3,
      size = 0.8
    ) +
    
    labs(
      
      title = paste0(
        tools::toTitleCase(
          tissue_name
        ),
        " XR-seq: Adult vs Pup"
      ),
      
      subtitle =
        paste0(
          toupper(
            REPAIR_CLASS
          ),
          " repair"
        ),
      
      x =
        "Delta mesor (Adult - Pup)",
      
      y =
        "Delta peak-to-trough amplitude (Adult - Pup)"
      
    ) +
    
    theme_classic(
      base_size = 13
    )
  
  
  ggsave(
    
    file.path(
      tissue_dir,
      "delta_amplitude_vs_delta_mesor.pdf"
    ),
    
    p_amp_mesor,
    
    width = 7,
    height = 6
    
  )
  
  
  # ========================================================
  # PLOT:
  #
  # Overall developmental rhythmicity difference
  # ========================================================
  
  if (
    "rms_diff_rhy" %in%
    names(diff_stats)
  ) {
    
    p_distance <- ggplot(
      
      diff_stats,
      
      aes(
        x = rms_diff_rhy
      )
      
    ) +
      
      geom_histogram(
        bins = 60
      ) +
      
      labs(
        
        title =
          paste0(
            tools::toTitleCase(
              tissue_name
            ),
            ": developmental change in XR-seq rhythm"
          ),
        
        subtitle =
          paste0(
            toupper(
              REPAIR_CLASS
            ),
            " repair"
          ),
        
        x =
          "RMS difference between Pup and Adult rhythmic curves",
        
        y =
          "Number of genes"
        
      ) +
      
      theme_classic(
        base_size = 13
      )
    
    
    ggsave(
      
      file.path(
        tissue_dir,
        "rms_differential_rhythmicity_distribution.pdf"
      ),
      
      p_distance,
      
      width = 7,
      height = 5
      
    )
    
  }
  
  
  # ========================================================
  # RETURN
  # ========================================================
  
  list(
    
    fit =
      fit,
    
    rhythm =
      rhythm_stats,
    
    differential =
      diff_stats,
    
    master =
      master,
    
    retained_genes =
      model_genes
    
  )
  
}


# ============================================================
# 14. RUN LIVER, KIDNEY, BRAIN
# ============================================================

tissues <- c(
  "liver",
  "kidney",
  "brain"
)


available_tissues <- intersect(
  
  tissues,
  
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
    "No liver/kidney/brain samples found."
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
  
  run_tissue_limorhyde
  
)


names(
  results
) <- available_tissues


# ============================================================
# 15. COMBINE TISSUES
# ============================================================

message("")
message(
  "Combining tissue results..."
)


all_rhythm <- rbindlist(
  
  lapply(
    results,
    `[[`,
    "rhythm"
  ),
  
  fill = TRUE
  
)


all_diff <- rbindlist(
  
  lapply(
    results,
    `[[`,
    "differential"
  ),
  
  fill = TRUE
  
)


all_master <- rbindlist(
  
  lapply(
    results,
    `[[`,
    "master"
  ),
  
  fill = TRUE
  
)


# ============================================================
# 16. WRITE COMBINED LIMORHYDE2 TABLES
# ============================================================

fwrite(
  
  all_rhythm,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_rhythm_stats.tsv"
  ),
  
  sep = "\t"
  
)


fwrite(
  
  all_diff,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_differential_rhythm_stats.tsv"
  ),
  
  sep = "\t"
  
)


fwrite(
  
  all_master,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_gene_tissue_limorhyde2_master.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 17. GLOBAL DIFFERENTIAL-RHYTHMICITY RANKING
# ============================================================

global_rank <- copy(
  all_diff
)


if (
  "rms_diff_rhy" %in%
  names(global_rank)
) {
  
  setorder(
    global_rank,
    -rms_diff_rhy
  )
  
} else if (
  "diff_rhy_dist" %in%
  names(global_rank)
) {
  
  setorder(
    global_rank,
    -diff_rhy_dist
  )
  
}


fwrite(
  
  global_rank,
  
  file.path(
    OUTDIR,
    "ALL_TISSUES_ranked_developmental_rhythm_change.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 18. ANALYSIS SUMMARY
# ============================================================

analysis_summary <- rbindlist(
  
  lapply(
    
    available_tissues,
    
    function(tissue_name) {
      
      result <-
        results[[tissue_name]]
      
      
      data.table(
        
        repair_class =
          REPAIR_CLASS,
        
        tissue =
          tissue_name,
        
        n_input_genes =
          nrow(
            rpkm_wide
          ),
        
        n_limorhyde2_genes =
          length(
            result$retained_genes
          ),
        
        n_pup_samples =
          nrow(
            metadata[
              tissue == tissue_name &
                age == "pup"
            ]
          ),
        
        n_adult_samples =
          nrow(
            metadata[
              tissue == tissue_name &
                age == "adult"
            ]
          ),
        
        n_pup_timepoints =
          uniqueN(
            metadata[
              tissue == tissue_name &
                age == "pup",
              time
            ]
          ),
        
        n_adult_timepoints =
          uniqueN(
            metadata[
              tissue == tissue_name &
                age == "adult",
              time
            ]
          )
        
      )
      
    }
    
  )
  
)


fwrite(
  
  analysis_summary,
  
  file.path(
    OUTDIR,
    "analysis_summary.tsv"
  ),
  
  sep = "\t"
  
)


print(
  analysis_summary
)


# ============================================================
# 19. SAVE ANALYSIS SETTINGS
# ============================================================

settings_table <- data.table(
  
  setting = c(
    
    "repair_class",
    
    "period",
    
    "n_knots",
    
    "min_rpkm",
    
    "min_samples_per_age",
    
    "pseudocount",
    
    "calculate_rms",
    
    "limorhyde2_method"
    
  ),
  
  value = c(
    
    REPAIR_CLASS,
    
    as.character(
      PERIOD
    ),
    
    as.character(
      N_KNOTS
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
      CALCULATE_RMS
    ),
    
    "trend"
    
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
# 20. SESSION INFO
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
  "ANALYSIS COMPLETE"
)
message(
  "===================================================="
)


message("")
message(
  "RAW RPKM matrix:"
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
  "Mean atlas (long):"
)

message(
  file.path(
    OUTDIR,
    paste0(
      "ALL_TISSUES_",
      REPAIR_CLASS,
      "_mean_atlas_long.tsv"
    )
  )
)


message("")
message(
  "Mean RPKM atlas (wide):"
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
  "Mean + SD + SEM atlas:"
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
  "Main limorhyde2 ML-ready table:"
)

message(
  file.path(
    OUTDIR,
    "ALL_TISSUES_gene_tissue_limorhyde2_master.tsv"
  )
)


message("")
message(
  "Differential rhythmicity table:"
)

message(
  file.path(
    OUTDIR,
    "ALL_TISSUES_differential_rhythm_stats.tsv"
  )
)


message("")
message(
  "Global developmental-rhythm ranking:"
)

message(
  file.path(
    OUTDIR,
    "ALL_TISSUES_ranked_developmental_rhythm_change.tsv"
  )
)


message("")

