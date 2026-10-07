# ============================================================
# XR-seq DEVELOPMENTAL RHYTHMICITY
#
# BUILD MODELING PHENOTYPE TABLES
#
# Input:
#   XRseq_consensus_master.tsv
#   XRseq_cross_tissue_consensus.tsv
#
# Output:
#
#   01_gene_tissue_modeling_dataset.tsv
#   02_complete_three_tissue_dataset.tsv
#   03_remodeling_binary_dataset.tsv
#   04_gain_vs_loss_dataset.tsv
#   05_tissue_specific_remodeling_dataset.tsv
#
# Plus QC/count tables.
#
# ============================================================
#
# IMPORTANT
# ------------------------------------------------------------
#
# These phenotype labels are derived from the SAME XR-seq
# measurements analyzed by dryR and limorhyde2.
#
# Therefore:
#
#   dryR model/BICW
#   limorhyde2 delta amplitude
#   limorhyde2 RMS difference
#   limorhyde2 phase difference
#
# are phenotype-defining / QC variables.
#
# DO NOT later use them as independent predictors of these
# same phenotype labels.
#
# Future predictors should instead be things such as:
#
#   promoter motifs
#   BMAL1/CLOCK occupancy
#   chromatin state
#   gene length
#   GC content
#   CpG density
#   baseline features chosen before defining the outcome
#   tissue
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  
  library(data.table)
  
})


# ============================================================
# 1. USER SETTINGS -- EDIT THESE IN RSTUDIO
# ============================================================

# Folder produced by:
#
# 03_XRseq_integrate_limorhyde2_dryR.R

INTEGRATION_DIR <- "limorhyde2_dryR_consensus"


# Output directory

OUTDIR <- file.path(
  INTEGRATION_DIR,
  "modeling_phenotypes"
)


dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. INPUT FILES
# ============================================================

GENE_TISSUE_FILE <- file.path(
  INTEGRATION_DIR,
  "XRseq_consensus_master.tsv"
)


CROSS_TISSUE_FILE <- file.path(
  INTEGRATION_DIR,
  "XRseq_cross_tissue_consensus.tsv"
)


if (!file.exists(GENE_TISSUE_FILE)) {
  
  stop(
    "Cannot find:\n",
    GENE_TISSUE_FILE
  )
  
}


if (!file.exists(CROSS_TISSUE_FILE)) {
  
  stop(
    "Cannot find:\n",
    CROSS_TISSUE_FILE
  )
  
}


message("")
message("Gene x tissue consensus:")
message(GENE_TISSUE_FILE)

message("")
message("Cross-tissue consensus:")
message(CROSS_TISSUE_FILE)

message("")
message("Output:")
message(OUTDIR)


# ============================================================
# 3. READ INPUT
# ============================================================

message("")
message("Reading consensus tables...")


gene_tissue <- fread(
  GENE_TISSUE_FILE
)


cross_tissue <- fread(
  CROSS_TISSUE_FILE
)


message(
  "Gene x tissue rows: ",
  nrow(gene_tissue)
)


message(
  "Cross-tissue genes: ",
  nrow(cross_tissue)
)


# ============================================================
# 4. REQUIRED COLUMN CHECKS
# ============================================================

required_gene_tissue <- c(
  
  "gene",
  "tissue",
  "consensus_broad",
  "consensus_subtype",
  "consensus_confidence",
  "chosen_model",
  "chosen_model_BICW",
  "dryR_high_confidence",
  "both_methods"
  
)


missing_gene_tissue <- setdiff(
  required_gene_tissue,
  names(gene_tissue)
)


if (
  length(missing_gene_tissue) > 0
) {
  
  stop(
    
    "Missing required gene-tissue columns:\n",
    
    paste(
      missing_gene_tissue,
      collapse = ", "
    )
    
  )
  
}


required_cross <- c(
  
  "gene",
  
  "state_liver",
  "state_kidney",
  "state_brain",
  
  "subtype_liver",
  "subtype_kidney",
  "subtype_brain",
  
  "cross_tissue_class",
  
  "n_tissues_observed",
  "n_remodeled_tissues"
  
)


missing_cross <- setdiff(
  required_cross,
  names(cross_tissue)
)


if (
  length(missing_cross) > 0
) {
  
  stop(
    
    "Missing required cross-tissue columns:\n",
    
    paste(
      missing_cross,
      collapse = ", "
    )
    
  )
  
}


# ============================================================
# 5. STANDARDIZE BASIC VARIABLES
# ============================================================

gene_tissue[
  ,
  gene :=
    as.character(gene)
]


gene_tissue[
  ,
  tissue :=
    tolower(
      as.character(tissue)
    )
]


cross_tissue[
  ,
  gene :=
    as.character(gene)
]


# ============================================================
# 6. DEFINE STRICT MODELING CLASSES
# ============================================================

CLEAR_CLASSES <- c(
  
  "Stable_nonrhythmic",
  
  "Stable_rhythmic",
  
  "Gain_in_adult",
  
  "Loss_in_adult",
  
  "Altered_rhythm"
  
)


STABLE_CLASSES <- c(
  
  "Stable_nonrhythmic",
  
  "Stable_rhythmic"
  
)


REMODELED_CLASSES <- c(
  
  "Gain_in_adult",
  
  "Loss_in_adult",
  
  "Altered_rhythm"
  
)


EXCLUDED_CLASSES <- c(
  
  "Ambiguous_low_BICW",
  
  "Discordant_loss_direction",
  
  "Discordant_gain_direction",
  
  "Incomplete_method_data",
  
  "Other"
  
)


# ============================================================
# 7. DETERMINE WHETHER EACH GENE x TISSUE ROW IS USABLE
# ============================================================

gene_tissue[
  ,
  usable_for_modeling :=
    consensus_broad %in%
    CLEAR_CLASSES &
    
    both_methods == TRUE &
    
    dryR_high_confidence == TRUE
]


# ============================================================
# 8. EXCLUSION REASON
# ============================================================

gene_tissue[
  ,
  exclusion_reason := fcase(
    
    usable_for_modeling,
    NA_character_,
    
    !both_methods,
    "Incomplete_method_data",
    
    !dryR_high_confidence,
    "Low_dryR_BICW",
    
    grepl(
      "^Discordant",
      consensus_broad
    ),
    "dryR_limorhyde2_direction_discordance",
    
    consensus_broad ==
      "Other",
    "Other_consensus_class",
    
    default =
      paste0(
        "Excluded_",
        consensus_broad
      )
    
  )
]


# ============================================================
# 9. WRITE EXCLUDED ROWS
# ============================================================

excluded_rows <- gene_tissue[
  usable_for_modeling == FALSE
]


fwrite(
  
  excluded_rows,
  
  file.path(
    OUTDIR,
    "QC_excluded_gene_tissue_rows.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 10. CREATE STRICT GENE x TISSUE MODELING DATASET
# ============================================================

modeling <- copy(
  
  gene_tissue[
    usable_for_modeling == TRUE
  ]
  
)


# ------------------------------------------------------------
# Stable vs remodeled target
# ------------------------------------------------------------

modeling[
  ,
  target_remodeled :=
    fifelse(
      consensus_broad %in%
        REMODELED_CLASSES,
      1L,
      0L
    )
]


modeling[
  ,
  target_remodeled_label :=
    fifelse(
      target_remodeled == 1L,
      "Remodeled",
      "Stable"
    )
]


# ------------------------------------------------------------
# Five-class phenotype
# ------------------------------------------------------------

modeling[
  ,
  phenotype_5class :=
    consensus_broad
]


modeling[
  ,
  phenotype_5class :=
    factor(
      
      phenotype_5class,
      
      levels = c(
        
        "Stable_nonrhythmic",
        
        "Stable_rhythmic",
        
        "Gain_in_adult",
        
        "Loss_in_adult",
        
        "Altered_rhythm"
        
      )
      
    )
]


# ------------------------------------------------------------
# Stable-state subtype
# ------------------------------------------------------------

modeling[
  ,
  stable_state :=
    fcase(
      
      consensus_broad ==
        "Stable_nonrhythmic",
      "Stable_nonrhythmic",
      
      consensus_broad ==
        "Stable_rhythmic",
      "Stable_rhythmic",
      
      default =
        NA_character_
      
    )
]


# ------------------------------------------------------------
# Remodeling subtype
# ------------------------------------------------------------

modeling[
  ,
  remodeling_type :=
    fcase(
      
      consensus_broad ==
        "Gain_in_adult",
      "Gain",
      
      consensus_broad ==
        "Loss_in_adult",
      "Loss",
      
      consensus_broad ==
        "Altered_rhythm",
      "Altered",
      
      default =
        NA_character_
      
    )
]


# ------------------------------------------------------------
# Group ID
#
# Later grouped cross-validation MUST keep all tissues
# belonging to the same gene together.
# ------------------------------------------------------------

modeling[
  ,
  group_id :=
    gene
]


# ============================================================
# 11. QC LABELS
#
# Rename several phenotype-defining variables so they are
# visibly NOT ordinary predictor columns.
# ============================================================

modeling[
  ,
  qc_dryR_model :=
    as.character(
      dryR_model
    )
]


modeling[
  ,
  qc_dryR_BICW :=
    chosen_model_BICW
]


modeling[
  ,
  qc_consensus_subtype :=
    consensus_subtype
]


# ============================================================
# 12. ORDER MAIN COLUMNS
# ============================================================

main_columns <- c(
  
  "gene",
  "tissue",
  "group_id",
  
  "phenotype_5class",
  
  "target_remodeled",
  "target_remodeled_label",
  
  "stable_state",
  "remodeling_type",
  
  "qc_consensus_subtype",
  "qc_dryR_model",
  "qc_dryR_BICW"
  
)


main_columns <- main_columns[
  main_columns %in%
    names(modeling)
]


setcolorder(
  
  modeling,
  
  c(
    
    main_columns,
    
    setdiff(
      names(modeling),
      main_columns
    )
    
  )
  
)


setorder(
  modeling,
  tissue,
  gene
)


fwrite(
  
  modeling,
  
  file.path(
    OUTDIR,
    "01_gene_tissue_modeling_dataset.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 13. GENE x TISSUE MODELING COUNTS
# ============================================================

modeling_counts <- modeling[
  ,
  .N,
  by = .(
    tissue,
    phenotype_5class
  )
]


setnames(
  modeling_counts,
  "N",
  "n_gene_tissue_rows"
)


modeling_counts[
  ,
  percent :=
    100 *
    n_gene_tissue_rows /
    sum(
      n_gene_tissue_rows
    ),
  by = tissue
]


fwrite(
  
  modeling_counts,
  
  file.path(
    OUTDIR,
    "QC_gene_tissue_class_counts.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 14. DEFINE COMPLETE THREE-TISSUE GENES
# ============================================================

is_clear_state <- function(x) {
  
  !is.na(x) &
    x %in%
    CLEAR_CLASSES
  
}


cross_tissue[
  ,
  liver_clear :=
    is_clear_state(
      state_liver
    )
]


cross_tissue[
  ,
  kidney_clear :=
    is_clear_state(
      state_kidney
    )
]


cross_tissue[
  ,
  brain_clear :=
    is_clear_state(
      state_brain
    )
]


cross_tissue[
  ,
  all_three_clear :=
    liver_clear &
    kidney_clear &
    brain_clear
]


complete_three <- copy(
  
  cross_tissue[
    all_three_clear == TRUE
  ]
  
)


# ============================================================
# 15. ADD SIMPLIFIED THREE-TISSUE LABELS
# ============================================================

complete_three[
  ,
  liver_remodeled :=
    as.integer(
      state_liver %in%
        REMODELED_CLASSES
    )
]


complete_three[
  ,
  kidney_remodeled :=
    as.integer(
      state_kidney %in%
        REMODELED_CLASSES
    )
]


complete_three[
  ,
  brain_remodeled :=
    as.integer(
      state_brain %in%
        REMODELED_CLASSES
    )
]


# Recompute for QC rather than trusting previous value.

complete_three[
  ,
  n_remodeled_tissues_recalculated :=
    liver_remodeled +
    kidney_remodeled +
    brain_remodeled
]


# ============================================================
# 16. HIGH-LEVEL CROSS-TISSUE CATEGORY
# ============================================================

complete_three[
  ,
  modeling_cross_tissue_group :=
    fcase(
      
      n_remodeled_tissues_recalculated == 0,
      "Stable_all_tissues",
      
      n_remodeled_tissues_recalculated == 1,
      "One_tissue_remodeled",
      
      n_remodeled_tissues_recalculated == 2,
      "Two_tissues_remodeled",
      
      n_remodeled_tissues_recalculated == 3,
      "All_three_remodeled",
      
      default =
        NA_character_
      
    )
]


# ============================================================
# 17. SHARED DIRECTION LABEL
# ============================================================

complete_three[
  ,
  shared_direction :=
    fcase(
      
      state_liver ==
        "Gain_in_adult" &
        state_kidney ==
        "Gain_in_adult" &
        state_brain ==
        "Gain_in_adult",
      "Shared_gain",
      
      state_liver ==
        "Loss_in_adult" &
        state_kidney ==
        "Loss_in_adult" &
        state_brain ==
        "Loss_in_adult",
      "Shared_loss",
      
      state_liver ==
        "Altered_rhythm" &
        state_kidney ==
        "Altered_rhythm" &
        state_brain ==
        "Altered_rhythm",
      "Shared_altered",
      
      n_remodeled_tissues_recalculated == 3,
      "Mixed_remodeling_all_three",
      
      default =
        NA_character_
      
    )
]


# ============================================================
# 18. OPPOSITE GAIN/LOSS FLAG
# ============================================================

complete_three[
  ,
  opposite_gain_loss :=
    (
      
      state_liver ==
        "Gain_in_adult" |
        
        state_kidney ==
        "Gain_in_adult" |
        
        state_brain ==
        "Gain_in_adult"
      
    ) &
    (
      
      state_liver ==
        "Loss_in_adult" |
        
        state_kidney ==
        "Loss_in_adult" |
        
        state_brain ==
        "Loss_in_adult"
      
    )
]


# ============================================================
# 19. COMPLETE THREE-TISSUE OUTPUT
# ============================================================

setorder(
  
  complete_three,
  
  -n_remodeled_tissues_recalculated,
  
  gene
  
)


fwrite(
  
  complete_three,
  
  file.path(
    OUTDIR,
    "02_complete_three_tissue_dataset.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 20. COMPLETE THREE-TISSUE SUMMARY
# ============================================================

complete_summary <- complete_three[
  ,
  .N,
  by = .(
    modeling_cross_tissue_group
  )
]


setnames(
  complete_summary,
  "N",
  "n_genes"
)


complete_summary[
  ,
  percent :=
    100 *
    n_genes /
    sum(n_genes)
]


fwrite(
  
  complete_summary,
  
  file.path(
    OUTDIR,
    "QC_complete_three_tissue_summary.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 21. DATASET 3:
# BINARY REMODELING
#
# Question:
#
# What predicts whether developmental XR-seq rhythmicity
# is remodeled at all?
#
# 0 = Stable
# 1 = Remodeled
#
# Stable:
#
#   Stable_nonrhythmic
#   Stable_rhythmic
#
# Remodeled:
#
#   Gain
#   Loss
#   Altered
#
# ============================================================

binary_dataset <- modeling[
  ,
  .(
    
    gene,
    
    tissue,
    
    group_id,
    
    phenotype_5class =
      as.character(
        phenotype_5class
      ),
    
    consensus_subtype,
    
    target_remodeled,
    
    target_label =
      target_remodeled_label,
    
    # ----------------------------------------------
    # QC / provenance only
    # DO NOT use these as predictors.
    # ----------------------------------------------
    
    qc_dryR_model =
      as.character(
        dryR_model
      ),
    
    qc_dryR_BICW =
      chosen_model_BICW
    
  )
]


setorder(
  binary_dataset,
  tissue,
  gene
)


fwrite(
  
  binary_dataset,
  
  file.path(
    OUTDIR,
    "03_remodeling_binary_dataset.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 22. BINARY TARGET COUNTS
# ============================================================

binary_counts <- binary_dataset[
  ,
  .N,
  by = .(
    tissue,
    target_label
  )
]


setnames(
  binary_counts,
  "N",
  "n_rows"
)


binary_counts[
  ,
  percent :=
    100 *
    n_rows /
    sum(n_rows),
  by = tissue
]


fwrite(
  
  binary_counts,
  
  file.path(
    OUTDIR,
    "QC_binary_remodeling_counts.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 23. DATASET 4:
# GAIN vs LOSS
#
# Altered_rhythm is intentionally excluded.
#
# Question:
#
# Among genes with directional developmental remodeling,
# what predicts whether rhythmicity is gained or lost?
#
# 0 = Loss
# 1 = Gain
#
# ============================================================

gain_loss <- modeling[
  consensus_broad %in%
    c(
      "Gain_in_adult",
      "Loss_in_adult"
    )
]


gain_loss[
  ,
  target_gain :=
    as.integer(
      consensus_broad ==
        "Gain_in_adult"
    )
]


gain_loss[
  ,
  target_direction :=
    fifelse(
      target_gain == 1L,
      "Gain",
      "Loss"
    )
]


gain_loss_dataset <- gain_loss[
  ,
  .(
    
    gene,
    
    tissue,
    
    group_id,
    
    target_gain,
    
    target_direction,
    
    consensus_subtype,
    
    # ----------------------------------------------
    # QC only
    # ----------------------------------------------
    
    qc_dryR_model =
      as.character(
        dryR_model
      ),
    
    qc_dryR_BICW =
      chosen_model_BICW
    
  )
]


setorder(
  gain_loss_dataset,
  tissue,
  gene
)


fwrite(
  
  gain_loss_dataset,
  
  file.path(
    OUTDIR,
    "04_gain_vs_loss_dataset.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 24. GAIN vs LOSS COUNTS
# ============================================================

gain_loss_counts <- gain_loss_dataset[
  ,
  .N,
  by = .(
    tissue,
    target_direction
  )
]


setnames(
  gain_loss_counts,
  "N",
  "n_rows"
)


gain_loss_counts[
  ,
  percent :=
    100 *
    n_rows /
    sum(n_rows),
  by = tissue
]


fwrite(
  
  gain_loss_counts,
  
  file.path(
    OUTDIR,
    "QC_gain_vs_loss_counts.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 25. DATASET 5:
# TISSUE-SPECIFIC REMODELING
#
# Restrict to genes:
#
#   - confidently classified in ALL THREE tissues
#   - remodeled in exactly 1 or 2 tissues
#
# Therefore each gene contributes:
#
#   at least one remodeled tissue
#   AND
#   at least one stable tissue
#
# This is crucial because it lets us compare tissue context
# within the SAME gene.
#
# ============================================================

tissue_specific_genes <- complete_three[
  n_remodeled_tissues_recalculated %in%
    c(
      1L,
      2L
    ),
  gene
]


message("")
message(
  "Complete genes with tissue-dependent remodeling: ",
  length(
    tissue_specific_genes
  )
)


# ============================================================
# 26. PULL GENE x TISSUE ROWS FOR THOSE GENES
# ============================================================

tissue_specific <- modeling[
  gene %in%
    tissue_specific_genes
]


# ------------------------------------------------------------
# Every selected gene MUST have exactly 3 rows.
# ------------------------------------------------------------

ts_gene_qc <- tissue_specific[
  ,
  .N,
  by = gene
]


if (
  any(
    ts_gene_qc$N != 3L
  )
) {
  
  bad <- ts_gene_qc[
    N != 3L
  ]
  
  
  fwrite(
    
    bad,
    
    file.path(
      OUTDIR,
      "ERROR_tissue_specific_gene_row_counts.tsv"
    ),
    
    sep = "\t"
    
  )
  
  
  stop(
    "Some complete-three-tissue genes do not have exactly ",
    "three usable gene-tissue rows."
  )
  
}


# ============================================================
# 27. ADD CROSS-TISSUE INFORMATION
# ============================================================

cross_annotation <- complete_three[
  gene %in%
    tissue_specific_genes,
  .(
    
    gene,
    
    cross_tissue_class,
    
    modeling_cross_tissue_group,
    
    n_remodeled_tissues =
      n_remodeled_tissues_recalculated,
    
    opposite_gain_loss
    
  )
]


tissue_specific <- merge(
  
  tissue_specific,
  
  cross_annotation,
  
  by = "gene",
  
  all.x = TRUE,
  
  sort = FALSE
  
)


# ============================================================
# 28. TISSUE-SPECIFIC TARGET
#
# 1 = this particular tissue is remodeled
# 0 = same gene is stable in this particular tissue
# ============================================================

tissue_specific[
  ,
  target_remodeled_in_tissue :=
    as.integer(
      consensus_broad %in%
        REMODELED_CLASSES
    )
]


tissue_specific[
  ,
  target_tissue_state :=
    fifelse(
      
      target_remodeled_in_tissue == 1L,
      
      "Remodeled",
      
      "Stable"
      
    )
]


# ============================================================
# 29. TISSUE-SPECIFIC PATTERN
# ============================================================

tissue_specific[
  ,
  tissue_specific_pattern :=
    fifelse(
      
      n_remodeled_tissues == 1L,
      
      "One_of_three_remodeled",
      
      "Two_of_three_remodeled"
      
    )
]


# ============================================================
# 30. SELECT MODELING-SAFE COLUMNS
# ============================================================

tissue_specific_dataset <- tissue_specific[
  ,
  .(
    
    gene,
    
    tissue,
    
    group_id =
      gene,
    
    target_remodeled_in_tissue,
    
    target_tissue_state,
    
    tissue_specific_pattern,
    
    n_remodeled_tissues,
    
    cross_tissue_class,
    
    opposite_gain_loss,
    
    phenotype_5class =
      as.character(
        phenotype_5class
      ),
    
    consensus_subtype,
    
    # ----------------------------------------------
    # QC only
    # ----------------------------------------------
    
    qc_dryR_model =
      as.character(
        dryR_model
      ),
    
    qc_dryR_BICW =
      chosen_model_BICW
    
  )
]


setorder(
  tissue_specific_dataset,
  gene,
  tissue
)


fwrite(
  
  tissue_specific_dataset,
  
  file.path(
    OUTDIR,
    "05_tissue_specific_remodeling_dataset.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 31. TISSUE-SPECIFIC TARGET COUNTS
# ============================================================

tissue_specific_counts <- tissue_specific_dataset[
  ,
  .N,
  by = .(
    tissue,
    target_tissue_state
  )
]


setnames(
  tissue_specific_counts,
  "N",
  "n_rows"
)


tissue_specific_counts[
  ,
  percent :=
    100 *
    n_rows /
    sum(n_rows),
  by = tissue
]


fwrite(
  
  tissue_specific_counts,
  
  file.path(
    OUTDIR,
    "QC_tissue_specific_target_counts.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 32. OVERALL DATASET SUMMARY
# ============================================================

dataset_summary <- data.table(
  
  dataset = c(
    
    "01_gene_tissue_modeling_dataset",
    
    "02_complete_three_tissue_dataset",
    
    "03_remodeling_binary_dataset",
    
    "04_gain_vs_loss_dataset",
    
    "05_tissue_specific_remodeling_dataset"
    
  ),
  
  n_rows = c(
    
    nrow(
      modeling
    ),
    
    nrow(
      complete_three
    ),
    
    nrow(
      binary_dataset
    ),
    
    nrow(
      gain_loss_dataset
    ),
    
    nrow(
      tissue_specific_dataset
    )
    
  ),
  
  n_unique_genes = c(
    
    uniqueN(
      modeling$gene
    ),
    
    uniqueN(
      complete_three$gene
    ),
    
    uniqueN(
      binary_dataset$gene
    ),
    
    uniqueN(
      gain_loss_dataset$gene
    ),
    
    uniqueN(
      tissue_specific_dataset$gene
    )
    
  )
  
)


fwrite(
  
  dataset_summary,
  
  file.path(
    OUTDIR,
    "00_dataset_summary.tsv"
  ),
  
  sep = "\t"
  
)


# ============================================================
# 33. LABEL SUMMARY
# ============================================================

# ------------------------------------------------------------
# Gene x tissue 5-class phenotype
# ------------------------------------------------------------

label_5class <- modeling[
  ,
  .N,
  by = .(
    tissue,
    label =
      as.character(
        phenotype_5class
      )
  )
]

label_5class[
  ,
  dataset :=
    "gene_tissue_5class"
]


# ------------------------------------------------------------
# Binary remodeling phenotype
# ------------------------------------------------------------

label_binary <- binary_dataset[
  ,
  .N,
  by = .(
    tissue,
    label =
      target_label
  )
]

label_binary[
  ,
  dataset :=
    "binary_remodeling"
]


# ------------------------------------------------------------
# Gain vs Loss phenotype
# ------------------------------------------------------------

label_gain_loss <- gain_loss_dataset[
  ,
  .N,
  by = .(
    tissue,
    label =
      target_direction
  )
]

label_gain_loss[
  ,
  dataset :=
    "gain_vs_loss"
]


# ------------------------------------------------------------
# Tissue-specific remodeling phenotype
# ------------------------------------------------------------

label_tissue_specific <- tissue_specific_dataset[
  ,
  .N,
  by = .(
    tissue,
    label =
      target_tissue_state
  )
]

label_tissue_specific[
  ,
  dataset :=
    "tissue_specific"
]


# ------------------------------------------------------------
# Combine all label summaries
# ------------------------------------------------------------

label_summary <- rbindlist(
  
  list(
    label_5class,
    label_binary,
    label_gain_loss,
    label_tissue_specific
  ),
  
  use.names = TRUE,
  fill = TRUE
  
)


setnames(
  label_summary,
  "N",
  "n_rows"
)


# ------------------------------------------------------------
# Put columns in intuitive order
# ------------------------------------------------------------

setcolorder(
  
  label_summary,
  
  c(
    "dataset",
    "tissue",
    "label",
    "n_rows"
  )
  
)


setorder(
  
  label_summary,
  
  dataset,
  tissue,
  label
  
)


fwrite(
  
  label_summary,
  
  file.path(
    OUTDIR,
    "00_label_summary.tsv"
  ),
  
  sep = "\t"
  
)
# ============================================================
# 34. WRITE ML LEAKAGE WARNING
# ============================================================

warning_text <- c(
  
  "XR-seq MODELING PHENOTYPE DATASETS",
  "",
  "IMPORTANT: TARGET LEAKAGE",
  "",
  "The phenotype labels in these files were generated using",
  "dryR and limorhyde2 from the XR-seq time-course data.",
  "",
  "Therefore the following MUST NOT automatically be used",
  "as independent predictors of these same labels:",
  "",
  "  dryR chosen_model",
  "  dryR chosen_model_BICW",
  "  limorhyde2 diff_peak_trough_amp",
  "  limorhyde2 diff_peak_phase",
  "  limorhyde2 diff_mesor",
  "  limorhyde2 diff_rhy_dist",
  "  limorhyde2 rms_diff_rhy",
  "",
  "These are phenotype-defining/QC variables.",
  "",
  "Recommended independent predictors include:",
  "",
  "  promoter sequence features",
  "  E-box / D-box / RRE motif features",
  "  BMAL1/CLOCK occupancy",
  "  chromatin accessibility",
  "  histone modifications",
  "  GC content",
  "  CpG density",
  "  gene length",
  "  promoter architecture",
  "  tissue identity",
  "",
  "For cross-validation:",
  "",
  "ALL rows belonging to the same gene must stay in the",
  "same fold. Never randomly split gene x tissue rows.",
  "",
  "Use grouped CV by gene and leave-one-tissue-out validation."
  
)


writeLines(
  
  warning_text,
  
  con = file.path(
    OUTDIR,
    "README_MODELING_TARGETS.txt"
  )
  
)


# ============================================================
# 35. SESSION INFO
# ============================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
  
)


# ============================================================
# 36. PRINT SUMMARY
# ============================================================

message("")
message(
  "===================================================="
)

message(
  "MODELING PHENOTYPE DATASETS COMPLETE"
)

message(
  "===================================================="
)

message("")

print(
  dataset_summary
)


message("")
message(
  "Main outputs:"
)

message("")

message(
  file.path(
    OUTDIR,
    "01_gene_tissue_modeling_dataset.tsv"
  )
)

message(
  file.path(
    OUTDIR,
    "02_complete_three_tissue_dataset.tsv"
  )
)

message(
  file.path(
    OUTDIR,
    "03_remodeling_binary_dataset.tsv"
  )
)

message(
  file.path(
    OUTDIR,
    "04_gain_vs_loss_dataset.tsv"
  )
)

message(
  file.path(
    OUTDIR,
    "05_tissue_specific_remodeling_dataset.tsv"
  )
)

message("")
message(
  "Done."
)