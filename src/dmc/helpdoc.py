#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Help text used by the CpGtools ``epical`` command-line interface."""

GENERAL_HELP = "Epical: a DNA methylation-based epigenetic age calculator."

FORMAT_HELP = (
    "Output plot format. Choose either 'pdf' or 'png'. "
    "Default: 'pdf'."
)

INPUT_HELP = (
    "Input tabular file containing DNA methylation data. The file must have "
    "a header row containing sample names or labels. The first column must "
    "contain CpG IDs, and the remaining columns must contain DNA methylation "
    "beta values as floating-point numbers between 0 and 1. Columns may be "
    "separated by tabs, commas, or another delimiter. Use 'NA' or 'NaN' for "
    "missing values. Plain-text and compressed files are supported, including "
    ".gz, .Z, .z, .bz, .bz2, and .bzip2 files."
)

OUTPUT_HELP = (
    "Prefix for output files. If omitted, the default prefix "
    "'clock_name_out' is used. Generated files may include: "
    "<PREFIX>.DNAm_age.tsv (predicted DNAm age); "
    "<PREFIX>.used_CpGs.tsv (CpGs used for prediction); "
    "<PREFIX>.missed_CpGs.txt (clock CpGs missing or excluded from input); "
    "<PREFIX>.coef.tsv (clock CpGs and coefficients, including whether each "
    "CpG was used); <PREFIX>.plots.R (R script for visualization); and "
    "<PREFIX>.coef_plot.pdf or .png (coefficient plot)."
)

NA_HELP = (
    "Maximum allowable fraction of missing clock CpGs. Default: 0.2 (20%%). "
    "If more than this fraction of required clock CpGs is missing, DNAm age "
    "estimation is not performed."
)

DEL_HELP = (
    "Column delimiter used in the input file, typically a tab or comma. "
    "If omitted, the delimiter is detected automatically."
)

DEBUG_HELP = "Print detailed diagnostic information for debugging."

META_HELP = (
    "Optional sample metadata file. The file must have a header row naming "
    "variables such as 'Sex' and 'Age', and the first column must contain "
    "sample IDs. If an 'Age' column is present, a scatter plot comparing "
    "chronological age with predicted DNAm age is generated."
)

EPM_META_HELP = (
    "Metadata file for EPM. The file must have a header row naming variables. "
    "An 'Age' column is required. A 'Designation' column may be used to label "
    "training and testing samples."
)

EPM_HELP = (
    "Epigenetic Pacemaker (EPM): a conditional expectation-maximization "
    "algorithm for modeling epigenetic states under an evolutionary framework. "
    "Unlike standard linear-regression clocks, EPM does not assume a linear "
    "relationship between epigenetic state and the trait of interest. "
    "Reference: Farrell C, et al. Bioinformatics (2020). "
    "PubMed: https://pubmed.ncbi.nlm.nih.gov/32573701/."
)

WLMT_HELP = (
    "WLMT is a whole-lifespan, multi-tissue mouse epigenetic age predictor "
    "based on 435 CpG sites identified from RRBS/WGBS data. Mouse clock CpG "
    "IDs use the form 'chrom_position', where position is the 1-based "
    "coordinate of cytosine. WLMT expects beta values on a 0-100 scale; input "
    "values on a 0-1 scale are converted automatically. Reference: Meer MV, "
    "et al. eLife (2018). PubMed: https://pubmed.ncbi.nlm.nih.gov/30427307/."
)

YOMT_HELP = (
    "YOMT is a multi-tissue mouse epigenetic age predictor based on 329 CpG "
    "sites identified from RRBS/WGBS data. Mouse clock CpG IDs use the form "
    "'chrom_position', where position is the 1-based coordinate of cytosine. "
    "Reference: Stubbs TM, et al. Genome Biology (2017). "
    "PubMed: https://pubmed.ncbi.nlm.nih.gov/28399939/."
)

MMLIVER_HELP = (
    "Mouse liver epigenetic age predictor based on 148 CpG sites. Mouse clock "
    "CpG IDs use the form 'chrom_position', where position is the 1-based "
    "coordinate of cytosine. Reference: Wang T, et al. Genome Biology (2017). "
    "PubMed: https://pubmed.ncbi.nlm.nih.gov/28351423/."
)

MMBLOOD_HELP = (
    "Mouse blood epigenetic age predictor based on 90 CpG sites. Mouse clock "
    "CpG IDs use the form 'chrom_position', where position is the 1-based "
    "coordinate of cytosine. Reference: Petkovich DA, et al. Cell Metabolism "
    "(2017). PubMed: https://pubmed.ncbi.nlm.nih.gov/28380383/."
)

EPM_OUTPUT_HELP = (
    "Prefix for EPM output files. If omitted, the default prefix 'EPM_out' is "
    "used. Outputs may include predicted EPM ages for training and testing "
    "samples, scatter plots comparing predicted and chronological ages, and "
    "selected CpG beta-value tables for training and testing samples."
)

LOG_HELP = (
    "Optional log file. If omitted, log messages are written to the terminal."
)


# Backward-compatible lowercase aliases used by existing CpGtools code.
general_help = GENERAL_HELP
format_help = FORMAT_HELP
input_help = INPUT_HELP
output_help = OUTPUT_HELP
na_help = NA_HELP
del_help = DEL_HELP
debug_help = DEBUG_HELP
meta_help = META_HELP
epm_meta_help = EPM_META_HELP
epm_help = EPM_HELP
WLMT_help = WLMT_HELP
YOMT_help = YOMT_HELP
mmLiver_help = MMLIVER_HELP
mmBlood_help = MMBLOOD_HELP
epm_output_help = EPM_OUTPUT_HELP
log_help = LOG_HELP
