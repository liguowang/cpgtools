#!/usr/bin/env python3
#
# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: "deep.omics.lab@gmail.com"
# Project: https://github.com/liguowang/cpgtools
#
# This file is part of CpGtools and is distributed under the MIT License.
# See the LICENSE.txt file in the project root for the full license text.

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

BIOLEARN_CLOCK_HELP = {
    "YingCausAge":
        "Causality-enriched blood DNAm clock estimating age from CpGs with evidence of causal effects on aging-related traits. (PMID: 38243142)",

    "YingDamAge":
        "Causality-enriched blood DNAm clock emphasizing CpGs whose age-related changes are predicted to be damaging and associated with adverse aging outcomes. (PMID: 38243142)",

    "YingAdaptAge":
        "Causality-enriched blood DNAm clock emphasizing CpGs whose age-related changes may represent adaptive or protective responses to aging-related damage. (PMID: 38243142)",

    "DunedinPoAm38":
        "A 46-CpG blood DNAm measure of Pace of Aging, trained against longitudinal change in 18 biomarkers of organ-system integrity measured at ages 26, 32, and 38; unlike age clocks, it estimates aging rate rather than accumulated biological age. (PMID: 32367804)",

    "GrimAgeV1":
        "Second-generation blood DNAm biomarker of mortality-associated biological age, combining DNAm surrogates for smoking pack-years and mortality-related plasma proteins to predict lifespan and healthspan. (PMID: 30669119)",

    "GrimAgeV2":
        "Updated GrimAge mortality-associated biological-age biomarker that extends GrimAgeV1 by adding DNAm surrogates for high-sensitivity C-reactive protein (logCRP) and hemoglobin A1C (logA1C). (PMID: 36516495)",

    "VidalBralo":
        "Minimal blood DNAm clock that estimates chronological age from methylation at 8 CpG sites, designed as a simplified assay for adult age estimation. (PMID: 27471517)",

    "AlcoholMcCartney":
        "Blood DNAm-based risk score trained to predict alcohol consumption as part of a family of methylation predictors for lifestyle and health-related traits. (PMID: 30257690)",

    "HRSInCHPhenoAge":
        "Principal-component version of DNAm PhenoAge designed to reduce technical noise and improve reliability for longitudinal and repeated-measure applications. (PMID: 36277076)",

    "EpiTOC1":
        "Mitotic DNAm clock based on 385 Polycomb target CpGs whose methylation increases with age, providing a relative measure of cumulative stem-cell divisions and tissue mitotic age. (PMID: 27716309)",

    "EpiTOC2":
        "Second-generation mitotic DNAm clock using 163 CpGs and an explicit DNAm transmission model to estimate the cumulative number and rate of stem-cell divisions in a tissue. (PMID: 32580750)",

    "SmokingMcCartney":
        "Blood DNAm-based smoking score calculated from 233 CpGs and trained as a methylation predictor of smoking exposure/status. (PMID: 30257690)",

    "DownSyndrome":
        "DNAm score capturing the genome-wide epigenetic signature of Down syndrome (trisomy 21), derived from CpGs differentially methylated between newborns with and without Down syndrome. (PMID: 33547282)",

    "StocZ":
        "Fully stochastic chronological-age clock built from simulated age-related DNAm trajectories at the 514 CpGs of the Zhang clock, designed to quantify the stochastic component underlying Zhang clock age prediction. (PMID: 38724732)",

    "StocP":
        "Fully stochastic chronological-age clock built from simulated age-related DNAm trajectories at the 513 CpGs of PhenoAge, designed to quantify the stochastic component underlying PhenoAge age prediction. (PMID: 38724732)",

    "StocH":
        "Fully stochastic chronological-age clock built from simulated age-related DNAm trajectories at the 353 CpGs of Horvath's multi-tissue clock, designed to quantify the stochastic component underlying Horvath clock age prediction. (PMID: 38724732)",

    "BMI_McCartney":
        "Blood DNAm-based risk score trained to predict body mass index (BMI) as part of a family of methylation predictors for lifestyle and health-related traits. (PMID: 30257690)",

    "EducationMcCartney":
        "Blood DNAm-based risk score trained to predict educational attainment as part of a family of methylation predictors for lifestyle and health-related traits. (PMID: 30257690)",

    "TotalCholesterolMcCartney":
        "Blood DNAm-based risk score trained to predict total cholesterol as part of a family of methylation predictors for cardiometabolic traits. (PMID: 30257690)",

    "HDLCholesterolMcCartney":
        "Blood DNAm-based risk score trained to predict HDL cholesterol as part of a family of methylation predictors for cardiometabolic traits. (PMID: 30257690)",

    "LDLCholesterolMcCartney":
        "Blood DNAm-based risk score trained to predict LDL plus remnant cholesterol as part of a family of methylation predictors for cardiometabolic traits. (PMID: 30257690)",

    "BodyFatMcCartney":
        "Blood DNAm-based risk score trained to predict percentage body fat as part of a family of methylation predictors for lifestyle and health-related traits. (PMID: 30257690)",

    "BMI_Reed":
        "Blood DNAm predictor of body mass index (BMI), developed to capture methylation patterns associated with obesity and adiposity. (PMID: 32228717)",

    "ProstateCancerKirby":
        "Prostate-tissue DNAm classifier developed to distinguish prostate cancer from non-cancer samples using cancer-associated methylation patterns. (PMID: 28412973)",

    "HepatoXu":
        "Circulating cell-free DNA methylation classifier developed for detection and diagnosis of hepatocellular carcinoma. (PMID: 29035356)",

    "CVD_Westerman":
        "Blood DNAm-based model developed to predict incident coronary heart disease and capture methylation signatures associated with cardiovascular risk. (PMID: 32308120)",

    "AD_Bahado-Singh":
        "Blood DNAm-based classifier developed to distinguish individuals with Alzheimer's disease from controls using disease-associated methylation signatures. (PMID: 33788842)",

    "DepressionBarbu":
        "Blood DNAm-based risk score developed from methylation differences associated with major depressive disorder and depression risk. (PMID: 32523041)",

    "Bocklandt":
        "Early age predictor based on methylation at two CpG sites associated with EDARADD and NPTX2, explaining about 73%% of chronological-age variation in the original study. (PMID: 21731603)",
}
    
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
