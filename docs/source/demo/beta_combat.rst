beta_combat.py
============================

Correct beta-value matrices for technical (batch) effects using the
`ComBat <https://pubmed.ncbi.nlm.nih.gov/16632515/>`_ algorithm.

Overview
--------

``beta_combat.py`` takes a CpG × sample beta matrix and a
sample-to-batch mapping, applies ComBat, and writes an adjusted matrix
plus before/after QC boxplots.

Synopsis
--------

::

   beta_combat.py \
     -i <beta_matrix.tsv[.gz]> \
     -g <batch_map.csv> \
     -o <output_prefix>

Input files
-----------

Beta matrix (TSV/TSV.GZ)
~~~~~~~~~~~~~~~~~~~~~~~~

- **Delimiter:** tab  
- **Header:** first row contains sample IDs  
- **Index:** first column contains CpG IDs  
- **Values:** beta values in [0, 1]  
- **Missing values:** rows containing any NA/empty cells are removed prior to ComBat  

Example::

   CpG_ID   Sample_01  Sample_02  Sample_03  Sample_04
   cg_001   0.831035   0.878022   0.794427   0.880911
   cg_002   0.249544   0.209949   0.234294   0.236680
   cg_003   0.845065   0.843957   0.840184   0.824286
   ...

Batch map (CSV)
~~~~~~~~~~~~~~~

- **Delimiter:** comma  
- **Columns:** ``Sample,Group``  
- **Sample IDs:** must match the sample IDs in the beta matrix header (case-sensitive)  
- **Grouping:** each sample belongs to exactly one batch group (e.g., plates, chips)  

Example::

   Sample,Group
   Sample_01,plate_1
   Sample_02,plate_1
   Sample_03,plate_2
   Sample_04,plate_2
   ...

Example input files
~~~~~~~~~~~~~~~~~~~

- `test_12_threebatch.beta.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_12_threebatch.beta.tsv.gz>`_
- `test_12_threebatch.beta.100K_NAs.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_12_threebatch.beta.100K_NAs.tsv.gz/>`_ (with 100,000 missign values)
- `test_12_threebatch.batch.csv <https://sourceforge.net/projects/cpgtools/files/test/test_12_threebatch.batch.csv>`_

Options
-------

::

   --version             show program's version number and exit
   -h, --help            show this help message and exit
   -i INPUT_FILE, --input_file=INPUT_FILE
                         Tab-separated data frame file containing beta values
                         with the 1st row containing sample IDs and the 1st
                         column containing CpG IDs.
   --axis=AXIS_CHOICE    How to do imputation (using the KNN algorithm) if the
                         input file has missing values. 1: search columns for k
                         nearest neighbours; 0: Search rows for k nearest
                         neighbours. default=1
   -g GROUP_FILE, --group=GROUP_FILE
                         Comma-separated group file defining the batch groups
                         of each sample.
   -o OUT_FILE, --output=OUT_FILE
                         The prefix of the output file.

Command example (input file has no missing values)
--------------------------------------------------

::

   $ beta_combat.py \
       -i test_12_threebatch.beta.tsv.gz \
       -g test_12_threebatch.batch.csv \
       -o output

Command example (input file has missing values)
--------------------------------------------------

If the input file has missing values, KNN will be used to impute the missing values first before batch effect correction.

::

   $ beta_combat.py \
       -i test_12_threebatch.beta.100K_NAs.tsv.gz \
       -g test_12_threebatch.batch.csv \
       -o output


Outputs (input file has no missing values)
-------------------------------------------

- ``<prefix>.combat.tsv`` — beta matrix after ComBat batch correction
- ``<prefix>.boxplot.png`` — distribution of beta values **before** batch effect correction  
- ``<prefix>.boxplot_combat.png`` — distribution of beta values **after** batch effect correction  

Outputs (input file with missing values)
------------------------------------------

- ``<prefix>.combat.tsv`` — beta matrix after ComBat batch correction (missing vlaues are predicted using KNN)
- ``<prefix>.combat_withNAs.tsv`` - beta matrix after ComBat batch correction (keep missing values)
- ``<prefix>.boxplot.png`` — distribution of beta values **before** batch effect correction  
- ``<prefix>.boxplot_combat.png`` — distribution of beta values **after** batch effect correction  

Figures
-------

.. image:: ../_static/output.boxplot.png
   :height: 400px
   :width: 600px
   :alt: Boxplot of beta values before ComBat

.. image:: ../_static/output.boxplot_combat.png
   :height: 400px
   :width: 600px
   :alt: Boxplot of beta values after ComBat

Notes & tips
------------

- Ensure all sample IDs in the beta matrix appear exactly once in the batch map.  
- Batch labels (``Group``) can be any strings (e.g., ``plate_1``, ``chip_B``), as long as they consistently identify batches.  
- If biological covariates should be adjusted for, handle them upstream before running this script (this wrapper applies basic ComBat only).  

Reference
---------

Johnson, W.E., Li, C., & Rabinovic, A. (2007).  
*Adjusting batch effects in microarray expression data using empirical Bayes methods.*  
Biostatistics, 8(1), 118–127. DOI: see `PubMed 16632515 <https://pubmed.ncbi.nlm.nih.gov/16632515/>`_.
