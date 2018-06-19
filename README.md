# Hematopoiesis2018

This repository contains all the needed code and metadata files to cluster and generate figures for Giladi A. & Paul F. et al. Single-cell characterization of haematopoietic progenitors and their trajectories in homeostasis and perturbed haematopoiesis, Nature Cell Biology 2018

In order to run the scripts, downloaded processed data from the GSE113495 needs to be added in specific folders:
Processed UMI-tab files (ABXXX.txt) from GSE92575 should be copied to the folder output/umi.tab
Processed CRISP-seq data from GSE113494	should be copied to annotations/crispseq_count.txt

To start analysis, run from the root directory:
Rscript run.r

The following R libraries should be installed:

Please send questions to Amir Giladi: amir.goldberg@weizmann.ac.il