BASEDIR = "#baseDir"
RESULTS_DIR = "#resultsDir"

##########################################################
# CHANGE TO PATH OF xhmm INSTALLATION:
##########################################################
XHMM_PATH = BASEDIR

XHMM_SCRIPTS_PATH = paste(XHMM_PATH, "/R_functions", sep="")
source(paste(XHMM_SCRIPTS_PATH, "/sourceDir.R", sep=""))
sourceDir(XHMM_SCRIPTS_PATH, trace = FALSE, recursive = TRUE)


##########################################################
# CHANGE TO DESIRED OUTPUT PATH FOR PLOTS:
##########################################################
PLOT_PATH = paste(RESULTS_DIR, "/plots", sep="")


##########################################################
# CHANGE TO FULL PATH AND PREFIX FOR ALL XHMM PIPELINE OUTPUT:
#
# (e.g., DATA.RD.txt, DATA.RD_PCA.PC.txt, 
# DATA.PCA_normalized.filtered.sample_zscores.RD.txt, etc.)
#
##########################################################
JOB_PREFICES = paste(RESULTS_DIR, "/DATA", sep="")


##########################################################
# CHANGE TO FULL PATH OF PLINK/Seq FORMAT EXON-TO-GENE LIST:
#
# THIS CAN BE GENERATED FOR EXOME.interval_list USING THE FOLLOWING PLINK/Seq COMMAND:
#
# pseq . loc-intersect --group refseq --locdb /path/to/locdb --file ./EXOME.interval_list --out ./annotated_targets.refseq
#
##########################################################
JOB_TARGETS_TO_GENES = "/groupvol/cgi/resources/annotations/Agilent_SureSelect_Human_All_Exon_V4+UTRs.targets.geneIDs.int"



##########################################################
# CAN BE (OPTIONALLY) POPULATED FROM EITHER CHOICE BELOW:
##########################################################
SAMPLE_FEATURES = NULL


##########################################################
# CHOICE A: CHANGE TO FULL PATH OF PLINK/Seq PEDIGREE FILE:
##########################################################
#PEDIGREE_FILE = "/path/to/pseq_pedigree_file"
#PEDIGREE_DATA = readPedigreeFile(PEDIGREE_FILE)
#SAMPLE_FEATURES = pedigreeDataToBinarySampleProperties(PEDIGREE_DATA)


##########################################################
# CHOICE B: OR, INSTEAD OF A PLINK/Seq PEDIGREE FILE, USE A PLINK/Seq SAMPLE PHENOTYPES FILE:
##########################################################
#PHENOTYPES_FILE = "/path/to/pseq_phenotypes_file"
#PHENOTYPES_DATA = readPhenotypesFile(PHENOTYPES_FILE)
#SAMPLE_FEATURES = phenotypeDataToBinarySampleProperties(PHENOTYPES_DATA)



XHMM_plots(PLOT_PATH, JOB_PREFICES, XHMM_PATH, JOB_TARGETS_TO_GENES, SAMPLE_FEATURES)
