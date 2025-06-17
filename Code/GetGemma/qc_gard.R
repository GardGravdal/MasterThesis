###################################################################################
# This is a script from Kenneth that generates inputGemma file. Make sure that qc.raw
# is generated from qc.R, GWAS$pheno.txt generated from h_dataPrep (at the end of script),
# and that BIM,BED,FAM and PED files are available prior to running.

# At the end of this script you find the system() commands to be ran in command prompt
# for the GEMMA analysis to be run through program downloaded locally.
###################################################################################

library(data.table)

# Set architecture type to run
{
arch <- 10

setwd("C:/Users/gard_/Documents/MasterThesis/Code/GetGemma")
do_qc <- function(fam_file, ## PATH TO .FAM file
                  ncores, # how many cores plink gets to use
                  mem, # how much memory plink gets to use
                  qc_filt, # list of quality control filters
                  plink_path # path to PLINK program
                  ) {
  
  pheno_path = paste0("C:/Users/gard_/Documents/MasterThesis/Code/GetGemma/GWAS_pheno/GWAS_pheno",arch,".txt")
  file_root <- gsub(pattern = ".fam", replacement = "", x = fam_file)
  fam <- fread(fam_file, select = c(1, 2), data.table = FALSE, header = FALSE)
  # Exclude samples with high heterozygosity
  fam_keep <- fam[!grepl(pattern = ".*HIGHHET.*", x = fam$V2), ]
  # Exclude samples with mismatches in sex
  fam_keep <- fam_keep[!grepl(pattern = ".*MISSEX.*", x = fam_keep$V2), ]

  # For inds. genotyped multiple times, keep last one
  fam_keep <- fam_keep[!duplicated(fam_keep$V1, fromLast = TRUE), ]
  
  write.table(fam_keep,
              file = paste0("keep.txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  exit_code <- # system2() calls the command line
    system2(plink_path,
            paste0("--bfile ", file_root, " ",
                   "--recode A ",
                   "--maf ", qc_filt$maf, " ", # Filter by minor allele frequency
                   "--geno ", qc_filt$genorate_snp, " ", # Filter SNPs by call rate
                   "--mind ", qc_filt$genorate_ind, " ", # Filter inds by call rate
                   "--chr-set 32 ", # Sparrow chromosomes
                   "--pheno ", pheno_path, " ",
                   "--pheno-name ", "pheno", " ",
                   "--make-bed", " ",
                   "--memory ", mem, " ",
                   "--keep keep.txt ",
                   "--threads ", ncores, " ",
                   "--out inputGEMMA",arch)) # Name of output file, feel free to change
                                

  
  if (exit_code != 0) {
    stop("Error in plink")
  }
  
  exit_code
}

qc_filt <- list(genorate_ind = 0.05,
                genorate_snp = 0.1,
                maf = 0.01)

do_qc(qc_filt = qc_filt,
      mem = 8000,
      ncores = 6,
      fam_file = "C:/Users/gard_/Documents/MasterThesis/Code/Data/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.fam", # FILL IN PATH TO FAM FILE,
      plink_path = "C:/Users/gard_/Documents/MasterThesis/Code/PLINK/plink_win64_20241013/plink" # FILL IN PATH TO PLINK PROGRAM
        )

}
# This is the process to run the gemma package:
# 1) Start Command prompt
# 2) Change directory (cd path/to/gemma-0.98.5) - Make sure that all the files are there
# 3) Start wsl by executing 'wsl' - Check that all the packages are there with command 'ls'
# 4) Give gemma permission to execute with command 'chmod +x gemma-0.98.5'
# 5) Run the system() commands in Command prompt:
#   i) ./gemma-0.98.5 -bfile inputGEMMA -gk 1 -o RelMat
#   ii) ./gemma-0.98.5 -bfile inputGEMMA -k output/RelMat.cXX.txt -lmm 2 -o GWASresults.lmm

# i)  Computes a genetic relationship matrix (GRM) using the provided genotype data
#-bfile inputGEMMA: Specifies the input binary PLINK file (inputGEMMA.bed, .bim, .fam)
#-gk 1: Computes a kinship matrix (GRM) using the identity-by-state (IBS) method
#-o RelMat: The output files will be named starting with RelMat
# output/RelMat.cXX.txt - The kinship matrix, which is required for mixed model association analysis in step ii)
system("gemma-0.98.5 -bfile inputGEMMA -gk 1 -o RelMat")

# ii) Runs a univariate LMM association analysis for GWAS
# -bfile inputGEMMA: Uses the same genotype dataset as in step i)
# -k output/RelMat.cXX.txt: Uses the kinship matrix computed in step i)
# -lmm 2: Runs the standard linear mixed model (LMM) using the restricted maximum likelihood (REML) method to 
#     correct for population structure and relatedness
# -o GWASresults.lmm: Outputs GWAS association results with this prefix
# output/GWASresults.lmm.assoc.txt - The main GWAS association results file, containing SNP effect sizes, p-values,
#     and other statistical metrics.
system("gemma-0.98.5 -bfile inputGEMMA -k output/RelMat.cXX.txt -lmm 2 -o GWASresults.lmm")

