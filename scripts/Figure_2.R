# Sodbo Sharapov, Yakov Tsepilov (c)

# Create directory for storing a figure

if(!dir.exists('../results'))
  dir.create('../results')

pdf('../results/figure_2.pdf')

# Load table with list of SNPs and metabolites for which we want to create figure 1

assoc_2_plot <- read.table("../data/BN_snps_traits_covariates.txt",
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           sep="\t",
                           fill = T)

# Load SNP infor with var(g)

snp_info <- read.table('../data/30_SNP_information.txt',
                       head = TRUE,
                       stringsAsFactors = FALSE)

var_snp <- snp_info$varg_1785

names(var_snp) <- snp_info$SNP

rm(snp_info)


library(corrplot)
library(ppcor)
par(mfrow=c(1,2))


# Plot 1

