# Sodbo Sharapov, Yakov Tsepilov (c)

# Create directory for storing a figure

if(!dir.exists('results'))
  dir.create('results')

# Load table with list of SNPs and metabolites for which we want to create figure 1

assoc_2_plot <- read.table("data/SNPs_for_suppl_tab1.txt",
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           #sep="\t",
                           fill = T)

# Load SNP infor with var(g)

snp_info <- read.table('data/30_SNP_information.txt',
                       head = TRUE,
                       stringsAsFactors = FALSE)

var_snp <- snp_info$varg_1785

names(var_snp) <- snp_info$SNP

rm(snp_info)

# Load table with correlation matrix for phenotypes

cor_matrix <- read.table('data/20171207_corr_matrix.txt')

xv <- NULL

assoc_2_plot$noise <- NA

assoc_2_plot$pleiotropic <- NA

for (i in 1:nrow(assoc_2_plot)){
  
  # Assign trait and SNP into variabels
  snp <- assoc_2_plot$SNP[i]
  
  trait <- assoc_2_plot$trait [i]
  
  # Load file with uGAS statistics
  ugas_file <- paste0('data/uGWAS_snps_from_paper/',trait,'.txt')
  
  ugas_stats <- read.table(ugas_file,
                           head = TRUE,
                           stringsAsFactors = FALSE)
  
  # Get SE and Chi2 from uGAS stats
  
  beta_snp_u <- ugas_stats$b[ugas_stats$SNP == snp]
  
  se_snp_u <- ugas_stats$se[ugas_stats$SNP == snp]
  
  chi2_snp_u <- (ugas_stats$Z[ugas_stats$SNP == snp])^2
  
  rm(ugas_file, ugas_stats)
  
  # Load file with biochemical network cGAS statistics
  bngas_file <- paste0('results/BN/',trait,'.txt')
  
  bngas_stats <- read.table(bngas_file,
                            head = TRUE,
                            stringsAsFactors = FALSE)
  
  # Get SE and Chi2 from biochemical network cGAS statistics
  
  se_snp_c <- bngas_stats$se[bngas_stats$SNP == snp]
  
  chi2_snp_c <- bngas_stats$chi2[bngas_stats$SNP == snp]
  
  
  # Get betas for SNP and covariates
  
  beta_snp_c <- bngas_stats$b[bngas_stats$SNP == snp]
  
  if(is.na(bngas_stats$b_cov[bngas_stats$SNP == snp]))
    next
  
  if(is.numeric(bngas_stats$b_cov[bngas_stats$SNP == snp])){
    
    beta_cov_c <- bngas_stats$b_cov[bngas_stats$SNP == snp]
    
  } else{
    
    beta_cov_c <- unlist(strsplit(bngas_stats$b_cov[bngas_stats$SNP == snp],";"))
    
    beta_cov_c <- as.numeric(beta_cov_c)
    
  }
  
  # Assign list of covariates
  
  covariates <- unlist(strsplit(bngas_stats$Covariates[i],";"))
  
  rm(bngas_file, bngas_stats)
  
  beta_cov_u <- NULL
  
  for(Covariate in covariates) {
    
    cov_file <- paste0('data/uGWAS_snps_from_paper/',Covariate,'.txt')
    
    cov_stats <- read.table(cov_file,
                            head = TRUE,
                            stringsAsFactors = FALSE)
    
    beta_cov_u <- c(beta_cov_u, cov_stats$beta[cov_stats$SNP == snp])
    
  }
  
  cor_trait_snp  <- beta_snp_u * sqrt(var_snp[snp])
  
  cor_cov_snp <- beta_cov_u * sqrt(var_snp[snp])
  
  rm(cov_file, Covariate, cov_stats)
  
  # Get correlation between trait and covariates
  cor_trait_cov <- as.numeric(cor_matrix[trait, covariates])
  
  f1 <- log10( (se_snp_u / se_snp_c )^2 )
  f2 <- log10( ( 1 - ( beta_cov_c %*% cor_cov_snp ) / cor_trait_snp )^2 )
  
  xv <- c(xv, log10 (chi2_snp_c / chi2_snp_u) )

  assoc_2_plot$noise[i] <-  f1
  
  assoc_2_plot$pleiotropic[i] <-  f2

}

write.table(assoc_2_plot[,c('noise','pleiotropic')],
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            dec = ',')