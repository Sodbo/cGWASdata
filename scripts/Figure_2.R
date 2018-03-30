# Sodbo Sharapov, Yakov Tsepilov (c)

# Create directory for storing a figure

if(!dir.exists('../results'))
  dir.create('../results')

png('../results/figure_2.png',
    height = 720,
    width = 1080)

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

# Load table with correlation matrix for phenotypes

cor_matrix <- read.table('../data/20171207_corr_matrix.txt')

# Create M_1 matrix for FADS1

create_M <- function(assoc_2_plot, snp, gene){
    
    index <- which(assoc_2_plot$SNP == snp & assoc_2_plot$Gene == gene)
    
    trait <- assoc_2_plot$Metabolite[index]
    
    covariates <- unlist(strsplit(assoc_2_plot$Covariates[index], split = ';'))
    
    M <- matrix(NA, nrow = 2 + length(covariates), ncol = 2 + length(covariates))
    
    rownames(M) <- colnames(M) <- c(trait,'SNP',covariates)
    
    # Fill correlations between trait and covariates
    
    M [1, covariates] <- as.numeric(cor_matrix[trait, covariates])
    
    M [covariates , covariates] <- as.matrix(cor_matrix[covariates, covariates])
    
    diag(M) <- 0
    
    M[lower.tri(M)] <- 0
    
    # Fill partial correlations
    
    bn_gas_file <- paste0('../results/BN/', trait, '.txt')
    
    bn_gas <- read.table(file = bn_gas_file,
                         head = TRUE,
                         stringsAsFactors = FALSE)
    
    betas <- c (bn_gas$b[bn_gas$SNP == snp])
    
    names(betas) <- 'SNP'
    
    b_cov <- bn_gas$b_cov[bn_gas$SNP == snp]
    
    if(!is.numeric(b_cov)){
      
      b_cov <- unlist(strsplit(bn_gas$b_cov[bn_gas$SNP == snp], split = ';'))
      
      b_cov <- as.numeric(b_cov)
      
    }
    
    names(b_cov) <- unlist(strsplit(bn_gas$Covariates[bn_gas$SNP == snp], split = ';'))
    
    betas <- c(betas, b_cov)
    
    M[names(betas),trait] <- betas
    
    # Fill correlations between SNP and trait
    
    u_gas_file <- paste0('../data/uGWAS_snps_from_paper/', trait, '.txt')
    
    u_gas <- read.table(file = u_gas_file,
                         head = TRUE,
                         stringsAsFactors = FALSE) 
    
    cor_trait_snp <- u_gas$beta[u_gas$SNP == snp] * sqrt(var_snp[snp])
    
    M[trait,'SNP'] <- cor_trait_snp
    
    # Fill correlations between SNP and covariates
      
    cor_cov_snp <- NULL
    
    for(covariate in covariates){
      
      u_gas_file <- paste0('../data/uGWAS_snps_from_paper/', covariate, '.txt')
      
      u_gas <- read.table(file = u_gas_file,
                          head = TRUE,
                          stringsAsFactors = FALSE) 
      
      cor_cov_snp <- c(cor_cov_snp, u_gas$beta[u_gas$SNP == snp] * sqrt(var_snp[snp]))
      
    }
    
    M['SNP',covariates] <- cor_cov_snp

    return(M)    
}
                                    
# Load libraries

library(corrplot)
library(ppcor)
par(mfrow=c(1,2))

M_1 <- create_M(assoc_2_plot,'rs174547','FADS1')

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M_1,cl.pos="b",na.label="X",method = "square",
         addCoef.col = "black",col = col(200),number.cex = 1, tl.cex=2,tl.col=c("blue","darkblue","red"))

M_2 <- create_M(assoc_2_plot,'rs8396','ETFDH')

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M_2,cl.pos="b",na.label="X",method = "square",
         addCoef.col = "black",col = col(200),number.cex = 1, tl.cex=2,tl.col=c("blue","darkblue","red","red"))

dev.off()
