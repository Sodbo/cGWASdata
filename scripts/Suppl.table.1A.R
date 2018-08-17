# Scritp for reproducing Supplementary Table 1A

# Load table with biochemical distances between 153 metabolites
# and restrict it to a subset of 105 metabolites for which at least the
# one-reaction-step immediate biochemical neighbors are known

bn_dist <- read.table(
  'data/20171207_biochemical_distances.txt',
  header = TRUE,
  row.names = 1
  )

bn_dist[bn_dist != 1] <- 0

# Get list of 105 metabolites

traits <- names(which(apply(bn_dist,2,sum)>0))

rm(bn_dist)

# Dots in the names !!!!!!!!!!!!!!
traits_exists <- sub('.txt','',list.files('data/uGWAS'))

traits <- intersect(traits,traits_exists)

rm(traits_exists)

snp_info <- data.table::fread('zcat data/SNP_information.txt.gz', data.table = FALSE)

minfreq=0.1
minP22=0
minR2=0.3
hw=1e-6
CR=0.95
#
gut_snps <- snp_info$SNP[which(snp_info[,"maf"]*as.numeric(snp_info[,"proper_info"])>=minfreq &
                                      as.numeric(snp_info[,"proper_info"])>=minR2 &
                                      snp_info[,"hw"]>=hw & 
                                      snp_info[,"CR_1785"]>=CR)]

rm(minfreq,minP22,minR2,hw,CR)

source('scripts/clumping_functions.R')

# Create locus table for uGAS

locus_table_ugas <- function_for_making_full_table_without_gcv('data/uGWAS/',
                                                          'SNP',
                                                          'beta',
                                                          'se',
                                                          'P',
                                                          traits, 
                                                          snp_info, 
                                                          thr=5e-8/151, 
                                                          delta = 5e5)

locus_table_ugas <- function_for_shlop_24_10_2013(
  locus_table_ugas,
  p_value="P-value",
  pos="Position",
  snp="SNP",
  delta=5e5,
  chr="Chromosome")

# Create locus table for bnGAS

locus_table_bngas <- function_for_making_full_table_without_gcv('results/BN/', 
                                                               'SNP',
                                                               'b',
                                                               'se',
                                                               'Pval_GC',
                                                               traits, 
                                                               snp_info, 
                                                               thr=5e-8/151, 
                                                               delta = 5e5)

locus_table_bngas <- function_for_shlop_24_10_2013(
  locus_table_bngas,
  p_value="P-value",
  pos="Position",
  snp="SNP",
  delta=5e5,
  chr="Chromosome")

tab_1A <- rbind(locus_table_ugas, locus_table_bngas)

tab_1A <- tab_1A[order(tab_1A$Chromosome, tab_1A$Position),]

tab_1A$code <- paste(tab_1A$SNP,tab_1A$trait,sep = '_')

tab_1A <- tab_1A[!duplicated(tab_1A$code),]

tab_1A <- tab_1A[,c('SNP','trait','Chromosome','Position')]

tab_1A$A1 <- NA
tab_1A$A2 <- NA
tab_1A$freq <- NA
tab_1A$proper_info <- NA

tab_1A$uGAS_b <- NA
tab_1A$uGAS_se <- NA
tab_1A$uGAS_p <- NA
tab_1A$uGAS_p_GC <- NA

tab_1A$bnGAS_b <- NA
tab_1A$bnGAS_se <- NA
tab_1A$bnGAS_p <- NA
tab_1A$bnGAS_p_GC <- NA

tab_1A$bnGAS_noise_comp <- NA
tab_1A$bnGAS_pleiotropic_comp <- NA

tab_1A$cGAS_b <- NA
tab_1A$cGAS_se <- NA
tab_1A$cGAS_p <- NA
tab_1A$cGAS_p_GC <- NA

tab_1A$cGAS_noise_comp  <- NA
tab_1A$cGAS_pleiotropic_comp <- NA

tab_1A$minup_GC  <- NA
tab_1A$minup_t  <- NA
tab_1A$minbdp_GC  <- NA
tab_1A$minbdp_t  <- NA
tab_1A$mincp_GC  <- NA
tab_1A$mincp_t  <- NA
tab_1A$Noise_component_diffrence <- NA
tab_1A$uChi2 <- NA
tab_1A$bdChi2 <- NA
tab_1A$cChi2 <- NA
tab_1A$Ratio_bdChi2_uChi2 <- NA	
tab_1A$Ratio_GGMChi2_uChi2	<- NA
tab_1A$N_bd_cvrts <- NA
tab_1A$bd_cvrts <- NA
tab_1A$N_c_cvrts <- NA
tab_1A$c_cvrts <- NA


# Load lambdas GC for uGWAS
lambda_ugas <- read.table(
  'data/uGWAS_gc_lambda.txt', 
  head = TRUE, 
  stringsAsFactors = FALSE,
  row.names = 1
  )

lambda_ugas <- lambda_ugas[traits,]

names(lambda_ugas) <- traits


# Load lambdas GC for GGM GWAS

lambda_ggm <- read.table(
  'data/GGM_cGWAS_gc_lambda.txt', 
  stringsAsFactors = FALSE,
  head = TRUE,
  row.names = 1
  )

lambda_ggm <- lambda_ggm[traits,]

names(lambda_ggm) <- traits

# Load lambdas GC for BN GWAS

lambda_bn <- read.table(
  'data/BN_cGWAS_gc_lambda.txt', 
  stringsAsFactors = FALSE,
  head = TRUE,
  row.names = 1
  )

lambda_bn <- lambda_bn[traits,]

names(lambda_bn) <- traits

for(index in 1:nrow(tab_1A)){
  
  snp <- tab_1A$SNP[index] 
    
  trait <- tab_1A$trait[index]
  
  tab_1A$A1[index] <- snp_info$A1[snp_info$SNP==snp]
  
  tab_1A$A2[index] <- snp_info$A2[snp_info$SNP==snp]
  
  tab_1A$freq[index] <- snp_info$freq[snp_info$SNP==snp]
  
  tab_1A$proper_info[index] <- snp_info$proper_info[snp_info$SNP==snp]
  
  metab_file_u <- paste0('data/uGWAS_snps_from_paper/',trait,'.txt')
  
  metab_u <- read.table(
    file = metab_file_u,
    head = TRUE,
    stringsAsFactors = FALSE
  )
  
  tab_1A$uGAS_b[index] <- metab_u$b[metab_u$SNP == snp]
  
  tab_1A$uGAS_se[index] <- metab_u$se[metab_u$SNP == snp]
  
  tab_1A$uGAS_p[index] <- pchisq((metab_u$Z[metab_u$SNP == snp])^2,
                                    df = 1,
                                    lower.tail = FALSE
  )
  
  tab_1A$uGAS_p_GC[index] <- pchisq((metab_u$Z[metab_u$SNP == snp])^2/lambda_ugas[trait],
                           df = 1,
                           lower.tail = FALSE
                          )
  
  metab_file_bn <- paste0('results/BN/',trait,'.txt')
  
  metab_bn <- read.table(
    file = metab_file_bn,
    head = TRUE,
    stringsAsFactors = FALSE
  )
  
  tab_1A$bnGAS_b[index] <- metab_bn$b[metab_bn$SNP == snp]
  
  tab_1A$bnGAS_se[index] <- metab_bn$se[metab_bn$SNP == snp]
  
  tab_1A$bnGAS_p[index] <- metab_bn$Pval[metab_bn$SNP == snp]
  
  tab_1A$bnGAS_p_GC[index] <- metab_bn$Pval_GC[metab_bn$SNP == snp]
  
  tab_1A$bnGAS_noise_comp[index]  <- log10((tab_1A$uGAS_se[index]  / tab_1A$bnGAS_se[index] )^2)
  
  tab_1A$bnGAS_pleiotropic_comp[index]  <- log10(metab_bn$chi2[metab_bn$SNP == snp]
                                                 / ((metab_u$Z[metab_u$SNP == snp])^2)
                                                 ) -
                                            tab_1A$bnGAS_noise_comp[index] 
  
  tab_1A$bd_cvrts[index] <- metab_bn$Covariates[metab_bn$SNP == snp]
  
  tab_1A$N_bd_cvrts[index] <- length(unlist(strsplit(unlist(tab_1A$bd_cvrts[index]),split = ';')))
  
  metab_file_ggm <- paste0('results/GGM/',trait,'.txt')
  
  metab_ggm <- read.table(
    file = metab_file_ggm,
    head = TRUE,
    stringsAsFactors = FALSE
  )
  
  tab_1A$cGAS_b[index] <- metab_ggm$b[metab_ggm$SNP == snp]
  
  tab_1A$cGAS_se[index] <- metab_ggm$se[metab_ggm$SNP == snp]
  
  tab_1A$cGAS_p[index] <- metab_ggm$Pval[metab_ggm$SNP == snp]
  
  tab_1A$cGAS_p_GC[index] <- metab_ggm$Pval_GC[metab_ggm$SNP == snp]
  
  tab_1A$cGAS_noise_comp[index]  <- log10((tab_1A$uGAS_se[index]  / tab_1A$cGAS_se[index] )^2)
  
  tab_1A$cGAS_pleiotropic_comp[index]  <- log10(metab_ggm$chi2[metab_ggm$SNP == snp] 
                                                 / ((metab_u$Z[metab_u$SNP == snp])^2)
  ) -
    tab_1A$cGAS_noise_comp[index]
  
  tab_1A$c_cvrts[index] <- metab_ggm$Covariates[metab_ggm$SNP == snp]
  
  tab_1A$N_c_cvrts[index] <- length(unlist(strsplit(unlist(tab_1A$c_cvrts[index]),split = ';')))
  
  # Find min uGWAS P and trait
  
  snp_u_p <- system(paste0('grep -w ', snp, ' data/uGWAS_snps_from_paper/*'),intern = TRUE)
  
  snp_u_p <- data.frame(matrix(unlist(sapply(snp_u_p,strsplit,split='\t')),byrow=TRUE,nrow=length(snp_u_p)), stringsAsFactors = FALSE)
  
  colnames(snp_u_p) <- c('trait','b','se','Z','P')
  
  snp_u_p$P <- as.numeric(snp_u_p$P)
  
  snp_u_p$Z <- as.numeric(snp_u_p$Z)
  
  snp_u_p$trait <- sub('data/uGWAS_snps_from_paper/','', snp_u_p$trait) 
  
  snp_u_p$trait <- sub(paste0('.txt:',snp),'', snp_u_p$trait)
  
  snp_u_p <- snp_u_p[snp_u_p$trait %in% traits,]
  
  snp_u_p <- snp_u_p[which.min(snp_u_p$P),,drop=FALSE]
  
  tab_1A$minup_GC[index] <- pchisq(snp_u_p$Z^2/lambda_ugas[snp_u_p$trait],df=1,lower.tail = FALSE)
    
  tab_1A$minup_t[index] <- snp_u_p$trait 
  
  # Find min bnGWAS P and trait
  
  snp_bn_p <- system(paste0('grep -w ', snp, ' results/BN/*'),intern = TRUE)
  
  snp_bn_p <- data.frame(matrix(unlist(sapply(snp_bn_p,strsplit,split='\t')),byrow=TRUE,nrow=length(snp_bn_p)), stringsAsFactors = FALSE)
  
  colnames(snp_bn_p) <- c('trait','b','se','Chi2','P')
  
  snp_bn_p$P <- as.numeric(snp_bn_p$P)
  
  snp_bn_p$Chi2 <- as.numeric(snp_bn_p$Chi2)
  
  snp_bn_p$trait <- sub('results/BN/','', snp_bn_p$trait) 
  
  snp_bn_p$trait <- sub(paste0('.txt:',snp),'', snp_bn_p$trait)
  
  snp_bn_p <- snp_bn_p[snp_bn_p$trait %in% traits,]
  
  snp_bn_p <- snp_bn_p[which.min(snp_bn_p$P),,drop=FALSE]
  
  tab_1A$minbdp_GC[index] <- pchisq(snp_bn_p$Chi2/lambda_bn[snp_bn_p$trait],df=1,lower.tail = FALSE)
  
  tab_1A$minbdp_t[index] <- snp_bn_p$trait
  
  # Find min cGWAS P and trait
  
  snp_ggm_p <- system(paste0('grep -w ', snp, ' results/GGM/*'),intern = TRUE)
  
  snp_ggm_p <- data.frame(matrix(unlist(sapply(snp_ggm_p,strsplit,split='\t')),byrow=TRUE,nrow=length(snp_ggm_p)), stringsAsFactors = FALSE)
  
  colnames(snp_ggm_p) <- c('trait','b','se','Chi2','P')
  
  snp_ggm_p$P <- as.numeric(snp_ggm_p$P)
  
  snp_ggm_p$Chi2 <- as.numeric(snp_ggm_p$Chi2)
  
  snp_ggm_p$trait <- sub('results/GGM/','', snp_ggm_p$trait) 
  
  snp_ggm_p$trait <- sub(paste0('.txt:',snp),'', snp_ggm_p$trait)
  
  snp_ggm_p <- snp_ggm_p[snp_ggm_p$trait %in% traits,]
  
  snp_ggm_p <- snp_ggm_p[which.min(snp_ggm_p$P),,drop=FALSE]
  
  tab_1A$mincp_GC[index] <- pchisq(snp_ggm_p$Chi2/lambda_ggm[snp_ggm_p$trait],df=1,lower.tail = FALSE)
  
  tab_1A$mincp_t[index] <- snp_ggm_p$trait
  
  tab_1A$Noise_component_diffrence[index] <- tab_1A$cGAS_noise_comp[index] - tab_1A$bnGAS_noise_comp[index]
  
  tab_1A$uChi2[index] <- qchisq(tab_1A$minup_GC[index], df=1, lower.tail=FALSE)
  
  tab_1A$bdChi2[index] <- qchisq(tab_1A$minbdp_GC[index], df=1, lower.tail=FALSE)
  
  tab_1A$cChi2[index] <- qchisq(tab_1A$mincp_GC[index], df=1, lower.tail=FALSE)
  
}

tab_1A$locus <- 1

for(index in 2:nrow(tab_1A)){
  
   ifelse(tab_1A$Chromosome[index-1]==tab_1A$Chromosome[index] & abs(tab_1A$Position[index-1] - tab_1A$Position[index])<5e5, tab_1A$locus[index] <- tab_1A$locus[index-1], tab_1A$locus[index] <- tab_1A$locus[index-1]+1)
}

for(locus in unique(tab_1A$locus)) {
  
  if(sum(tab_1A$locus==locus)==2){
  
    tab_1A[tab_1A$locus==locus,c('uChi2','bdChi2','cChi2')][1,] <- pmax(tab_1A[tab_1A$locus==locus,c('uChi2','bdChi2','cChi2')][1,],tab_1A[tab_1A$locus==locus,c('uChi2','bdChi2','cChi2')][2,])
    
    tab_1A[tab_1A$locus==locus,c('uChi2','bdChi2','cChi2')][2,] <- NA
    
  }
  
  tab_1A[tab_1A$locus==locus,c('Ratio_bdChi2_uChi2')][1] <- tab_1A[tab_1A$locus==locus,c('bdChi2')][1] / tab_1A[tab_1A$locus==locus,c('uChi2')][1]
  
  tab_1A[tab_1A$locus==locus,c('Ratio_GGMChi2_uChi2')][1] <- tab_1A[tab_1A$locus==locus,c('cChi2')][1] / tab_1A[tab_1A$locus==locus,c('uChi2')][1]

}

data.table::fwrite(tab_1A,
  sep = '\t',
  file = 'results/Suppl_table_1A.tsv'
)

# Remove locus on chromosome 2, See supplementary Note 2 for the fetails

tab_1A <- tab_1A[-c(3:4),]

# The average ratio of the maximum test statistic between BN-cGAS and uGAS
print(mean(tab_1A$Ratio_bdChi2_uChi2,na.rm = TRUE))

# A paired-sample Wilcoxon test
print(wilcox.test(tab_1A$uChi2-tab_1A$bdChi2)$p.value)