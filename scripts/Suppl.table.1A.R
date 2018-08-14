# Scritp for reproducing Supplementary Table 1A

<<<<<<< HEAD
=======
Pval_thr <- 5e-8/105

>>>>>>> 8459dbb901c1d8619fedb490813ee2ba031febcf
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

<<<<<<< HEAD
rm(bn_dist)

# Dots in the names !!!!!!!!!!!!!!
=======
>>>>>>> 8459dbb901c1d8619fedb490813ee2ba031febcf
traits_exists <- sub('.txt','',list.files('data/uGWAS'))

traits <- intersect(traits,traits_exists)

rm(traits_exists)

<<<<<<< HEAD
snp_info <- data.table::fread('zcat data/SNP_information.txt.gz', data.table = FALSE)
=======
snp_info <- data.table::fread('zcat data/SNP_information.txt.zip', data.table = FALSE)

snp_info$maf <- pmin(1 - snp_info$freq, snp_info$freq)
>>>>>>> 8459dbb901c1d8619fedb490813ee2ba031febcf

minfreq=0.1
minP22=0
minR2=0.3
hw=1e-6
CR=0.95
#
<<<<<<< HEAD
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
=======
gut_snps <- snp_info$SNP[which(snp_info[,"maf"]*as.numeric(snp_info[,"R2_impute_info"])>=minfreq &
                                      as.numeric(snp_info[,"R2_impute_info"])>=minR2 &
                                      snp_info[,"hw"]>=hw)]

source('scripts/clumping_functions.R')

locus_table <- function_for_making_full_table_without_gcv('data/uGWAS/', 
                                                          traits, 
                                                          snp_info, 
                                                          thr=5e-8, 
                                                          delta = 5e5)

locus_table_2 <- function_for_shlop_24_10_2013(
  locus_table,
>>>>>>> 8459dbb901c1d8619fedb490813ee2ba031febcf
  p_value="P-value",
  pos="Position",
  snp="SNP",
  delta=5e5,
  chr="Chromosome")

<<<<<<< HEAD
tab_1A <- rbind(locus_table_ugas, locus_table_bngas)

tab_1A <- tab_1A[order(tab_1A$Chromosome, tab_1A$Position),]

tab_1A$code <- paste(tab_1A$SNP,tab_1A$trait,sep = '_')

tab_1A <- tab_1A[!duplicated(tab_1A$code),]

tab_1A <- tab_1A[,c('SNP','trait','Chromosome','Position')]

# Remove locus on chromosome 2, See supplementary Note 2 for the fetails

tab_1A <- tab_1A[-c(3:4),]

tab_1A$A1 <- NA
tab_1A$A2 <- NA
tab_1A$freq <- NA
tab_1A$proper_info <- NA

tab_1A$uGAS_b <- NA
tab_1A$uGAS_se <- NA
tab_1A$uGAS_p <- NA

tab_1A$bnGAS_b <- NA
tab_1A$bnGAS_se <- NA
tab_1A$bnGAS_p <- NA

tab_1A$bnGAS_noise_comp	<- NA
tab_1A$bnGAS_pleiotropic_comp <- NA

tab_1A$cGAS_b <- NA
tab_1A$cGAS_se <- NA
tab_1A$cGAS_p <- NA

tab_1A$cGAS_noise_comp	<- NA
tab_1A$cGAS_pleiotropic_comp <- NA

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
  
  tab_1A$uGAS_p[index] <- pchisq((metab_u$Z[metab_u$SNP == snp])^2/lambda_ugas[trait],
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
  
  tab_1A$bnGAS_p[index] <- metab_bn$chi2[metab_bn$SNP == snp] / lambda_bn[trait]
  
  tab_1A$bnGAS_p[index] <- pchisq(tab_1A$bnGAS_p[index], df=1 , lower.tail = FALSE)
  
  tab_1A$bnGAS_noise_comp[index]  <- log10((tab_1A$uGAS_se[index]  / tab_1A$bnGAS_se[index] )^2)
  
  tab_1A$bnGAS_pleiotropic_comp[index]  <- log10(metab_bn$chi2[metab_bn$SNP == snp] / lambda_bn[trait] 
                                                 / ((metab_u$Z[metab_u$SNP == snp])^2/lambda_ugas[trait])
                                                 ) -
                                            tab_1A$bnGAS_noise_comp[index] 
  
  metab_file_ggm <- paste0('results/GGM/',trait,'.txt')
  
  metab_ggm <- read.table(
    file = metab_file_ggm,
    head = TRUE,
    stringsAsFactors = FALSE
  )
  
  tab_1A$cGAS_b[index] <- metab_ggm$b[metab_ggm$SNP == snp]
  
  tab_1A$cGAS_se[index] <- metab_ggm$se[metab_ggm$SNP == snp]
  
  tab_1A$cGAS_p[index] <- metab_ggm$chi2[metab_ggm$SNP == snp] / lambda_ggm[trait]
  
  tab_1A$cGAS_p[index] <- pchisq(tab_1A$cGAS_p[index], df=1 , lower.tail = FALSE)
  
  tab_1A$cGAS_noise_comp[index]  <- log10((tab_1A$uGAS_se[index]  / tab_1A$cGAS_se[index] )^2)
  
  tab_1A$cGAS_pleiotropic_comp[index]  <- log10(metab_ggm$chi2[metab_ggm$SNP == snp] / lambda_ggm[trait] 
                                                 / ((metab_u$Z[metab_u$SNP == snp])^2/lambda_ugas[trait])
  ) -
    tab_1A$cGAS_noise_comp[index] 
  
  
}

tab_1A[,c('bnGAS_noise_comp','bnGAS_pleiotropic_comp','cGAS_noise_comp','cGAS_pleiotropic_comp')]

=======
locus_table_2[locus_table_2$`P-value`< 5e-8/151,]

	"metab_file_bn <- paste0('results/BN/',trait,'.txt')
	metab_bn <- read.table(
		file = metab_file_bn,
		head = TRUE,
		stringsAsFactors = FALSE
		)

	metab_bn <- metab_bn[metab_bn$Pval < Pval_thr,]

	locus_tab_bn <- rbind(locus_tab_bn,metab_bn) "

	"tab1$uGAS_b[i] <- metab_u$beta[metab_u$SNP == snp]

	tab1$uGAS_se[i] <- metab_u$se[metab_u$SNP == snp]

	tab1$uGAS_p[i] <- pchisq((metab_u$Z[metab_u$SNP == snp])^2/lambda_ugas$gc_lambda[lambda_ugas$trait == trait],
		df = 1,
		lower.tail = FALSE
		)

	tab1$bdGAS_b[i] <- metab_bn$b[metab_bn$SNP == snp]

	tab1$bdGAS_se[i] <- metab_bn$se[metab_bn$SNP == snp]

	tab1$bdGAS_p[i] <- metab_bn$chi2[metab_bn$SNP == snp] / lambda_bn[trait]
	
	tab1$bdGAS_p[i] <- pchisq(tab1$bdGAS_p[i], df=1 , lower.tail = FALSE)

	tab1$bd_noise_comp <- log10((tab1$uGAS_se / tab1$bdGAS_se)^2)
	
	tab1$cGAS_b[i] <- metab_ggm$b[metab_ggm$SNP == snp]

	tab1$cGAS_se[i] <- metab_ggm$se[metab_ggm$SNP == snp]

	tab1$cGAS_p[i] <- metab_ggm$chi2[metab_ggm$SNP == snp] / lambda_ggm[trait]

	tab1$cGAS_p[i] <- pchisq(tab1$cGAS_p[i], df=1 , lower.tail = FALSE)"
>>>>>>> 8459dbb901c1d8619fedb490813ee2ba031febcf

"write.table(tab1, 
	quote = FALSE,
	row.names = FALSE,
	sep = '\t',
#	dec = ',',
	file = 'results/Suppl_table_1A.csv'
)"


"# Load lambdas GC for uGWAS
lambda_ugas <- read.table(
	'data/uGWAS_gc_lambda.txt', 
	head = TRUE, 
	stringsAsFactors = FALSE,
	row.names = 1
	)


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

names(lambda_bn) <- traits"