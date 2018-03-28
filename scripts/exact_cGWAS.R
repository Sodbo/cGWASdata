########### Functions for exact estimations of cGWAS from uGWAS
########### Yakov Tsepilov, Sodbo Sharapov


# The main function. It uses uGWAS results stored in .RData format. 
# CovM - covariance matrix of response vairable and covariates
# all_varg - a named numeric vector with variances of studied SNPs
# response - a character, name of response variable
# N - sample size
# covariates - a character vector with list of covariates
# cn_* - column names
# output_threshold - the threshold for p-value filter for storing of outputs
# good_snps - the list of SNPs that should be used in anlisys (in case you want to filter them)
# all_CR - call rate for each SNP
# correction - should be GC correction applied?
# path_uGWAS - path to uGWAS in .RData format

exact_cGWAS <- function(
	CovM,
	all_varg,
	response,
	N,
	covariates = NULL,
	cn_b = "beta_SNP",
	cn_snp = "SNP",
	cn_se = "se_SNP",
	output_threshold = 1e-6,
	good_snps = NULL,
	all_CR = 1,
	correction = FALSE,
	path_uGWAS = './'){
	
	#checking the data

	cat("Loading data ...","\n")

	# In this list sum stats for responce and covariates will be collected

	sum_stats <- list()

	# Reading sum_stats and pulling into the list

	for(trait in c(response,covariates)) {

		trait_file <- paste0(path_uGWAS,'/',trait, '.txt')

		sum_stats[[length(sum_stats)+1]] <- read.table(trait_file, 
			head = TRUE, 
			stringsAsFactors = FALSE)

		rm(trait)

	}

	names(sum_stats) <- c(response, covariates)

	all_traits <- names(sum_stats)

	cat("Filtering for NAs ...","\n")

	for(trait in all_traits){

		sum_stats[[trait]] <- na.omit(sum_stats[[trait]])

		rm(trait)

	}
	  
	cat("Searching for overlapping SNPs ...","\n")

	snps_over <- sum_stats[[1]][,cn_snp]

	for(trait in all_traits) {

		snps_over <- intersect(snps_over, sum_stats[[trait]][,cn_snp])

		rm(trait)

	}

	# Intersect with good SNPs

	if (!is.null(good_snps)) 

		snps_over <- intersect(good_snps, snps_over)

	if (is.null(snps_over) | length(snps_over)==0){

		warning("No overlapping SNPs between inputs")
		opt <- options(show.error.messages=FALSE) 
		on.exit(options(opt)) 
		stop() 

	}	
	
	snps <- unique(snps_over)

	rm(snps_over)

	all_varg <- all_varg[snps]

	for(trait in all_traits){

		keep_snps <- match(sum_stats[[trait]][,cn_snp],snps, nomatch = FALSE)

		sum_stats[[trait]] <- sum_stats[[trait]][keep_snps,]

		rm(trait)

	}

	Nsnps <- length(snps)

	cat(Nsnps,"overlapping SNPs in total...","\n")
	
	if (is.null(covariates) | length(covariates) == 0){
	
		cat("No covariates were specified. No cGWAS was performed ...","\n")
		cat("Returnign uGWAS sum_stats for the ", names(sum_stats)[1],"\n")
		
		out <- data.frame(
			b = sum_stats[[1]][,cn_b], 
			se = sum_stats[[1]][,cn_se],
			chi2 = (sum_stats[[1]][,cn_b] / sum_stats[[1]][,cn_se])^2
			)
		
		rownames(out) <- snps

		# Pval is calculated in the very end after GC
		
	} else {
	
		# Creating S matrix

		source('../code/utilities.R')

		Nprd <- length(all_traits)

		S <- array(NA, c(Nprd + 1, Nprd + 1))

		colnames(S) <- rownames(S) <- c(names(sum_stats), "g")

		S[all_traits, all_traits] <- as.matrix(CovM[all_traits, all_traits])
	
		# calculation

		out <- array(NA,c(Nsnps, 4))

		colnames(out) <- c("b","se","chi2","Pval")

		rownames(out) <- snps
		
		# beta matrix
		betas <- array(NA, c(Nsnps, length(all_traits)))
		rownames(betas) <- snps
		colnames(betas) <- all_traits
		
		for (prd in all_traits){

			betas[,prd] <- sum_stats[[prd]][,cn_b]

		}
		cat("Starting cGWAS...","\n")
		
		#pb <- txtProgressBar(style=3)

		for (i in 1:Nsnps){

			#setTxtProgressBar(pb, value = i/Nsnps)

			S["g","g"] <- varg <- all_varg[i]

			for (prd in all_traits){

				b1 <- betas[i,prd]

				S["g",prd] <- S[prd,"g"] <- b1*varg

			}

			out[i,"b"] <-.b(response=1, pred=2:(Nprd+1),S=S)["g",]

            out[i,"se"] <-.seb(response=1, pred=2:(Nprd+1),S = S, N = N*all_CR[i])["g"]

		}

		out[,"chi2"] <- (out[,"b"]/out[,"se"])^2

		out <- as.data.frame(out)

		#close(pb)
		
		#out=as.data.frame(out)
		#out=cbind(SNP=snps,out)
		#snps=snps[!is.na(out[,"chi2"])]
		#out=out[!is.na(out[,"chi2"]),]	
		cat("cGWAS was calculated for",dim(out)[1],"SNPs","\n")
	}	
	
	# Combinig output 

	output <- list()
	
	output$gc_lambda <- median(out$chi2,na.rm=T) / qchisq(0.5, 1, lower.tail=F)
	
	cat("GC Lambda of all",dim(out)[1],"SNPs is",output$gc_lambda,"\n")
	
	if (correction){
		cat("Performing genomic control correction ...","\n")
		
		out$chi2 <- out$chi2 / output$gc_lambda
		output$status="Corrected for GC lambda"

	} else{

		output$status <- "Not corrected for GC lambda"

	}
	
	out$Pval <- pchisq(out$chi2, 1, lower.tail = FALSE)
		
	output$response <- response

	output$covariates <- paste(covariates, collapse=";")
	
	cat("Filtering results using threshold = ",output_threshold,"\n")

	out <- out[out$Pval <=output_threshold ,]
	
	cat(dim(out)[1],"SNPs with p-value <=",output_threshold,"\n")
	
	#colnames(out)=c("b","se","chi2","Pval")
	
	output$results <- out
	
	output$snps <- rownames(out)
	
	cat("Finished. Enjoy=)","\n")
	
	return(output)

}