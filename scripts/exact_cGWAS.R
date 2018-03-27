########### Functions for exact estimations of cGWAS from uGWAS
########### Yakov Tsepilov, Sodbo Sharapov


# The main function. It uses uGWAS results stored in .RData format. 
# CovM - covariance matrix of response vairable and covariates
# all_varg - a numeric vector with variances of studied SNPs
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
	covariates=NULL,
	cn_b = "beta_SNP",
	cn_snp = "SNP",
	cn_se = "se_SNP",
	output_threshold=1e-6,
	good_snps = NULL,
	all_CR = 1,
	correction = TRUE,
	path_uGWAS = NULL){
	
	#checking the data

	cat("Load data...","\n")

	Nprd <- length(covariates)

	#eval(text=paste("snps=",response",[,cn_snp]",sep=""))

	snps=names(all_varg)

	resp_pred<-c(response,covariates)
	
    if (is.null(path_uGWAS)){
        prd <- resp_pred[1]
        for (prd in resp_pred){
            load(paste(path_uGWAS,prd,".RData",sep=""))
        }
	} else {
        prd=resp_pred[1]
        for (prd in resp_pred){
            load(paste(path_uGWAS,prd,".RData",sep=""))
        }
    }
        
	cat("Filtering for NAs and searching overlapping SNPs...","\n")
	prd=resp_pred[1]
	for (prd in resp_pred){
        eval(parse(text=paste("ind=which(is.na(",prd,"[,cn_b]))",sep="")))
        if (length(ind)>0){
            #eval(parse(text=paste(prd,"=",prd,"[!is.na(",prd,"[,cn_b]),]",sep="")))
            eval(parse(text=paste(prd,"=",prd,"[-ind,]",sep="")))
        }
	}
	
	prd=resp_pred[1]
	for (prd in resp_pred){
		eval(parse(text=paste("snps_prd=",prd,"[,cn_snp]",sep="")))
		snps=intersect(snps,snps_prd)
	}
	
	if (!is.null(gut_snps)) snps=intersect(gut_snps,snps)
	#if (!is.null(all_CR)) snps=intersect(names(all_CR),snps)
	
	if (is.null(snps)|length(snps)==0){
		warning("No overlapping SNPs between inputs")
		opt <- options(show.error.messages=FALSE) 
		on.exit(options(opt)) 
		stop() 
	#stop("No overlapping SNPs between inputs")
	}	
	
	snps=unique(snps)
	Nsnps=length(snps)
		
	all_varg=all_varg[snps]
	#all_CR=all_CR[snps]
	prd=resp_pred[1]
	for (prd in resp_pred){
		eval(parse(text=paste("index=match(snps,",prd,"[,cn_snp])",sep="")))
		eval(parse(text=paste(prd,"=",prd,"[index,]",sep="")))
	}
	
	cat(Nsnps,"overlapping SNPs in total...","\n")
	
	if (is.null(covariates)|length(covariates)==0){
	
		cat("No covariates were specified. No cGWAS was performed...","\n")
		
		out=array(NA,c(Nsnps,4))
		colnames(out)=c("b","se","chi2","Pval")
		rownames(out)=snps
		
		prd=resp_pred[1]
		eval(parse(text=paste("out[,\"b\"]<-",prd,"[,cn_b]",sep="")))
		eval(parse(text=paste("out[,\"se\"]<-",prd,"[,cn_se]",sep="")))
		out[,"chi2"]<-(out[,"b"]/out[,"se"])^2
		
		snps=snps[!is.na(out[,"chi2"])]
		out=out[!is.na(out[,"chi2"]),]
		
	} else{
	
		#creating of S matrix
		S=array(NA,c(Nprd+2,Nprd+2))
		colnames(S)=rownames(S)=c(response,covariates,"g")
		S[resp_pred,resp_pred] <- CovM[resp_pred,resp_pred]
		
		#calculation
		out=array(NA,c(Nsnps,4))
		colnames(out)=c("b","se","chi2","Pval")
		rownames(out)=snps
		
		#beta matrix
		betas=array(NA,c(Nsnps,length(resp_pred)))
		rownames(betas)=snps
		colnames(betas)=resp_pred
		
		for (prd in resp_pred){
			eval(parse(text=paste("betas[,prd]<-",prd,"[,cn_b]",sep="")))
		}
		cat("Starting cGWAS...","\n")
		
		i=1
		pb <- txtProgressBar(style=3)
		for (i in (1:Nsnps)){
			setTxtProgressBar(pb, value=i/Nsnps)
			S["g","g"]<-varg<-all_varg[i]
			prd=resp_pred[1]
			for (prd in resp_pred){
				b1=betas[i,prd]
				S["g",prd]<-S[prd,"g"]<-b1*varg
			}
			b<-.b(response=1,pred=2:(Nprd+2),S=S)["g",]
			if (length(all_CR)==1){
                se<-.seb(response=1,pred=2:(Nprd+2),S=S,N=N*all_CR)["g"]
            } else{
                se<-.seb(response=1,pred=2:(Nprd+2),S=S,N=N*all_CR[i])["g"]
            }

			out[i,"b"]<-b
			out[i,"se"]<-se
			#out[i,"chi2"]<-(b/se)^2
		}
		out[,"chi2"]<-(out[,"b"]/out[,"se"])^2
		close(pb)
		
		#out=as.data.frame(out)
		#out=cbind(SNP=snps,out)
		#snps=snps[!is.na(out[,"chi2"])]
		#out=out[!is.na(out[,"chi2"]),]	
		cat("cGWAS was calculated for",dim(out)[1],"SNPs","\n")
	}	
	
	#output
	output=list()
	
	output$gc_lambda=median(out[,"chi2"],na.rm=T)/qchisq(0.5,1,lower.tail=F)
	cat("GC Lambda of all",dim(out)[1],"SNPs is",output$gc_lambda,"\n")
	
	if (correction){
		cat("Correcting...","\n")
		out[,"chi2"]=out[,"chi2"]/output$gc_lambda
		output$status="corrected for GC lambda"
	} else{
		output$status="Not corrected for GC lambda"
	}
	
	out[,"Pval"]=pchisq(out[,"chi2"],1,lower.tail=F)
		
	output$response=response
	output$covariates=paste(covariates,collapse=",")
	
	cat("Filtering results using threshold =",output_threshold,"\n")
	
	snps=snps[out[,"Pval"]<=output_threshold]
	out=out[out[,"Pval"]<=output_threshold,]
	
	cat(dim(out)[1],"SNPs with p-value <=",output_threshold,"\n")
	
	#colnames(out)=c("b","se","chi2","Pval")
	output$results=out
	output$snps=snps
	
	cat("Finished. Enjoy=)","\n")
	return(output)
}


