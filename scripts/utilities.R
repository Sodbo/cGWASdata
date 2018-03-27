#core functions

#betas
	
.b=function(response,pred,S){
	out <- solve(S[pred,pred])%*%S[response,pred]
	return(out)
}

# standart errors of betas

.seb=function(response,pred,S,N){
	sigma_joint=(S[response,response]-t(.b(response,pred,S=S))%*%S[response,pred])
	#sigma_joint=S[response,response]
	out=as.vector(sigma_joint)*solve(S[pred,pred])/(N-length(pred)-1)
	
	out=diag(out)
	out=sqrt(out)
	return(out)
}

# 
.manova=function(response,pred,S,N){
	sigma_joint=(S[response,response]-t(.b(response,pred,S=S))%*%S[response,pred])
	Fst=((S[response,response]-as.vector(sigma_joint))/as.vector(sigma_joint))*((N-length(pred)-1)/length(pred))
	return(Fst)
}

#T statistics (not T^2)

.t=function(response,pred,S,N){
	out=.b(response=response,pred=pred,S=S)/.seb(response=response,pred=pred,S=S,N=N)
	return(out)
}

ReadCovariates=function(covariates,path2data,suffix=".txt"){
	Nprd=length(covariates)
	prd=covariates[1]
	for (prd in covariates){
	if (!exists(prd)){
		eval(text=paste(prd,"<-read.table(paste(path2data,prd,suffix,sep=\"\"),header=T,stringsAsFactors=F)",sep=""))
	}
	}
}
