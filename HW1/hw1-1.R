#log determinant
logdet <- function(R)
{
	return(sum(log(eigen(R)$values)))	
}

#log marginal likelihood
logmarglik <- function(data,A)
{
	#the response is the first column 
	D1 = data[,1];	
	#get the colums associated with the explanatory variables A
	DA = data[,A];
	
	MA = diag(length(A)) + t(DA) %*% DA;
	return(lgamma((nrow(DA)+ncol(DA)+2)/2)-lgamma((ncol(DA)+2)/2) - 0.5*logdet(MA) - 0.5*(nrow(DA)+ncol(DA)+2)*log(1+ t(D1)%*%D1 - t(D1)%*%DA%*%solve(MA)%*%t(DA)%*%D1))
}

#test to see that the function does
data = as.matrix(read.table("erdata.txt"));
A = c(2,5,10);

logmarglik(data, A)