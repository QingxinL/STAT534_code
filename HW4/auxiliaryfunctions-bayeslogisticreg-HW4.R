#returns the MLEs of the coefficients of a logistic regression
getcoefglm <- function(response,explanatory,data)
{
  #check if the regression has no explanatory variables
  if(0==length(explanatory))
  {
     #regression with no explanatory variables
     return(coef(glm(data[,response] ~ 1,family=binomial(link=logit))));
  }  
  return(coef(glm(data[,response] ~ data[,explanatory],family=binomial(link=logit))));
}

#log determinant
logdet <- function(R)
{
  return(sum(log(eigen(R)$values)))	
}

#the inverse of the logit function
inverseLogit <- function(x)
{
  return(exp(x)/(1+exp(x))); 
}

#function for the computation of the Hessian
inverseLogit2 <- function(x)
{
  return(exp(x)/(1+exp(x))^2); 
}

#computes pi_i = P(y_i = 1 | x_i)
getPi <- function(x,beta)
{
  x0 = cbind(rep(1,length(x)),x);
  return(inverseLogit(x0%*%beta));
}

#another function for the computation of the Hessian
getPi2 <- function(x,beta)
{
  x0 = cbind(rep(1,length(x)),x);
  return(inverseLogit2(x0%*%beta));
}

#logistic log-likelihood
logisticLoglik <- function(y,x,beta)
{
  Pi = getPi(x,beta);
  return(sum(y*log(Pi))+sum((1-y)*log(1-Pi)));
}

#calculates l^*(\beta_0,\beta_1)
lStar <- function(y,x,beta)
{
  return(-mean(beta^2)+logisticLoglik(y,x,beta));
}

#obtain the gradient for Newton-Raphson
getGradient <- function(y,x,beta)
{
  gradient = matrix(0,2,1);
  Pi = getPi(x,beta);
  
  gradient[1,1] = sum(y-Pi)-beta[1];
  gradient[2,1] = sum((y-Pi)*x)-beta[2];
  
  return(gradient);
}

#obtain the Hessian for Newton-Raphson
getHessian <- function(y,x,beta)
{
  hessian = matrix(0,2,2);
  Pi2 = getPi2(x,beta);
  
  hessian[1,1] = sum(Pi2)+1;
  hessian[1,2] = sum(Pi2*x);
  hessian[2,1] = hessian[1,2];
  hessian[2,2] = sum(Pi2*x^2)+1;
  
  return(-hessian);
}

#this function implements our own Newton-Raphson procedure
getcoefNR <- function(response,explanatory,data)
{
  epsilon = 1e-10;
  
  #2x1 matrix of coefficients`
  beta = matrix(0,2,1);
  y = data[,response];
  x = data[,explanatory];
  
  #current value of lstar
  currentLStar = lStar(y,x,beta);
  
  #infinite loop unless we stop it someplace inside
  iteration = 0;
  while(1)
  {
	iteration = iteration + 1;
    newBeta = beta - solve(getHessian(y,x,beta))%*%getGradient(y,x,beta);
    newLStar = lStar(y,x,newBeta);
    
    #at each iteration the log-likelihood must increase
    if(newLStar<currentLStar)
    {
		cat('CODING ERROR!! :: ',iteration,' :: ',newLStar,' ',currentLStar,'\n');
      break;
    }
    #stop if the log-likelihood does not improve by too much
    #if(newLStar-currentLStar<1e-6)
    if(max(abs(newBeta-beta))<epsilon)
    {
      break; 
    }
    beta = newBeta;
    currentLStar = newLStar;
  }
  
  return(as.vector(beta));
}

#performs one iteration of the Metropolis-Hastings algorithm
mhLogisticRegression <- function(y,x,beta,invNegHessian)
{
  betaCandidate = mvrnorm(n=1,mu=as.vector(beta),Sigma=invNegHessian);
  
  #current value of lstar
  currentLStar = lStar(y,x,beta);
  
  #values of lstar associated with betaCandidate
  candidateLStar = lStar(y,x,betaCandidate);
  
  if(candidateLStar>=currentLStar)
  {
    return(betaCandidate); 
  }
  
  u = runif(n=1,min=0,max=1);
  if(u<=exp(candidateLStar-currentLStar))
  {
    #accept the move
    return(betaCandidate);
  }
  #reject the move and stay at the current state
  return(beta);
}

#Problem 1: approximate the marginal likelihood (3) using the Laplace approximation
#'betaMode' is the posterior mode obtained using Newton-Raphson
getLaplaceApprox <- function(response,explanatory,data,betaMode)
{
  y = data[,response];
  x = data[,explanatory];
  
  #current value of log-likelihood
  maxLogLik = logisticLoglik(y,x,betaMode);
  
  logmarglik = -mean(betaMode^2)+maxLogLik-0.5*logdet(-getHessian(y,x,betaMode));
  
  return(logmarglik);
}

#Numerically calculate the marginal likelihood (3) using Monte Carlo (numerically unstable version)
getMonteCarloRaw <- function(response,explanatory,data,NumberOfIterations)
{
  y = data[,response];
  x = data[,explanatory];

  logmarglik = 0;
  for(i in 1:NumberOfIterations)
  {
    logmarglik = logmarglik+exp(logisticLoglik(y,x,rnorm(n=2)));
  }
  
  return(log(logmarglik)-log(NumberOfIterations));
}

#Numerically calculate the marginal likelihood (3) using Monte Carlo (numerically stable version)
getMonteCarlo <- function(response,explanatory,data,NumberOfIterations)
{
  y = data[,response];
  x = data[,explanatory];

  loglikVec = vector(mode="numeric",length=NumberOfIterations);
  for(i in 1:NumberOfIterations)
  {
    loglikVec[i] = logisticLoglik(y,x,rnorm(n=2));
  }
  maxloglikVec = max(loglikVec);
  loglikVec = exp(loglikVec-maxloglikVec);
  logmarglik = log(mean(loglikVec))+maxloglikVec;
  
  return(logmarglik);
}

#Problem 2 calculate the posterior means of beta0 and beta1 by
#sampling from the joint posterior distribution of beta0 and beta1
#with the Metropolis-Hastings algorithm
getPosteriorMeans <- function(response,explanatory,data,betaMode,NumberOfIterations)
{
   y = data[,response];
   x = data[,explanatory];
 
   betaBayes = matrix(0,nrow=2,ncol=1);
   betaCurrent = betaMode;
 
   #the negative of the inverse Hessian matrix evaluated at the posterior mode
   invNegHessian = -solve(getHessian(y,x,betaMode));
   
   for(iteration in 1:NumberOfIterations)
   {
      betaCurrent =  mhLogisticRegression(y,x,betaCurrent,invNegHessian);
      betaBayes = betaBayes+betaCurrent;
   }
   return(betaBayes/NumberOfIterations);
}