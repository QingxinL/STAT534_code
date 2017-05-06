mleLogistic <- function(explanatory,response,data)
{
  #check if the regression has no explanatory variables
  if(0==length(explanatory))
  {
    #regression with no explanatory variables
    return(as.numeric(coef(glm(data[,response] ~ 1,family=binomial(link=logit)))));
  }  
  return(as.numeric(coef(glm(data[,response] ~ data[,explanatory],family=binomial(link=logit)))));
}



# calculate the first derivative
logitDeri <- function()
{
  
}

# calculate the Hessian matriax of l*
hessian <- function(explanatory, response, data, beta0, beta1)
{
  ytemp = beta0 + data[, explanatory]*beta1;
  logitFirstDeri = exp(ytemp)/((1+exp(ytemp))^2);
  a = 1 - sum(logitFirstDeri);
  b = - sum(logitFirstDeri*data[, explanatory]);
  c = b;
  d = - 1 - sum(logitFirstDeri*(data[, explanatory]^2));
  
  hessianMat = matrix(c(a, c ,b ,d), ncol=2);
  return(hessianMat);
}

# calculate the delta(l*)
deriOfLikeliHood <- function(explanatory, response, data, beta0, beta1)
{
  pi_i = inv.logit(beta0 + beta1*data[, explanatory]);
  deri0 = sum(data[, response] - pi_i);
  deri1 = sum(data[, response]*data[, explanatory] - pi_i*data[, explanatory]);
  deri = c(deri0, deri1);
  return(deri);
}

# use the Newton-Raphson algorithm to perform the estimation
# iterLimit  = 0.0001
# use the Newton-Raphson algorithm to perform the estimation
# iterLimit  = 0.0001
newtonRaphson <- function(explanatory, response, data, iterLimit)
{
  stopFlag = FALSE;
  beta_last = c(0, 0);
  beta_new = c(0, 0);
  k = 0;
  while(!stopFlag)
  {
    k = k + 1;
    hessianMat = hessian(explanatory, response, data, beta_last[1], beta_last[2]);
    deriLH = deriOfLikeliHood(explanatory, response, data, beta_last[1], beta_last[2]);
    beta_new = beta_last - solve(hessianMat) %*% deriLH;
    error0 = beta_new[1] - beta_last[1];
    error1 = beta_new[2] - beta_last[2];
    beta_last = beta_new;
    if ((error0 < iterLimit) && (error1 < iterLimit)) 
      stopFlag = TRUE;
    # print(k);
    # print(error0);
    # print(error1);
  }
  return(beta_new);
}

# this function calculate the likelyhood of l
likeliHood <- function(explanatory, response, data, beta0, beta1)
{
  num = length(explanatory);
  sumLikeHood = 0;
  # for (i in c(1:num))
  # {
  #   pi_i = inv.logit(beta0 + beta1*data[i, explanatory]);
  #   sumLikeHood = sumLikeHood + (data[i, response] * log(pi_i) + (1 - data[i, response])*log(1 - pi_i));
  # }
  pi_i = inv.logit(beta0 + beta1*data[, explanatory]);
  sumLikeHood = sum(data[, response]*log(pi_i) + (1 - data[, response])*log(1 - pi_i));
  return(sumLikeHood);
}
# l = likeliHood(1, 61, data, 0, 0)

# this function calculate the estimated likelyhood l*
estLikeliHood  <- function(explanatory, response, data, beta0, beta1)
{
  estLH = 0;
  estLH = - log(2*pi) - 0.5*(beta0^2 + beta1^2) + likeliHood(explanatory, response, data, beta0, beta1);
  return(estLH);
}

# this function calculate the logarithm of the marginal likelihood log(p(D))
laplaceLogLik <- function(explanatory, response, data)
{
  # num = length(explanatory);
  marginLH = 0;
  iterLimit = 0.0001;
  beta = newtonRaphson(explanatory, response, data, iterLimit);
  estLH = estLikeliHood(explanatory, response, data, beta[1], beta[2]);
  hessianMat = hessian(explanatory, response, data, beta[1], beta[2]);
  # print(beta);
  # print(estLH);
  # print(hessianMat);
  marginLH = log(2*pi) + (estLH) - 0.5*log(det(-hessianMat));
  
  return(marginLH);
}

# sample u from uniform (0,1) distribution and make the decision
sampleU <- function(lhCandidate, lhCurrent)
{
  u = runif(1, min = 0, max = 1);
  if (log(u) <= (lhCandidate - lhCurrent))
  {
      return(TRUE)  
  }
  else
  {
      return(FALSE)
  }
}

# implement the Metropolis-Hastings algorithm
metropolisHastings <- function(explanatory, response, data, iterations)
{
  iterLimit = 0.0001;
  betaInit = newtonRaphson(explanatory, response, data, iterLimit);
  betaCurrent = betaInit;
  betaNext = c(0, 0);
  betaCandidate = c(0, 0);
  sumBeta = c(0, 0);
  
  # uptate the current state of the Markov chain 
  for (i in c(1:iterations))
  {
    hessianMat = hessian(explanatory, response, data, betaInit[1], betaInit[2]);
    sigma = -solve(hessianMat);
    betaCandidate = mvrnorm(n=1, betaInit, sigma);
    lhCandidate = estLikeliHood(explanatory, response, data, betaCandidate[1], betaCandidate[2]);
    lhCurrent = estLikeliHood(explanatory, response, data, betaCurrent[1], betaCurrent[2]);
    
    if (lhCandidate >= lhCurrent)
    {
      betaNext = betaCandidate;
    }
    else if (sampleU(lhCandidate, lhCurrent))
    {
      betaNext = betaCandidate;
    }
    else if(!sampleU(lhCandidate, lhCurrent))
    {
      betaNext = betaCurrent;
    }
    betaCurrent = betaNext;
    sumBeta = sumBeta + betaCurrent;
  }
  
  meanBeta = sumBeta/iterations;
  return(meanBeta);
  
}
  
bayesLogistic <- function(apredictor,response,data,NumberOfIterations)
{
  #source('bayeslogisticreg-startcode-HW4.R');
  print(appredictor);l
  explanatory = apredictor;
  beta_mean = metropolisHastings(explanatory, response, data, NumberOfIterations);
  beta_mle = mleLogistic(explanatory,response,data);
  logmarglik = laplaceLogLik(explanatory, response, data);
  result = list('apredictor' = apredictor, 'logmarglik' = logmarglik, 'beta0bayes' = beta_mean[1], 'beta1bayes' = beta_mean[2],
              'beta0mle' = beta_mle[1], 'beta1mle' = beta_mle[2]);
  return(result);
}

#PARALLEL VERSION
#datafile = the name of the file with the data
#NumberOfIterations = number of iterations of the Metropolis-Hastings algorithm
#clusterSize = number of separate processes; each process performs one or more
#univariate regressions
main <- function(datafile,NumberOfIterations,clusterSize)
{
  #read the data
  data = read.table(datafile,header=FALSE);
  
  #the sample size is 148 (number of rows)
  #the explanatory variables are the first 60 columns for '534binarydata.txt'
  #the last column is the binary response
  response = ncol(data);
  lastPredictor = ncol(data)-1;
  
  #initialize a cluster for parallel computing
  cluster <- makeCluster(clusterSize, type = "SOCK")
  
  library(MASS)
  # import function to the cluster
  clusterExport(cluster, list("mleLogistic", "hessian", "deriOfLikeliHood", "newtonRaphson", "likeliHood", 
                              "estLikeliHood", "laplaceLogLik", "sampleU", "metropolisHastings", "bayesLogistic",
                              "mvrnorm" ))
  #run the MC3 algorithm from several times
  results = clusterApply(cluster, 1:lastPredictor, bayesLogistic,
                         response,data,NumberOfIterations);
  
  #print out the results
  for(i in 1:lastPredictor)
  {
    cat('Regression of Y on explanatory variable ',results[[i]]$apredictor,
        ' has log marginal likelihood ',results[[i]]$logmarglik,
        ' with beta0 = ',results[[i]]$beta0bayes,' (',results[[i]]$beta0mle,')',
        ' and beta1 = ',results[[i]]$beta1bayes,' (',results[[i]]$beta1mle,')',
        '\n');    
  }
  
  #destroy the cluster
  stopCluster(cluster);  
}
#install.packages("gtools");
require(gtools);
require(MASS);
setwd("~/Course/STAT534/STAT534_code/HW4");

#NOTE: YOU NEED THE PACKAGE 'SNOW' FOR PARALLEL COMPUTING
#install.packages("snow");
require(snow);

#this is where the program starts
main('534binarydata.txt',10000,2)

