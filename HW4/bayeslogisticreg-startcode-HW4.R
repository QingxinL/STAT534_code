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
                              "mvrnorm", "inv.logit" ));
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
main('534binarydata.txt',10000,10)


#result 3
# Regression of Y on explanatory variable  1  has log marginal likelihood  -84.13831  with beta0 =  -0.770417  ( -0.7845448 )  and beta1 =  1.020589  ( 1.042284 ) 
# Regression of Y on explanatory variable  2  has log marginal likelihood  -91.47557  with beta0 =  -0.7517921  ( -0.8117805 )  and beta1 =  -0.57384  ( -0.683414 ) 
# Regression of Y on explanatory variable  3  has log marginal likelihood  -89.23678  with beta0 =  -0.8018777  ( -0.8174069 )  and beta1 =  0.7551913  ( 0.7676812 ) 
# Regression of Y on explanatory variable  4  has log marginal likelihood  -94.71905  with beta0 =  -0.7491684  ( -0.7581619 )  and beta1 =  0.3773244  ( 0.3837138 ) 
# Regression of Y on explanatory variable  5  has log marginal likelihood  -92.75307  with beta0 =  -0.7723428  ( -0.7836879 )  and beta1 =  0.5404678  ( 0.5479422 ) 
# Regression of Y on explanatory variable  6  has log marginal likelihood  -94.41787  with beta0 =  -0.7489349  ( -0.7575917 )  and beta1 =  0.4003865  ( 0.4039104 ) 
# Regression of Y on explanatory variable  7  has log marginal likelihood  -95.37601  with beta0 =  -0.7439262  ( -0.7514321 )  and beta1 =  0.3187571  ( 0.3231848 ) 
# Regression of Y on explanatory variable  8  has log marginal likelihood  -87.98998  with beta0 =  -0.7979223  ( -0.8103876 )  and beta1 =  0.7988585  ( 0.8091978 ) 
# Regression of Y on explanatory variable  9  has log marginal likelihood  -95.41706  with beta0 =  -0.7232373  ( -0.7355123 )  and beta1 =  0.3254644  ( 0.3120068 ) 
# Regression of Y on explanatory variable  10  has log marginal likelihood  -94.80952  with beta0 =  -0.7482577  ( -0.7567287 )  and beta1 =  0.3800449  ( 0.3830711 ) 
# Regression of Y on explanatory variable  11  has log marginal likelihood  -91.95151  with beta0 =  -0.7798327  ( -0.7920907 )  and beta1 =  0.581219  ( 0.5904112 ) 
# Regression of Y on explanatory variable  12  has log marginal likelihood  -96.1376  with beta0 =  -0.7251639  ( -0.7441575 )  and beta1 =  -0.2143411  ( -0.2387211 ) 
# Regression of Y on explanatory variable  13  has log marginal likelihood  -92.10862  with beta0 =  -0.7710163  ( -0.7841235 )  and beta1 =  0.5656302  ( 0.5733513 ) 
# Regression of Y on explanatory variable  14  has log marginal likelihood  -95.36066  with beta0 =  -0.7423282  ( -0.751308 )  and beta1 =  0.3204262  ( 0.3237944 ) 
# Regression of Y on explanatory variable  15  has log marginal likelihood  -91.2233  with beta0 =  -0.8063615  ( -0.8197562 )  and beta1 =  0.7008603  ( 0.7177387 ) 
# Regression of Y on explanatory variable  16  has log marginal likelihood  -93.45905  with beta0 =  -0.7675068  ( -0.7760113 )  and beta1 =  0.4948745  ( 0.5017284 ) 
# Regression of Y on explanatory variable  17  has log marginal likelihood  -91.5717  with beta0 =  -0.7497155  ( -0.8022312 )  and beta1 =  -0.5568936  ( -0.6565506 ) 
# Regression of Y on explanatory variable  18  has log marginal likelihood  -89.05816  with beta0 =  -0.8007475  ( -0.8123199 )  and beta1 =  0.735877  ( 0.7476675 ) 
# Regression of Y on explanatory variable  19  has log marginal likelihood  -90.39235  with beta0 =  -0.7632465  ( -0.775153 )  and beta1 =  0.7524889  ( 0.7655354 ) 
# Regression of Y on explanatory variable  20  has log marginal likelihood  -94.83722  with beta0 =  -0.7451763  ( -0.7564703 )  and beta1 =  0.3693161  ( 0.3733819 ) 
# Regression of Y on explanatory variable  21  has log marginal likelihood  -85.71701  with beta0 =  -0.8377385  ( -0.8553298 )  and beta1 =  0.9088393  ( 0.9283823 ) 
# Regression of Y on explanatory variable  22  has log marginal likelihood  -84.00894  with beta0 =  -0.8510534  ( -0.8712209 )  and beta1 =  0.9858494  ( 1.009002 ) 
# Regression of Y on explanatory variable  23  has log marginal likelihood  -79.4728  with beta0 =  -0.9208662  ( -0.9460136 )  and beta1 =  1.266981  ( 1.307546 ) 
# Regression of Y on explanatory variable  24  has log marginal likelihood  -94.2046  with beta0 =  -0.7335245  ( -0.764407 )  and beta1 =  -0.3923911  ( -0.4373921 ) 
# Regression of Y on explanatory variable  25  has log marginal likelihood  -96.5519  with beta0 =  -0.7294977  ( -0.7387917 )  and beta1 =  0.1662188  ( 0.169028 ) 
# Regression of Y on explanatory variable  26  has log marginal likelihood  -91.70005  with beta0 =  -0.7554251  ( -0.7676342 )  and beta1 =  0.5818631  ( 0.5907662 ) 
# Regression of Y on explanatory variable  27  has log marginal likelihood  -90.64171  with beta0 =  -0.7549472  ( -0.8105548 )  and beta1 =  -0.6028632  ( -0.6907203 ) 
# Regression of Y on explanatory variable  28  has log marginal likelihood  -95.18164  with beta0 =  -0.7450695  ( -0.7536521 )  and beta1 =  0.3388201  ( 0.3426751 ) 
# Regression of Y on explanatory variable  29  has log marginal likelihood  -91.9578  with beta0 =  -0.7797926  ( -0.7900873 )  and beta1 =  0.5802032  ( 0.5901491 ) 
# Regression of Y on explanatory variable  30  has log marginal likelihood  -90.77594  with beta0 =  -0.803072  ( -0.8176779 )  and beta1 =  0.6935169  ( 0.7076027 ) 
# Regression of Y on explanatory variable  31  has log marginal likelihood  -90.85716  with beta0 =  -0.7900831  ( -0.805714 )  and beta1 =  0.6521466  ( 0.6639053 ) 
# Regression of Y on explanatory variable  32  has log marginal likelihood  -91.10308  with beta0 =  -0.7613303  ( -0.8358587 )  and beta1 =  -0.6119996  ( -0.7765949 ) 
# Regression of Y on explanatory variable  33  has log marginal likelihood  -95.70375  with beta0 =  -0.7427332  ( -0.748413 )  and beta1 =  0.2866353  ( 0.2897193 ) 
# Regression of Y on explanatory variable  34  has log marginal likelihood  -89.58609  with beta0 =  -0.7662572  ( -0.8418246 )  and beta1 =  -0.6676208  ( -0.804528 ) 
# Regression of Y on explanatory variable  35  has log marginal likelihood  -96.65826  with beta0 =  -0.7243931  ( -0.7376418 )  and beta1 =  0.1449705  ( 0.1481233 ) 
# Regression of Y on explanatory variable  36  has log marginal likelihood  -96.70551  with beta0 =  -0.7283896  ( -0.7371812 )  and beta1 =  0.1354769  ( 0.1368646 ) 
# Regression of Y on explanatory variable  37  has log marginal likelihood  -83.39354  with beta0 =  -0.8616591  ( -0.8820737 )  and beta1 =  1.031753  ( 1.055489 ) 
# Regression of Y on explanatory variable  38  has log marginal likelihood  -93.60809  with beta0 =  -0.7420307  ( -0.779148 )  and beta1 =  -0.4378426  ( -0.5022609 ) 
# Regression of Y on explanatory variable  39  has log marginal likelihood  -87.46206  with beta0 =  -0.8275944  ( -0.842756 )  and beta1 =  0.8534551  ( 0.8762895 ) 
# Regression of Y on explanatory variable  40  has log marginal likelihood  -90.97866  with beta0 =  -0.7895535  ( -0.7994708 )  and beta1 =  0.6521353  ( 0.6628068 ) 
# Regression of Y on explanatory variable  41  has log marginal likelihood  -91.38773  with beta0 =  -0.7443477  ( -0.7917664 )  and beta1 =  -0.5613535  ( -0.6275432 ) 
# Regression of Y on explanatory variable  42  has log marginal likelihood  -86.50954  with beta0 =  -0.7769127  ( -0.8700183 )  and beta1 =  -0.795791  ( -0.9696883 ) 
# Regression of Y on explanatory variable  43  has log marginal likelihood  -94.20598  with beta0 =  -0.7357722  ( -0.7649961 )  and beta1 =  -0.3901267  ( -0.4361109 ) 
# Regression of Y on explanatory variable  44  has log marginal likelihood  -95.88853  with beta0 =  -0.7370846  ( -0.7460102 )  and beta1 =  0.2623494  ( 0.2667709 ) 
# Regression of Y on explanatory variable  45  has log marginal likelihood  -95.1518  with beta0 =  -0.7266736  ( -0.7562733 )  and beta1 =  -0.3175296  ( -0.3564728 ) 
# Regression of Y on explanatory variable  46  has log marginal likelihood  -86.5804  with beta0 =  -0.8423751  ( -0.8607597 )  and beta1 =  0.9349891  ( 0.9587251 ) 
# Regression of Y on explanatory variable  47  has log marginal likelihood  -91.59112  with beta0 =  -0.7770919  ( -0.7922888 )  and beta1 =  0.599  ( 0.6084323 ) 
# Regression of Y on explanatory variable  48  has log marginal likelihood  -91.4343  with beta0 =  -0.7821331  ( -0.7922874 )  and beta1 =  0.6120124  ( 0.6209276 ) 
# Regression of Y on explanatory variable  49  has log marginal likelihood  -90.38113  with beta0 =  -0.7818578  ( -0.7941998 )  and beta1 =  0.6556444  ( 0.6646216 ) 
# Regression of Y on explanatory variable  50  has log marginal likelihood  -92.52367  with beta0 =  -0.7385583  ( -0.7799963 )  and beta1 =  -0.4986917  ( -0.559328 ) 
# Regression of Y on explanatory variable  51  has log marginal likelihood  -95.32861  with beta0 =  -0.727822  ( -0.7515387 )  and beta1 =  -0.3011734  ( -0.3315456 ) 
# Regression of Y on explanatory variable  52  has log marginal likelihood  -93.76152  with beta0 =  -0.7601196  ( -0.7692825 )  and beta1 =  0.4620891  ( 0.4670527 ) 
# Regression of Y on explanatory variable  53  has log marginal likelihood  -94.00387  with beta0 =  -0.7403301  ( -0.7741983 )  and beta1 =  -0.4091363  ( -0.4711203 ) 
# Regression of Y on explanatory variable  54  has log marginal likelihood  -91.22528  with beta0 =  -0.7921889  ( -0.8066512 )  and beta1 =  0.6539473  ( 0.6691851 ) 
# Regression of Y on explanatory variable  55  has log marginal likelihood  -90.65106  with beta0 =  -0.7919779  ( -0.8059926 )  and beta1 =  0.6638834  ( 0.6758964 ) 
# Regression of Y on explanatory variable  56  has log marginal likelihood  -90.92952  with beta0 =  -0.7727167  ( -0.7845683 )  and beta1 =  0.6243401  ( 0.6304744 ) 
# Regression of Y on explanatory variable  57  has log marginal likelihood  -95.2672  with beta0 =  -0.7436732  ( -0.7522831 )  and beta1 =  0.3276178  ( 0.3336013 ) 
# Regression of Y on explanatory variable  58  has log marginal likelihood  -92.41541  with beta0 =  -0.761188  ( -0.774532 )  and beta1 =  0.5434923  ( 0.5533604 ) 
# Regression of Y on explanatory variable  59  has log marginal likelihood  -96.99462  with beta0 =  -0.7227983  ( -0.7339938 )  and beta1 =  0.01172908  ( 0.01186692 ) 
# Regression of Y on explanatory variable  60  has log marginal likelihood  -94.03109  with beta0 =  -0.73104  ( -0.7647065 )  and beta1 =  -0.4025531  ( -0.4500869 ) 
