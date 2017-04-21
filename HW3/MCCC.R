getLogisticAIC <- function(response,explanatory,data)
{
  #check if the regression has no explanatory variables
  if(0==length(explanatory))
  {
    #regression with no explanatory variables
    deviance = glm(data[,response] ~ 1,family=binomial(link=logit))$deviance;
  }
  else {
    #regression with at least one explanatory variable
    deviance = glm(data[,response] ~ as.matrix(data[,as.numeric(explanatory)]),
                   
                   family=binomial(link=logit))$deviance;
  }  
  return(deviance+2*(1+length(explanatory)));
  
}

#this function checks whether all fitted values are strictly greater
#than 0 and strictly smaller than 1
#this is an empirical way to check if the MLEs exist
isValidLogistic <- function(response,explanatory,data)
{
  epsilon = 1e-20;
  
  if(0==length(explanatory))
  {
    #regression with no explanatory variables
    fittedValues = fitted(glm(data[,response] ~ 1,family=binomial(link=logit)));
  }
  else
  {
    #regression with at least one explanatory variable
    fittedValues = fitted(glm(data[,response] ~ as.matrix(data[,as.numeric(explanatory)]),family=binomial(link=logit)));
  }
  
  if(all(fittedValues>epsilon) && all(fittedValues<1-epsilon))
  {
    return(TRUE); #MLES are well determined 
  }
  return(FALSE); #MLES are not well determined
}

#this is the version of the 'isValidLogistic' function
#based on Charles Geyers RCDD package
isValidLogisticRCDD <- function(response,explanatory,data)
{  
  logisticreg = suppressWarnings(glm(data[,response] ~ as.matrix(data[,as.numeric(explanatory)]),family=binomial(link=logit),x=TRUE));
  tanv = logisticreg$x;
  tanv[data[,response] == 1, ] <- (-tanv[data[,response] == 1, ]);
  vrep = cbind(0, 0, tanv);
  #with exact arithmetic; takes a long time
  #lout = linearity(d2q(vrep), rep = "V");
  
  lout = linearity(vrep, rep = "V");
  return(length(lout)==nrow(data));
}

lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

# identifay the valid neighbors of the current logistic regression
nbdValid <- function(response, A_r, explanatory, data)
{
  nbdA_r = list();
  nbdInd = 0;
  # add one
  addNbdSet = setdiff(explanatory, A_r);
  for (addNbd in addNbdSet)
  {
    newSet = c(A_r, addNbd);
    if (isValidLogisticRCDD(response, newSet, data)) 
    {
      nbdInd = nbdInd +1 ;
      nbdA_r[[nbdInd]] = newSet;
    };
  }
  # delete one
  for (deleteNbd in A_r)
  {
    newSet = setdiff(A_r, deleteNbd);
    if (isValidLogisticRCDD(response, newSet, data)) 
    {
      nbdInd = nbdInd +1 ;
      nbdA_r[[nbdInd]] = newSet;
    };
  }
  #return
  return(nbdA_r)
}

modelSelectionMC3 <- function(response, explanatory, data, iterations)
{
  # first iteration
  k = sample.int(length(explanatory), 1);
  regValid = FALSE;
  while(!regValid)
  {
    A_0 = sample(explanatory, k);
    regValid = isValidLogisticRCDD(response, A_0, data);
  }
  bestModel = A_0;
  
  # iteration
  
  A_r = A_0;
  iter = 0;
  while(iter < iterations)
  {
    iter = iter + 1;
    cat('\niter = ', iter);
    # step 1,2
    nbd_Ar = nbdValid(response, A_r, explanatory, data);
    # step 3
    A_p = unlist(sample(nbd_Ar))
    # step 4
    nbd_Ap = nbdValid(response, A_p, explanatory, data);
    
    # step 5, 6
    p_Ap = -getLogisticAIC(response, A_p, data) - log(length(nbd_Ap));
    p_Ar = -getLogisticAIC(response, A_r, data) - log(length(nbd_Ar)); 
    
    # step 7
    if (p_Ap > p_Ar) 
    {
      A_r = A_p;
      if (getLogisticAIC(response, A_p, data) < getLogisticAIC(response, bestModel, data))
      {
        bestModel = A_p;
      }
      else
      {
        bestModel = bestModel;
      }
    }
    else #step 8
    {
      u = runif(1);
      if (log(u) < (p_Ap - p_Ar)) 
      {
        A_r = A_p;
        if (getLogisticAIC(response, A_p, data) < getLogisticAIC(response, bestModel, data))
        {
          bestModel = A_p;
        }
        else
        {
          bestModel = bestModel;
        }
      }
      else
      {
        A_r = A_r;
        bestModel;
      }
    }
    
  }
  AICfinal = getLogisticAIC(response, bestModel, data);
  return(list(model = bestModel, AIC = AICfinal));
}

# Main ----------------------------------------------------------------


main <- function(datafile)
{
  #read the data
  data = read.table(datafile, header=FALSE);
  
  #the sample size is 148 (number of rows)
  #the explanatory variables are the first 60 columns
  #the last column is the binary response
  response = ncol(data);
  lastPredictor = ncol(data) - 1;
  
  
  result = modelSelectionMC3(response, explanatory, data, 25);
  
  cat('\n best model is ', result$model);
  cat('\n\nBest AIC = ', result$AIC);
  # Problem 1
  # for the small dataset, the result is 
  # 
  # Problem 2
  for (i in c(1:9))
  {
    result = modelSelectionMC3(response, explanatory, data, 25);
    
    cat('\n best model is ', result$model);
    cat('\n\nBest AIC = ', result$AIC);
  }
  
}
# explanatory  = c(1:3);
# isValidLogisticRCDD(response, c(1:3), data)

setwd("~/Course/STAT534/STAT534_code/HW3")
main('534binarydatasmall.txt')