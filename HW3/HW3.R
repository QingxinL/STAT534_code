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
    AIC_Ap = getLogisticAIC(response, A_p, data);
    AIC_Ar = getLogisticAIC(response, A_r, data);
    AIC_Best = getLogisticAIC(response, bestModel, data);
    p_Ap = - AIC_Ap - log(length(nbd_Ap));
    p_Ar = - AIC_Ar - log(length(nbd_Ar));

    # step 7
    if (p_Ap > p_Ar)
    {
      A_r = A_p;
      if (AIC_Ap < AIC_Best)
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
        if (AIC_Ap < AIC_Best)
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

  explanatory = c(1:lastPredictor);
  result = modelSelectionMC3(response, explanatory, data, 25);

  cat('\n\nbest model is ', result$model);
  cat('\nBest AIC = ', result$AIC);
  # Problem 1
  # for the small dataset, the result is

  for (i in c(1:19))
  {
    result = modelSelectionMC3(response, explanatory, data, 25);

    cat('\n\nbest model is ', result$model);
    cat('\nBest AIC = ', result$AIC);
  }

}
tempdir()

setwd("~/Course/STAT534/STAT534_code/HW3")
main('534binarydata.txt')

##############################################################
####### PROBLEM 1
##############################################################
# for small dataset
best model is  10 2 9 7 8 5 1 4 6 3
Best AIC =  112.8834

best model is  5 2 7 6 10 3 4 9 1 8
Best AIC =  112.8834
The results are very similiar to greedy search.
#
# Problem 2
# The effiency of this algorithm is not very good, so I got results for 5 instances.
# The results of MCCC algorithm are different due to its random itialize parameters and sample.
# For this dataset, the algorithm performs better than greedy search
# best model is  12 55 27 26 58 35 39 28 30 29 47 44 41 4 50 49 3 56 2 52 53 10 40 20 15 34
# Best AIC =  83.97523
#
# best model is  15 19 1 44 20 25 3 59 17 21 28 36 30 45
# Best AIC =  90.71816
#
# best model is  12 8 52 48 53 15 19 45 4 36 57 16 22 31 26 24 38 43 14 21 55 25 30 59 39 9 10 58
# Best AIC =  131.0038
#
# best model is  18 54 60 6 47 53 25 49 14 59 17 8 1 42 50 38 33 21 55 9 22 39 30 4 24
# Best AIC =  111.8044
#
# best model is  50 20 8 47 16 37 29 33 13 25 21 17 49 24 54 36 19 2
# Best AIC =  100.8705
