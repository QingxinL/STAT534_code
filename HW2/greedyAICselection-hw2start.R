#HOMEWORK 2, PROBLEM 1
#this function uses ’glm’ to fit a logistic regression
#and returns the BIC = deviance + log(SampleSize)*NumberOfRegressionCoefficients
getLogisticBIC <- function(response,explanatory,data)
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
  return(deviance+log(nrow(data))*(1+length(explanatory)));
}


#this function uses 'glm' to fit a logistic regression
#and returns the AIC = deviance + 2*NumberOfCoefficients 
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
# AIC ----------------------------------------------------------------

#HOMEWORK 2, PROBLEM 2: Forward greedy search
forwardSearchAIC <- function(response,data,lastPredictor)
{
  
  #start with the empty regression with no predictors
  bestRegression = NULL;
  #calculate the AIC of the empty regression
  bestRegressionAIC = getLogisticAIC(response,bestRegression,data);
  cat('\n\n\n\nforwardSearch :: The empty logistic regression has AIC = ',bestRegressionAIC,'\n');
  
  #vector that keeps track of all the variables
  #that are currently NOT in the model
  VariablesNotInModel = 1:lastPredictor;
  
  #add variables one at a time to the current regression
  #and retain that variable that gives the smallest values of AIC associated
  #Make the model that includes that variable to be the current model
  #if it gives a smaller AIC than the AIC of the current regression
  
  #stop when there are no variables that can be included in the model
  stepNumber = 0;
  while(length(VariablesNotInModel)>=1)
  {
    #record the number of steps performed
    stepNumber = stepNumber + 1;
    
    #create a vector that records the AIC values of the regressions
    #we are examining; the number of these regressions is equal
    #with the number of variables that are not in the model
    regAIC = vector('numeric',length(VariablesNotInModel));
    
    #take each variable that is not in the model
    #and include it in the model
    ind = 0;
    for (i in VariablesNotInModel){
      # add a variable to current regression model
      ind = ind + 1;
      addRegression = c(bestRegression, i);
      regAIC[ind] =  getLogisticAIC(response,addRegression,data);
    }
    #find mininum regAIC, add to current regression model
    addToModel = VariablesNotInModel[which.min(regAIC)];
    newRegression = c(bestRegression, addToModel);
    
    currentAIC = getLogisticAIC(response,bestRegression,data);
    newAIC = getLogisticAIC(response,newRegression,data);
    
    #if the new variable can reduce the AIC 
    if (newAIC < currentAIC){
      bestRegression = newRegression;
      bestRegressionAIC = newAIC;
      
      #remove the variable in VariablesNotInModel
      VariablesNotInModel = VariablesNotInModel[!VariablesNotInModel %in% addToModel];
    }
    else break;
  }
  
  return(list(aic=bestRegressionAIC,reg=bestRegression));
}

##HOMEWORK 2, PROBLEM 3: Backward greedy search
backwardSearchAIC <- function(response,data,lastPredictor)
{
  #start with the full regression that includes all the variables
  bestRegression = 1:lastPredictor;
  #calculate the AIC of the full regression
  bestRegressionAIC = getLogisticAIC(response,bestRegression,data);
  cat('\n\n\n\nbackwardSearch :: The full logistic regression has AIC = ',bestRegressionAIC,'\n');
  
  #sequentially delete one variable from the current regression
  #and retain that variable that gives the smallest AIC; make the model
  #in which that variable is deleted to be the current model if
  #this leads to a current model with a smaller AIC
  stepNumber = 0;
  while(length(bestRegression)>=1)
  {
    #record the number of steps performed
    stepNumber = stepNumber + 1;
    
    #create a vector that records the AIC values of the regressions
    #we are examining; the number of these regressions is equal
    #with the number of variables that are in the model
    regAIC = vector('numeric',length(bestRegression));
    
    ind = 0;
    for (i in bestRegression){
      # remove a variable to current regression model
      ind = ind + 1;
      deleteRegression = bestRegression[!bestRegression %in% i];
      regAIC[ind] =  getLogisticAIC(response,deleteRegression,data);
    }
    #find mininum regAIC, remove from current regression model
    delToModel = bestRegression[which.min(regAIC)];
    newRegression = bestRegression[!bestRegression %in% delToModel];
    
    currentAIC = getLogisticAIC(response,bestRegression,data);
    newAIC = getLogisticAIC(response,newRegression,data);
    
    #if the new regression model can reduce the AIC 
    if (newAIC < currentAIC){
      bestRegression = newRegression;
      bestRegressionAIC = newAIC;
    }
    else break;
  }
  
  return(list(aic=bestRegressionAIC,reg=bestRegression));
}

# BIC---------------------------------------------------
###

#HOMEWORK 2, PROBLEM 2: Forward greedy search
forwardSearchBIC <- function(response,data,lastPredictor)
{
  
  #start with the empty regression with no predictors
  bestRegression = NULL;
  #calculate the BIC of the empty regression
  bestRegressionBIC = getLogisticBIC(response,bestRegression,data);
  cat('\n\n\n\nforwardSearch :: The empty logistic regression has BIC = ',bestRegressionBIC,'\n');
  
  #vector that keeps track of all the variables
  #that are currently NOT in the model
  VariablesNotInModel = 1:lastPredictor;
  
  #add variables one at a time to the current regression
  #and retain that variable that gives the smallest values of BIC associated
  #Make the model that includes that variable to be the current model
  #if it gives a smaller BIC than the BIC of the current regression
  
  #stop when there are no variables that can be included in the model
  stepNumber = 0;
  while(length(VariablesNotInModel)>=1)
  {
    #record the number of steps performed
    stepNumber = stepNumber + 1;
    
    #create a vector that records the BIC values of the regressions
    #we are examining; the number of these regressions is equal
    #with the number of variables that are not in the model
    regBIC = vector('numeric',length(VariablesNotInModel));
    
    #take each variable that is not in the model
    #and include it in the model
    ind = 0;
    for (i in VariablesNotInModel){
      # add a variable to current regression model
      ind = ind + 1;
      addRegression = c(bestRegression, i);
      regBIC[ind] =  getLogisticBIC(response,addRegression,data);
    }
    #find mininum regBIC, add to current regression model
    addToModel = VariablesNotInModel[which.min(regBIC)];
    newRegression = c(bestRegression, addToModel);
    
    currentBIC = getLogisticBIC(response,bestRegression,data);
    newBIC = getLogisticBIC(response,newRegression,data);
    
    #if the new variable can reduce the BIC 
    if (newBIC < currentBIC){
      bestRegression = newRegression;
      bestRegressionBIC = newBIC;
      
      #remove the variable in VariablesNotInModel
      VariablesNotInModel = VariablesNotInModel[!VariablesNotInModel %in% addToModel];
    }
    else break;
  }
  
  return(list(aic=bestRegressionBIC,reg=bestRegression));
}

##HOMEWORK 2, PROBLEM 3: Backward greedy search
backwardSearchBIC <- function(response,data,lastPredictor)
{
  #start with the full regression that includes all the variables
  bestRegression = 1:lastPredictor;
  #calculate the BIC of the full regression
  bestRegressionBIC = getLogisticBIC(response,bestRegression,data);
  cat('\n\n\n\nbackwardSearch :: The full logistic regression has BIC = ',bestRegressionBIC,'\n');
  
  #sequentially delete one variable from the current regression
  #and retain that variable that gives the smallest BIC; make the model
  #in which that variable is deleted to be the current model if
  #this leads to a current model with a smaller BIC
  stepNumber = 0;
  while(length(bestRegression)>=1)
  {
    #record the number of steps performed
    stepNumber = stepNumber + 1;
    
    #create a vector that records the BIC values of the regressions
    #we are examining; the number of these regressions is equal
    #with the number of variables that are in the model
    regBIC = vector('numeric',length(bestRegression));
    
    ind = 0;
    for (i in bestRegression){
      # remove a variable to current regression model
      ind = ind + 1;
      deleteRegression = bestRegression[!bestRegression %in% i];
      regBIC[ind] =  getLogisticBIC(response,deleteRegression,data);
    }
    #find mininum regBIC, remove from current regression model
    delToModel = bestRegression[which.min(regBIC)];
    newRegression = bestRegression[!bestRegression %in% delToModel];
    
    currentBIC = getLogisticBIC(response,bestRegression,data);
    newBIC = getLogisticBIC(response,newRegression,data);
    
    #if the new regression model can reduce the BIC 
    if (newBIC < currentBIC){
      bestRegression = newRegression;
      bestRegressionBIC = newBIC;
    }
    else break;
  }
  
  return(list(aic=bestRegressionBIC,reg=bestRegression));
}


#we structure the R code as a C program; this is a choice, not a must
main <- function(datafile)
{
  #read the data
  data = read.table(datafile,header=FALSE);
  
  #the sample size is 148 (number of rows)
  #the explanatory variables are the first 60 columns
  #the last column is the binary response
  response = ncol(data);
  lastPredictor = ncol(data)-1;
  
  #perform a forward "greedy" search for the best logistic regression
  #i.e., the logistic regression with the smallest AIC
  forwardResults = forwardSearchAIC(response,data,lastPredictor);
  
  #perform a backward "greedy" search for the best logistic regression
  backwardResults = backwardSearchAIC(response,data,lastPredictor);
  
  #output the results of our searches
  cat('\n\nForward search gives regression with ',length(forwardResults$reg),'explanatory variables [');
  if(length(forwardResults$reg)>=1)
  {
    for(i in 1:length(forwardResults$reg)) cat(' ',forwardResults$reg[i]);
  }
  cat('] with AIC = ',forwardResults$aic,'\n');
  
  cat('\n\nBackward search gives regression with ',length(backwardResults$reg),'explanatory variables [');
  if(length(backwardResults$reg)>=1)
  {
    for(i in 1:length(backwardResults$reg)) cat(' ',backwardResults$reg[i]);
  }
  cat('] with AIC = ',backwardResults$aic,'\n');
  
  
  # BIC 
  #perform a forward "greedy" search for the best logistic regression
  #i.e., the logistic regression with the smallest BIC
  forwardResults = forwardSearchBIC(response,data,lastPredictor);
  
  #perform a backward "greedy" search for the best logistic regression
  backwardResults = backwardSearchBIC(response,data,lastPredictor);
  
  #output the results of our searches
  cat('\n\nForward search gives regression with ',length(forwardResults$reg),'explanatory variables [');
  if(length(forwardResults$reg)>=1)
  {
    for(i in 1:length(forwardResults$reg)) cat(' ',forwardResults$reg[i]);
  }
  cat('] with BIC = ',forwardResults$aic,'\n');
  
  cat('\n\nBackward search gives regression with ',length(backwardResults$reg),'explanatory variables [');
  if(length(backwardResults$reg)>=1)
  {
    for(i in 1:length(backwardResults$reg)) cat(' ',backwardResults$reg[i]);
  }
  cat('] with BIC = ',backwardResults$aic,'\n');
  
}

# setwd("~/Course/STAT534/STAT534_code/HW2")
main('534binarydata.txt')

# result ----------------------------------------------------------------
# forwardSearch :: The empty logistic regression has AIC =  188.5055 
# backwardSearch :: The full logistic regression has AIC =  122 
# Forward search gives regression with  10 explanatory variables [  23  3  22  21  49  4  20  34  53  46] with AIC =  22 
# Backward search gives regression with  10 explanatory variables [  1  3  9  12  20  23  25  34  41  46] with AIC =  22 
# 
# forwardSearch :: The empty logistic regression has BIC =  191.5027 
# backwardSearch :: The full logistic regression has BIC =  304.8299 
# Forward search gives regression with  10 explanatory variables [  23  3  22  21  49  4  20  34  53  46] with BIC =  54.96934 
# Backward search gives regression with  10 explanatory variables [  1  3  9  12  20  23  25  34  41  46] with BIC =  54.96934 

# Problem 4 ------------------------------------------------------------

# The regression models I got in Problem 2 and Problem 3 are not same, but they have same AIC.
# Same situation for BIC. And the regression models I got from BIC are as same as those got 
# from AIC.
