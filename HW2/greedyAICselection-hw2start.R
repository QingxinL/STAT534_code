#HOMEWORK 2, PROBLEM 1
#this function uses 'glm' to fit a logistic regression
#and returns the AIC = deviance + 2*NumberOfCoefficients 
getLogisticAIC <- function(response,explanatory,data)
{

}

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
    

  }
  
  return(list(aic=bestRegressionAIC,reg=bestRegression));
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
}

main('534binarydata.txt');