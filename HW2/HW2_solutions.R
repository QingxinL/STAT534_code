
# compute the BIC for a single model, given preditor indices vars, response 
# index y, and a data frame "data"

# returns the model BIC, or "NA" when glm fails to converge.

getLogisticBIC <- function(y, vars, data)
{	
	# check if the regression has no explanatory variables
	if(0 == length(vars))
	{
		out = glm(data[,y] ~ 1, family = binomial);
	}
	
	# regression with at least one explanatory variable
	else
	{
		out = suppressWarnings(glm(data[,y] ~ as.matrix(data[,as.numeric(vars)]),
				family = binomial));
	}
	
	# compute the BIC
	BIC = out$deviance+log(nrow(data))*(1+length(vars))

	# we must check whether the glm algorithm properly converged,
	# as judged by the algorithm.
	
	converged = out$converged
	
	# when glm fails to converge return NA
	return(ifelse(converged, BIC, NA))
}

##############################################################
####### PROBLEM 1
##############################################################

# compute the AIC for a single model.

getLogisticAIC <- function(y, vars, data)
{
	# check if the regression has no explanatory variables
	if(0 == length(vars))
	{
		out = glm(data[,y] ~ 1, family = binomial);
	}
	
	# regression with at least one explanatory variable. We suppress
	# convergence warnings to reduce screen clutter. 
	else
	{
		out = suppressWarnings(glm(data[,y] ~ as.matrix(data[,as.numeric(vars)]),
				family = binomial));
	}
	
	# compute the AIC
	AIC = out$deviance+2*(1+length(vars))
	
	# we must check whether the glm algorithm properly converged,
	# as judged by the algorithm.
	
	converged = out$converged
	
	# when glm fails to converge return NA, otherwise return AIC
	return(ifelse(converged, AIC, NA))
}


#####################################################
# PROBLEM 2
#####################################################

# "addOneVarAIC" takes the indicies of variables both in, and not in the model
# For each variable not in the model, it is added to the current variable pool
# and the AIC is evaluated for this new model. 

# Returns a vector of these added-variable AIC values
addOneVarAIC = function(varsNotIn, data, y, varsIn)
{
	AIC = vector(length = length(varsNotIn))
	
	# cycle through variables not in current model, computing AIC 
	for(i in 1:length(varsNotIn))
	{
		AIC[i] = getLogisticAIC(y, var = c(varsNotIn[i], varsIn), data)
	}
	return(AIC)
}

#### MAIN FORWARD AIC FUNCTION

# takes a data frame, and the index corresponding to the response.
# returns a list with element "vars" giving indices of selected variables
# and "AIC" giving AIC of forward-selected model.

forwardAIC = function(data, y)
{
	varsIn = NULL
	curBestAIC = getLogisticAIC(y, vars = varsIn, data)
	
	# check whether we manage to start the algorithm
	if(is.na(curBestAIC))
	{
		cat("Failed to fit null model \n")
		return(NULL)
	}
	
	# move through model sizes 1 through p, if needed.
	num.pred = 1
	while(num.pred <= ncol(data) - 1)
	{
		# find the variables not in the model
		varsNotIn = setdiff(1:ncol(data), c(y, varsIn))
		
		# get the updated AICs
		updateAIC = addOneVarAIC(varsNotIn, data, y, varsIn)
		
		# check whether any of the updated models have proper AIC values
		if(all(is.na(updateAIC)))
		{
			cat("All added-variable model fits failed to converge.\n",
					"Returning best AIC model prior to iteration with total convergence failure: \n")
			break
		}
		
		# check whether we improved model fit. If we failed, we exit the loop
		# and return the best model.
		if(min(updateAIC, na.rm = TRUE) > curBestAIC)
		{
			break
		}
		
		# update the model variables and the current AIC
		varsIn = c(varsIn, varsNotIn[which.min(updateAIC)])
		curBestAIC = min(updateAIC, na.rm = TRUE)
		num.pred = num.pred + 1
	}
	
	cat("The set of variables selected by forward AIC selection is \n", sort(varsIn),
			"\n with AIC = ",curBestAIC,"\n")
	return(list(AIC = curBestAIC, vars = sort(varsIn)))
}

############################################################
#### PROBLEM 3
###########################################################

# "subtractOneVarAIC" takes the data, and the indicies of the current values 
# in the model. Model variables are removed one at a time and the AIC evaluated
# for the reduced models.

# returns a vector of these AIC values
subtractOneVarAIC = function(data, y, varsIn)
{
	AIC = vector(length = length(varsIn))
	for(i in 1:length(varsIn))
	{
		# update the explanatory variables
		AIC[i] = getLogisticAIC(y, vars = varsIn[-i], data)
	}
	return(AIC)
}


#### MAIN BACKWARD AIC FUNCTION
# takes a data frame, and the index corresponding to the response.
# returns a list with element "vars" giving indices of selected variables
# and "AIC" giving AIC of backward-selected model.

backwardAIC = function(data, y)
{
	# total number of variables
	p = ncol(data)
	
	#start with all explanatory variables
	varsIn = setdiff(1:p, y)
	
	# get AIC for full model. If this model fails to converge, we return NULL
	# and give a handy suggestion.
	curBestAIC = getLogisticAIC(y, vars = varsIn , data)
	if(is.na(curBestAIC))
	{
		cat("Failed to fit full model. Exiting algorithm\n")
		return(NULL)
	}

	# move through model sizes
	for(i in 1:(p-1))
	{
		# get new AIC values for all possible 1-removed models and we
		# check whether convergence failed for any of them.
		
		updateAIC = subtractOneVarAIC(data, y, varsIn)
		if(all(is.na(updateAIC)))
		{
			cat("All reduced-variable model fits failed to converge.\n",
					"Returning best AIC model prior to iteration with total convergence failure: \n")
			break
		}
		
		# check whether we failed to improve model fit. If so, we break.
		if(min(updateAIC, na.rm = TRUE) > curBestAIC)
		{
			break
		}
		
		# take out variable corresponding to largest decrease in AIC, and
		# move to next iteration
		varsIn = setdiff(varsIn, which.min(updateAIC))
		curBestAIC = min(updateAIC, na.rm = TRUE)
	}
	
	# print the best vector (perhaps conditional on convergence, in a sense)
	cat("The set of variables selected by backward AIC selection is \n", sort(varsIn),
			"\n with AIC = ", curBestAIC, "\n")
	
	return(list(AIC = curBestAIC, vars = varsIn))
}

####################################################################################
#### PROBLEM 4
####################################################################################

# This function takes the indicies of variables both in, and not in the model
# For each variable not in the model, it is added to the current variable pool and the BIC
# is evaluated for this new model. Returns a vector of these BIC values
addOneVarBIC = function(varsNotIn, data, y, varsIn)
{
	BIC = vector(length = length(varsNotIn))
	for(i in 1:length(varsNotIn))
	{
		BIC[i] = getLogisticBIC(y, vars = c(varsNotIn[i], varsIn), data)
	}
	return(BIC)
}

# This function takes the data, and the indicies of the current values in
# the model. Each variable in the model is removed and the BIC evaluated
# for the reduced model.

# returns a vector of these BIC values
subtractOneVarBIC = function(data, y, varsIn)
{
	BIC = vector(length = length(varsIn))
	for(i in 1:length(varsIn))
	{
		# update the explanatory variables
		BIC[i] = getLogisticBIC(y, vars = varsIn[-i], data)
	}
	return(BIC)
}

#### MAIN BIC BACKWARD ELIMINATION FUNCTION

#essentially the same as backwardAIC
backwardBIC = function(data, y)
{
	# total number of variables
	p = ncol(data)
	
	#start with all explanatory variables
	varsIn = setdiff(1:p, y)
	
	# get BIC for full model. If this model fails to converge, we return NULL
	# and give a handy suggestion.
	curBestBIC = getLogisticBIC(y, vars = varsIn , data)
	if(is.na(curBestBIC))
	{
		cat("Failed to fit full model. Exiting algorithm\n")
		return(NULL)
	}
	
	# move through model sizes
	for(i in 1:(p-1))
	{
		# get new BIC values for all possible 1-removed models and we
		# check whether convergence failed for any of them.
		
		updateBIC = subtractOneVarBIC(data, y, varsIn)
		if(all(is.na(updateBIC)))
		{
			cat("All reduced-variable model fits failed to converge.\n",
					"Returning best BIC model prior to iteration with total convergence failure: \n")
			break
		}
		
		# check whether we failed to improve model fit. If so, we break.
		if(min(updateBIC, na.rm = TRUE) > curBestBIC)
		{
			break
		}
		
		# take out variable corresponding to largest decrease in BIC, and
		# move to next iteration
		varsIn = setdiff(varsIn, which.min(updateBIC))
		curBestBIC = min(updateBIC, na.rm = TRUE)
	}
	
	# print the best vector (perhaps conditional on convergence, in a sense)
	cat("The set of variables selected by backward BIC selection is \n", sort(varsIn),
			"\n with BIC = ", curBestBIC, "\n")
	
	return(list(BIC = curBestBIC, vars = varsIn))
}

#### MAIN FORWARD BIC FUNCTION
forwardBIC = function(data, y)
{
	varsIn = NULL
	curBestBIC = getLogisticBIC(y, vars = varsIn, data)
	
	# check whether we manage to start the algorithm
	if(is.na(curBestBIC))
	{
		cat("Failed to fit null model \n")
		return(NULL)
	}
	
	# move through model sizes 1 through p, if needed.
	num.pred = 1
	while(num.pred <= ncol(data) - 1)
	{
		# find the variables not in the model
		varsNotIn = setdiff(1:ncol(data), c(y, varsIn))
		
		# get the updated BICs
		updateBIC = addOneVarBIC(varsNotIn, data, y, varsIn)
		
		# check whether any of the updated models have proper BIC values
		if(all(is.na(updateBIC)))
		{
			cat("All added-variable model fits failed to converge.\n",
					"Returning best BIC model prior to iteration with total convergence failure: \n")
			break
		}
		
		# check whether we improved model fit. If we failed, we exit the loop
		# and return the best model.
		if(min(updateBIC, na.rm = TRUE) > curBestBIC)
		{
			break
		}
		
		# update the model variables and the current BIC
		varsIn = c(varsIn, varsNotIn[which.min(updateBIC)])
		curBestBIC = min(updateBIC, na.rm = TRUE)
		num.pred = num.pred + 1
	}
	
	cat("The set of variables selected by forward BIC selection is \n", sort(varsIn),
			"\n with BIC = ",curBestBIC,"\n")
	return(list(BIC = curBestBIC, vars = sort(varsIn)))
}


data = read.table("534binarydata.txt", header = FALSE)
 
forwardAIC(data, 61)
backwardAIC(data, 61)
forwardBIC(data, 61)
backwardBIC(data, 61)

#################################################################################
# With this exercise we learn not to always trust the output from functions
# that we rely on. When "glm" fails to converge, we have no guarantee that 
# its output is sensible. In order to avoid dealing with junky numbers, we put checks
# in our functions so that they recognize when glm fails. We do something
# reasonable when "glm" does fail, which is as good as we can get without
# resorting to linear programming like that seen in Geyer's paper. 

# The backward AIC and BIC algorithms abort due to convergence failure
# in "glm". The convergence failure in the full AIC and BIC models means that 
# our backwards elimination algorithms will not produce meaningful output, hence
# we cannot compare the results with the forward routine. The two forward 
# algorithms return the same model, with explanatory variables

# 3 4 20 21 22 23 34 49 50 53
#################################################################################



