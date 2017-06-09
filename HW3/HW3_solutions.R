# Function to implement MC3 Stochastic search
# over model space of logistic linear regressions
###############################################################################
library(rcdd)

## This is the main function to execute the stochastic search. 
# response: column index corresponding to the response variable
# data: data frame or matrix with all variables
# n_iter: number of search iterations to perform

# returns a list containing the AIC of the best model encountered
# as well as the indicies of variables in the best model encountered.

MC3search = function(response, data, n_iter = 25)
{
	# All possible variables
	vars = setdiff(seq_len(ncol(data)), response)
	
	# iteration zero: get initial model
	starting_fit = getStartingModel(vars, response, data)

	# glm actually gives us the AIC as part of its output.
	cur_AIC = bestModAIC = starting_fit$model$aic
	cur_mod_vars = bestModVars = starting_fit$vars_in

	## get information about the current model's neighbors.
	cur_nbrs = getValidNeighbors(cur_mod_vars, vars, response, data)
	cur_num_nbrs = cur_nbrs$n_validNbr
	
	best_AIC = Inf
	best_AIC_vars = c()
	
	for(r in 1:n_iter)
	{
		# sample a proposal model from current model's neighbors
		ind = sample.int(cur_num_nbrs, 1)
		prop_mod_vars = cur_nbrs$models[[ind]]$vars_nbr
		prop_AIC = cur_nbrs$models[[ind]]$aic
		
		## get information about the proposal model's valid neighbors
		prop_nbrs = getValidNeighbors(prop_mod_vars, vars, response, data)
		prop_num_nbrs = prop_nbrs$n_validNbr
	
		## accept or reject the proposal: do we MOVE?
		move = moveCheck(cur_AIC, cur_num_nbrs, prop_AIC, prop_num_nbrs)
		
		# if we move, we can take all the information about the neighbors
		# of the proposal model and reuse it, to avoid redundant calculation.
		if(move)
		{
			cur_mod_vars = prop_mod_vars
			cur_AIC = prop_AIC
			cur_num_nbrs = prop_num_nbrs
			cur_nbrs = prop_nbrs
		}

		if(cur_AIC < best_AIC)
		{
			best_AIC_vars = cur_mod_vars
			best_AIC = cur_AIC
		}
	}
	cat("Chain finds a minimum AIC of", best_AIC, "\n from model with variables" , 
			sort(best_AIC_vars), "\n \n")
	
	return(list(bestAIC = best_AIC, bestAICvars = sort(best_AIC_vars)))
}

# check whether we accept the new proposal for our chain
moveCheck = function(cur_AIC, cur_num_nbrs, prop_AIC, prop_num_nbrs)
{
	p.cur = -cur_AIC - log(cur_num_nbrs)
	p.new = -prop_AIC - log(prop_num_nbrs)
	
	if(p.new > p.cur)
		#if this happens we accept the move automatically
		return(TRUE)
	else
	return(log(runif(1)) < p.new - p.cur)
}

# randomly sample model size, then model predictors.
# vars: list of indices of explanatory variables
initModelSample = function(vars)
{
	# sample number of predictors in model
	k = sample(0:length(vars), 1)
	
	# sample variables in model of size k
	sample(vars, k, replace = FALSE)
}

# rather than recompute the fit for the model, pass the original
# output from the function that calls it.
isValidLogisticRCDD <- function(fit, response, data)
{
	if(fit$nvar == 1)
	{
		return(TRUE)
	}

	tanv = fit$x
	tanv[data[,response] == 1, ] <- (-tanv[data[,response] == 1, ])
	vrep = cbind(0, 0, tanv)
	
	#with exact arithmetic; takes a long time
	#lout = linearity(d2q(vrep), rep = "V")
	
	lout = linearity(vrep, rep = "V")
	return(length(lout)==nrow(data))
}

# with this function we randomly choose starting models until we get a 
# valid place to start. Return the starting model and the variable list that
# fit it.

# vars: vector of indicies of explanatory variables
# data: matrix or data frame of all variables
# response: index of response variable
getStartingModel = function(vars, response, data)
{
	# keep sampling models until we get a valid fit.
	# technically we might worry that all models will fail,
	# though the empty model will always work unless the response
	# values are all one or all zero, an unlikely case.
	while(1)
	{
		#sample from pool of valid models
		vars_in = initModelSample(vars)

		#fit the sampled model
		if(length(vars_in) == 0)
		{		
			vars_in = NULL
			fit = glm(data[,response] ~ 1, family=binomial(link=logit),
					x = TRUE)
			fit$nvar = 1
		}
		else
		{
			fit = glm(data[,response] ~ as.matrix(data[,as.numeric(vars_in)]), 
					family=binomial(link=logit), x = TRUE)
			fit$nvar = length(vars_in) + 1
		}
		
		#check whether the sampled model is valid.		
		if(isValidLogisticRCDD(fit, response, data))
		{
			return(list(model = fit, vars_in = vars_in))
		}
	}	
}

# from all neighbors, returning number of valid neighbors, a list
# containing the neighbor model variables, and a list containing the
# fits of all neighbor model variables.
getValidNeighbors = function(vars_in, all_vars, response, data)
{
	# compute all regressions for neighbors. Neighbors either add one
	# variable not in the model or remove a variable that is in the model. 
	regressions = lapply(all_vars, function(i)
			{
				#for each variable "i", either remove it or add it depending on whether
				#it is already in the model.
				test = i %in% vars_in
				if(test)
					vars_nbr = setdiff(vars_in, i)
				else
					vars_nbr = c(vars_in,i)
				
				if(length(vars_nbr) == 0)
				{
					fit = glm(data[,response] ~ 1, family = binomial, x = TRUE)
					fit$nvar = 1
					fit$vars_nbr = vars_nbr
				}
				else
				{
					fit = glm(data[,response] ~ as.matrix(data[,vars_nbr]), 
							family = binomial, x = TRUE)
					fit$nvar = length(vars_nbr) + 1
					fit$vars_nbr = vars_nbr
				}
				return(fit)
			})
	
	# validate each neighbor regression with RCDD
	valid.check = vector(length = length(all_vars))
	for(i in seq_len(length(all_vars)))
	{
		valid.check[i] = isValidLogisticRCDD(regressions[[i]], response, data)
	}
	# get number of valid neighbors
	n_validNbr = sum(as.numeric(valid.check))
	
	# return valid neighbors, their fits, and the number of them.
	return(list(models = regressions[valid.check], 
					n_validNbr = n_validNbr))
}

set.seed(123)
data = read.table("534binarydatasmall.txt", header = FALSE)
nChain = 10
results = vector("list", length = nChain)
for(i in 1:nChain)
{
	results[[i]] = MC3search(response = 11, data, n_iter = 25)
}


#####################
# Problem 2 comments, in brief.
#####################

# A few of the chains arrived at the same model, but 
# clearly there is quite a lot of variability in the model output.
# Given that this is a stochastic search, this is an expected outcome
# for chains run only 25 iterations. See sample output below:

#Chain finds a minimum AIC of 117.3657 
#from model with variables 1 2 3 5 6 7 8 10 

#Chain finds a minimum AIC of 133.9951 
#from model with variables 1 2 3 8 9 

#Chain finds a minimum AIC of 137.4509 
#from model with variables 1 2 4 5 7 9 10 

#Chain finds a minimum AIC of 111.2943 
#from model with variables 1 2 3 5 6 7 8 9 10 

#Chain finds a minimum AIC of 133.3354 
#from model with variables 2 3 6 7 8 9 10 

#Chain finds a minimum AIC of 119.1088 
#from model with variables 1 3 4 5 6 7 8 9 10 

#Chain finds a minimum AIC of 112.8261 
#from model with variables 1 2 3 5 7 8 9 10 

#Chain finds a minimum AIC of 111.2943 
#from model with variables 1 2 3 5 6 7 8 9 10 

#Chain finds a minimum AIC of 111.4978 
#from model with variables 1 2 3 6 7 8 9 10 

#Chain finds a minimum AIC of 111.2943 
#from model with variables 1 2 3 5 6 7 8 9 10 


