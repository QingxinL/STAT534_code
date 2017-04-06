logdet <- function(mat1){
  # This funciton computes the logarithm of the determinant of a matrix mat1
  eigen_values <- eigen(mat1)[1]
  log_values <- log(unlist(eigen_values))
  log_det <-sum(log_values)
  return(log_det)
}

logmarglik <- function(data_in, A){
  # this funciton computes the logarithm of the marginal likelihood
  n <- dim(data_in)[1]
  num_A <- length(A)
  D1 <- matrix(unlist(data_in[1]), ncol = 1, byrow = FALSE)
  D_A <- matrix(unlist(data_in[A]), ncol = num_A, byrow = FALSE)  
  #print(D_A[1:5,])
  D_At <- t(D_A)
  M_A <- diag(num_A) + D_At %*% D_A
  #print(M_A)
  last_term <- 1 + t(D1)%*%D1 - t(D1)%*%D_A%*%solve(M_A)%*%D_At%*%D1
  
  log_marglik = lgamma((n + num_A + 2)/2) - lgamma((num_A + 2)/2) - 0.5 * logdet(M_A) - (n+num_A+2)/2 * log(last_term)
  return(log_marglik)
}

#test the function

m_t = matrix(c(2,1,3,4), ncol=2)
logdet(m_t)

setwd("~/Course/STAT534/STAT534_code/HW1")
mydata = read.table("./erdata.txt") 

logmarglik(mydata,c(2,5,10))



