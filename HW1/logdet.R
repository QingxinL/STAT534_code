logdet <- function(mat1){
  # This funciton computes the logarithm of the determinant of a matrix mat1
  eigen_values <- eigen(mat1)[1]
  log_values <- log(unlist(eigen_values))
  log_det <-sum(log_values)
  return(log_det)
}

logmarglik <- function(data, A){
  
}

mydata = read.table("~/Course/STAT534/STAT534_code/erdata.txt")