gramschmidt <- function(x) {
  # Get the number of rows and columns of the matrix
  n <- ncol(x)
  m <- nrow(x)
  
  # Initialize matrices Q and R
  Q <- matrix(0, m, n)
  R <- matrix(0, n, n)
  
  # Take care of the first one (base case)
  u1 <- x[,1]
  e1 <- u1/sqrt(sum(u1^2))
  
  Q[,1] <- e1
  R[1,1] <- t(x[,1]) %*% e1
    
  # Gram-Schmidt process
  # For each column (j) from 2 to n
  for (j in 2:n) {
    # temp store vector we need from x
    aj = x[,j]
    
    # For each row entry (i) from 1:j in this column j
    for (i in 1:j) {
      
      # Define each ij entry of R as the product of Qi and our vector from x
      R[i,j] <- t(aj) %*% Q[,i]  
      
      # Now we subtract from the x vector (aj) the dot product of the (product 
      # of Qi and the x vector (equivalent to Rij)) times (Qi)
      aj <- aj - (R[i,j] * Q[,i])

      
      }      
    
    # On the diagonal entry Rjj we want product
    # R[j,j] <- sqrt(sum(uj^2))
    # Q[,j] <- uj / R[j,j]
   
    Q[,j] <- aj / sqrt(sum(aj^2))
    R[j,j] <- t(x[,j] )%*% Q[,j]
  }
  
  # Return matrices Q and R in a list
  QRdecomp <- list('Q'=Q, 'R'=R)
  return(QRdecomp)
}




A <- matrix(sample(10,9), 3,3)
A

gramschmidt(A)
