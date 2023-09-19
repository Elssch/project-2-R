library(data.table)

path <- readline(prompt = "Give the path to file: ")
data <- fread(path)

select <- as.numeric(readline(prompt = "What element do you need:\n1 - mixed-effects model,\n2 - parameters of variance,\n3 - both elements: "))

y_col = readline(prompt = "Give the explained variable (column name): ")
X_col = readline(prompt = "Give the fixed effects: ")
Z_col = readline(prompt = "Give the random effects: ")
A_col = readline(prompt = "Give the covariances of random effects: ")
sigma_a <- as.numeric.readline(prompt = "Podaj wartość sigma_a: ")
sigma_e <- as.numeric.readline(prompt = "Podaj wartość sigma_e: ")

y <- as.matrix(data[[y_col]])
X <- as.matrix(data[[X_col]])
Z <- as.matrix(data[[Z_col]])
A <- as.matrix(data[[A_col]])

mme <- function(y, X, Z, A, sigma_a, sigma_e) {
  G1 <- ginv(A) * (sigma_e / sigma_a)
  C <- rbind(cbind(t(X) %*% X, t(X) %*% Z),
             cbind(t(Z) %*% X, t(Z) %*% Z + G1))
  r <- rbind(t(X) %*% y, t(Z) %*% y)
  mme <- ginv(C) %*% r
  list(C = C, mme = mme)
}
mme_result <- mme(y, X, Z, A, sigma_a, sigma_e)

C = as.matrix(mme_result$C)
mme_est = mme_result$mme

num_ef_los <- nrow(Z)
indeks_los <- (ncol(X) + 1):(ncol(X) + num_ef_los)
C22 <- C[indeks_los, indeks_los]
r2 = diag(1 - C22 * 2)
r = sqrt(r2)

if (opcja == 1 || opcja == 3) {
  print("Estimators of fixed effects:")
  print(mme_est[1:(ncol(X))])
  print("Estimators of random effects:")
  print(mme_est[ncol(X) + 1):length(mme_est)])
}

if (select == 2 || select == 3) {
  EM = function(y, X, Z, A, sigma_a, sigma_e) {
    n = nrow(X)
    p = ncol(X) 
    q = nrow(A) 
    
    t = 1 
    tmp = 0.1 
    
    while (tmp > 0.00001) {
      mme_new = mme(y, X, Z, A, sigma_a, sigma_e)
      C_new = ginv(mme_new$C)
      Ck = C_new[(p+1):(p+q), (p+1):(p+q)]
      mme2 = mme_new$mme
      
      a = as.matrix(mme2[(p+1):(p+q)])
      sigma_a_new = (t(a)%*%ginv(A)%*%a + sum(diag(ginv(A)%*%Ck))*c(sigma_e))/q
      
      res = as.matrix(y-X%*%as.matrix(mme2[1:p]) - Z%*%as.matrix(mme2[(p+1):(p+q)]))
      X.tmp1 = cbind(X,Z) %*% C_new
      X.tmp2 = t(cbind(X,Z))
      sigma_e_new = (t(res)%*%res + sum(diag(X.tmp1%*%X.tmp2))*c(sigma_e))/n
      
      tmp = max(abs(sigma_a - sigma_a_new), abs(sigma_e - sigma_e_new))
      sigma_a = sigma_a_new
      sigma_e = sigma_e_new
      
      t = t + 1
    }
    list(t = t, sigma_a = sigma_a, sigma_e = sigma_e)
  }
  result = EM(y, X, Z, A, sigma_a, sigma_e)
  print("Estimators of parametres of variance")
  print(result$sigma_a+result$sigma_e)
}