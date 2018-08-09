##################################
#Post MCMC analysis
#including test for eta
#get fitted values
##################################



delete_burn_in <-function(ret, burn_in){
  total_iter = length(ret$alpha)-1
  new_alpha = ret$alpha[burn_in:total_iter]
  new_eta = ret$eta[burn_in:total_iter]
  new_tau2 = ret$tau2[burn_in:total_iter]
  return(list(alpha=new_alpha, eta=new_eta, tau2=new_tau2))
}

get_fitted_y<-function(N, alpha_hat, eta_hat, neighbor){
  muhat = rep(0, N) #initial values is alpha_hat-1 to avoid getting 0
  #A Gibbs process
  iters = 5
  for (t in 1:iters){
    print(t)
    for (i in 1:N) {
      
      muhat[i] = alpha_hat+eta_hat*(neighbor[i,]%*%(muhat-rep(alpha_hat,N)))
    }
  }
  yhat = exp(muhat)
  return(yhat)
}

print_param_estimates<-function(ret){
  alpha_test = t.test(ret$alpha)
  cat(paste(round(mean(ret$alpha),5), " (", round(alpha_test$conf.int[1],5), ", ", round(alpha_test$conf.int[2],5), ") ", sep=""))
  cat("\n")
  eta_test = t.test(ret$eta)
  cat(paste(round(mean(ret$eta),5), " (", round(eta_test$conf.int[1],5), ", ", round(eta_test$conf.int[2],5), ") ", sep=""))
  cat("\n")
  tau2_test = t.test(ret$tau2)
  cat(paste(round(mean(ret$tau2),5), " (", round(tau2_test$conf.int[1],5), ", ", round(tau2_test$conf.int[2],5), ") ", sep=""))
  cat("\n")
  #options(digits=10)
}

goodness_of_fit_test<-function(y, y_hat, bin_num) {
  return(FALSE)
}