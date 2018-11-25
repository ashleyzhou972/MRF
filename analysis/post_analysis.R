##################################
#Post MCMC analysis
#including test for eta
#get fitted values
##################################



delete_burn_in <-function(ret, burn_in, neighbor){
  total_iter = length(ret$alpha)
  new_w = ret$w[,(burn_in:total_iter)]
  new_alpha = ret$alpha[burn_in:total_iter]
  new_eta = ret$eta[burn_in:total_iter]
  new_tau2 = ret$tau2[burn_in:total_iter]
  return(list(w=new_w, alpha=new_alpha, eta=new_eta, tau2=new_tau2, neighbor=neighbor))
}

get_fitted_y<-function(ret){
  w=ret$w
  yhat=rowMeans(exp(w))
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

RMSE<-function(y_ob, y_fitted) {
  naid = c(which(is.na(y_ob)), which(is.na(y_fitted)))
  rmse=sqrt(mean((y_ob-y_fitted)^2, na.rm=T))
  return(rmse)
}

# RRSE<-function(y_ob, y_fitted, ret) {
#   #ret should be deleted of burn-in
#   fitted_mean = mean(exp(ret$alpha))
#   rrse = sqrt(sum((y_ob-y_fitted)^2,na.rm=T)/sum((y_fitted-fitted_mean)^2,na.rm=T))
#   return(rrse)
# }

RRSE1<-function(y_ob, y_fitted, alpha) {
  #ret should be deleted of burn-in
  fitted_mean = mean(exp(alpha))
  rrse = sqrt(sum((y_ob-y_fitted)^2,na.rm=T)/sum((y_fitted-fitted_mean)^2,na.rm=T))
  return(rrse)
}

RRSE2<-function(y_ob, y_fitted) {
  rrse = sqrt(sum((y_ob-y_fitted)^2,na.rm=T)/sum((y_ob-mean(y_ob, na.rm=T))^2,na.rm=T))
  return(rrse)
}

