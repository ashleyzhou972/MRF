
library(stats4)


profile_ll<-function(tau2){
  op2<-mle(logl_2,start=list(alpha =1, eta = 1),fixed = list(tau2 = tau2))
  #summary(op2)
  r = logLik(op2)
  return(list(r,op2))
}

profile_alpha<-function(op2){
  #op2<-mle(logl_2,start=list(alpha =1, eta = 1),fixed = list(tau2 = tau2))
  #summary(op2)
  r = op2@coef[1]
  return(r)
}


profile_eta<-function(op2){
  #op2<-mle(logl_2,start=list(alpha =1, eta = 1),fixed = list(tau2 = tau2))
  #summary(op2)
  r = op2@coef[2]
  return(r)
}



  
  




