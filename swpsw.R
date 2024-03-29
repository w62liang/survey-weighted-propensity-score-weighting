################################################################
#           Algorithms for the paper titled                    # 
# "Propensity Score Weighting with Post-Treatment Survey Data" #
################################################################
# Author: Wei Liang and Changbao Wu                            #
# Date: 2023-09-24                                             #
################################################################


### Survey-weighted covariate balancing score function

#-------------------inputs----------------------
# subsampleat: indices of the treated individuals
# subsampleac: indices of the control individuals
# x: initial value
# weights: survey weights
# X: design matrix of covariates
# the values of a and b rely on the estimand
#-------------------------------------------------

score_bal <- function(x, subsamplet, subsamplec, weights, X, a = -1, b = -1){
  
  temp <- drop(1/(1+exp(X%*%(-x))))
  return(drop(t(X[subsamplet,])%*%drop(weights[subsamplet]*(temp^a*(1-temp)^(b+1))[subsamplet])) - 
         drop(t(X[subsamplec,])%*%drop(weights[subsamplec]*(temp^(a+1)*(1-temp)^b)[subsamplec])))
}



### survey-weighted propensity score weighting

#-----------------------------------inputs------------------------------------
# swt: survey weights
# dat: a dataframe which includes the response, the treatment variable and  
# the covariates to be adjusted.
# x_name: names of the covariates.
# y_name: response's name.
# t_name: treatment's name; the treatment must be binary and takes values 
# 0 and 1.
# subgroup: indicators for the sub-population, e.g., (1,1,0,...,1,0).
# estimand: SFATE, SFATT, or SFATO.
# method: "GLM" or "CBPS".
#-----------------------------------------------------------------------------

#------------------------------------outputs----------------------------------
# estimate: point estimate of the finite-population weighted average treatment
# effect.
# std.er: estimated standard error of the point estimator.
# wtnl1 and wtnl0: estimated and normalized weights for re-weighting the 
# outcomes.
#------------------------------------------------------------------------------

svypsw <- function(swt, dat, x_name, y_name, t_name, estimand = "SFATE", 
                       method = "CBPS", subgroup){
  
  # estimated population size and sample size
  Nb <- sum(swt)
  n <- length(swt)
  
  # response and treatment variables
  y <- dat[,y_name]
  treat <- (dat[,t_name]==1)
  subtreat <- (treat*subgroup)==1
  subcontrol <- ((1-treat)*subgroup)==1
  samplet <- (1:n)[treat]
  samplec <- (1:n)[!treat]
  subsamplet <- (1:n)[subtreat]
  subsamplec <- (1:n)[subcontrol]
  
  # design matrix
  X <- as.matrix(data.frame(model.matrix(~., data = dat[,x_name])))
  
  # estimands
  if(estimand == "SFATE"){
    a <- -1
    b <- -1
  }
  if(estimand == "SFATT"){
    a <- 0
    b <- -1
  }
  if(estimand == "SFATO"){
    a <- 0
    b <- 0
  }
  
  # logistic regression
  forml <- as.formula(paste(t_name,"~", paste(x_name, collapse = "+")))
  fit <- glm(forml, family = binomial(link = "logit"), data = dat, weights = swt/Nb)
  alpha_ml <- fit$coefficients
  
  # estimated propensity score
  if(method == "GLM"){
  e_fit <- fit$fitted.values
  }
  
  if(method == "CBPS"){
  alpha_bal <- nleqslv(fn = score_bal, subsamplet = subsamplet, 
                         subsamplec = subsamplec, weights = swt/Nb,
                         X = X, x = alpha_ml,a = a, b = b)$x
  e_fit <- drop(1/(1+exp(X%*%(-alpha_bal))))
  }
  
  # estimated weights
  ps <- e_fit
  wt1 <- subtreat*ps^a*(1-ps)^(b+1)
  wt0 <- subcontrol*ps^(a+1)*(1-ps)^b
  wtnl1 <- wt1*swt/sum(wt1*swt)
  wtnl0 <- wt0*swt/sum(wt0*swt)
  
  # estimates of the treatment effects
  est1 <- drop(y%*%wtnl1)
  est0 <- drop(y%*%wtnl0)
  est <- est1 - est0
  
  ## variance estimation 
  wt <- wt1 - wt0
  
  if(method=="GLM"){
    wtps <- treat - ps
    fdr_ps <- -ps*(1-ps)
  }
  if(method=="CBPS"){
    wtps <- wt
      if((a == -1)&&(b == -1)){
        fdr_ps <- -(1-ps)/ps*subtreat - ps/(1-ps)*subcontrol
      }
      if((a == 0)&&(b == -1))
      {
        fdr_ps <- -ps/(1-ps)*subcontrol
      }
      if((a == 0)&&(b == 0)){
        fdr_ps <- -ps*(1-ps)*(subgroup==1)
      }
  }
  
  pia_var <- swt*(swt - 1)
  
  # naive estimating equation
  psi <- rbind(t(wtps*X),
               ((y - est1)*wt1), 
               ((y - est0)*wt0), 0)
  
  len_psi <- length(psi[,1])
  
  # covariance of psi
  cov_psi <- psi%*%(pia_var*t(psi))/Nb^2
  
  # partial derivative of psi with respect to gamma
  if((a == -1)&&(b == -1)){
    ppsi_gamma <- rbind(t(fdr_ps*swt*X)%*%X,
                        drop(((-(1-ps)/ps)*(y-est1)*swt)[subsamplet]%*%X[subsamplet,]),
                        drop(((ps/(1-ps))*(y-est0)*swt)[subsamplec]%*%X[subsamplec,]),
                        0)/Nb
  }
  if((a == 0)&&(b == -1)){
    ppsi_gamma <- rbind(t(fdr_ps*swt*X)%*%X,0,
                        drop(((ps/(1-ps))*(y-est0)*swt)[subsamplec]%*%X[subsamplec,]),
                        0)/Nb
  }
  if((a == 0)&&(b == 0)){
    ppsi_gamma <- rbind(t(fdr_ps*swt*X)%*%X,
                        drop(((-(1-ps)*ps)*(y-est1)*swt)[subsamplet]%*%X[subsamplet,]),
                        drop(((ps*(1-ps))*(y-est0)*swt)[subsamplec]%*%X[subsamplec,]),
                        0)/Nb
  }
  
  # partial derivative of psi with respect to tau1, tau0, and tau
  ppsi_tau1 <- rep(0,len_psi)
  ppsi_tau0 <- ppsi_tau1
  ppsi_tau <- ppsi_tau1
  ppsi_tau1[len_psi-2] <- -sum(wt1*swt)/Nb
  ppsi_tau0[len_psi-1] <- -sum(wt0*swt)/Nb
  ppsi_tau1[len_psi] <- -1
  ppsi_tau0[len_psi] <- 1
  ppsi_tau[len_psi] <- 1
  
  # partial derivative of psi with respect to eta
  ppsi_eta <- cbind(ppsi_gamma,ppsi_tau1,ppsi_tau0,ppsi_tau)
  ppsi_eta_inv <- solve(ppsi_eta)
  
  # covariance of eta_hat
  cov_eta <- (ppsi_eta_inv)%*%(cov_psi)%*%t(ppsi_eta_inv)
  
  # variance of tau_hat
  var_tau <- cov_eta[len_psi,len_psi]
  
  return(list(estimate = est, std.er = sqrt(var_tau), wtnl1 = wtnl1, wtnl0 = wtnl0))
}



#----------- a real data example ----------------#
#----------------------------------------------------------------------------------------------
# The the Medical Expenditure Panel Survey data is publicly available at 
# https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-224
#----------------------------------------------------------------------------------------------

# imputation function
impute_meidan <- function(x){
  x[is.na(x)] <- median(x, na.rm=TRUE)
  return(x)
}

# load packages
library(nleqslv)

# load data
h224.na <- read.csv("h224_na.csv", row.names = 1)

# impute missing values with the median
dat_h224 <- data.frame(apply(h224.na, MARGIN = 2, FUN = impute_meidan))

# factorize discrete variables
continuous_variables <- c("mental_score","physical_score","BMI","age",
                         "health_care_exp","personal_weights","race")
discrete_variables <- colnames(dat_h224)[!colnames(dat_h224)%in%continuous_variables]
dat_h224[,discrete_variables] <- lapply(dat_h224[,discrete_variables], factor)

# race
#------------------
# 1 - hispanic
# 2 - white
# 3 - black
# 4 - asian
# 5 - other
#------------------

# pick up hispanics and non-hispanic whites
dat_ <- dat_h224[((dat_h224$race==1)|(dat_h224$race==2)),]
dat_$race[dat_$race==1] <- 0
dat_$race[dat_$race==2] <- 1

# sub-population: 18<=age<=55
ind_sub <- (dat_$age>=18)&(dat_$age<=55)

# estimate the SFATE by survey-weghited GLM method
xname <- colnames(dat_)[!colnames(dat_)%in%c("health_care_exp",
                                             "personal_weights","race")]

swglm_ate <- svypsw(swt = dat_$personal_weights, dat = dat_, x_name = xname,
       y_name = "health_care_exp", t_name = "race", method = "GLM",
       estimand = "SFATE", subgroup = ind_sub)

cat("estimate:",swglm_ate$estimate," standard error:", swglm_ate$std.er)

# estimate the SFATT by survey-weghited CBPS method
swcbps_att <- svypsw(swt = dat_$personal_weights, dat = dat_, x_name = xname,
                y_name = "health_care_exp", t_name = "race", method = "CBPS",
                estimand = "SFATT", subgroup = ind_sub)

cat("estimate:",swcbps_att$estimate," standard error:", swcbps_att$std.er)

