## relevant packages

# Implementation of the copula analysis.
install.packages("VineCopula") 
library(VineCopula)

# Implementation of the Anderson-Darling Goodness-of-Fit test.
install.packages("goftest") 
library(goftest)

# Implementation of the (Lilliefors-Corrected) Kolmogorov-Smirnov Goodness-of-Fit test.
install.packages("KScorrect")
library(KScorrect)

# Implementation of the time series analysis.
install.packages("fGarch")
library(fGarch)

library(quantmod)
library(zoo)

install.packages("nortest")
library(nortest)

set.seed(1928)

install.packages("sn")
library(sn)
##############new edit:1) removed net 2) change logic on NA removal: remove rows where either SSE or SPX is empty 

#load data

sse <- read.csv("SSE.csv",header=TRUE)
nk <- read.csv("NIKKEI.csv",header=TRUE)

#prepare data
sse_ret <- sse[,c("Date","Close")]
nk_ret <- nk[,c("Date","Close")]

port <- merge(sse_ret,nk_ret,by="Date",all=FALSE)
names(port) <- c("Date","SSE","NK")
#transform to numerical value
port$SSE <- as.numeric(port$SSE)
port$NK <- as.numeric(port$NK)

w1<-0.5
w2<-0.5

port$Total <- w1*port$NK + w2*port$SSE
head(port)


###remove rows if either SPX or SSE has missing value
#port <- port[!(is.na(port$SPX) | is.na(port$SSE)),] 

port <- port[!(is.na(port$NK) | is.na(port$SSE)),] 


port$SSE_log <- diff(c(NA,log(port$SSE))) #take log return of SSE
port$NK_log <- diff(c(NA,log(port$NK))) #take log return of SSE

port$Total_log <- diff(c(NA,log(port$Total))) #take log return of Total

head(port)
names(port) <- c("Date","SSE_Close","NK_Close","Port_Close","SSE_log","NK_log","port_log")

port <- port[-1, ]  # remove the first row as log would be NA

head(port)
dim(port)

#plot the data
plot(port$NK_log,port$SSE_log)


#Fit the marginal distribution====================================

ret_1 <- port$NK_log
ret_2 <- port$SSE_log

# Building AR models: the Box - Jenkins approach
# Step 1: Identification
par(mfrow=c(1,1))
acf(ret_1, col="green", lwd=2)
pacf(ret_1, col="green", lwd=2,lag.max = 10)  # p = 1, or p = 7
########new edit(Alex): add lag.max to make the plot easier to interpret.
########new edit(Alex): Check AR1 vs AR7, 1) see if it passes the test 2) if so, use AIC 
acf(ret_1^2, col="green", lwd=2)
########we plot the square of log return to check variance. 
########In our case, square ret shows auto correlating, it means the conditional variance varies with time. (meaning we need garch/arch models, although it does not tell which model to use)  

acf(ret_2, col="green", lwd=2)
pacf(ret_2, col="green", lwd=2,lag.max=10)   # p = 4   ########new edit: add lag.max
acf(ret_2^2, col="green", lwd=2)


#Step 2: Estimation

# new edit: (Alex) Fit model1 through a loop (1) order between order  (2) Garch(1,1) to Garch(3,3), p = c(1,2,3),q=c(1,2,3) (3) among different distribution

aic <- list()
dist <- c( "norm","snorm", "sstd", "std","snig", "QMLE","ged","sged") # define distribution

# Iterate over ARMA and GARCH orders
for(d in dist){
  for(arma_order in c(1)){
    for(q in 1:3){
      for(p in 1:3){
        # Dynamically construct the formula string
        formula_string <- paste("~arma(", arma_order, ",0)+garch(", p, ",", q, ")", sep = "")
        # Convert the string to a formula object
        formula_obj <- as.formula(formula_string)
        # fit the GARCH model 
        model <- tryCatch({
          garchFit(formula = formula_obj, data = ret_1, trace = F, cond.dist = d)
        }, error = function(e) {
          return(NA)  # Properly return NA on error
        })
        # define model id 
        model_id <- paste("ARMA", arma_order, "GARCH", p, q, "Distribution", d, sep="_")
        # Store the AIC value if the model was successfully fitted
        if(!is.na(model)){
          aic[[model_id]] <- unname(model@fit$ics["AIC"])  # Corrected to use 'model' object
        } else {
          aic[[model_id]] <- NA
        }
      }
    }
  }
}

# select optimal model with lowest AIC
aic_df <- data.frame(model = names(aic),AIC = unlist(aic),row.names=NULL)
dim(aic_df)
# Find the index of the model with the minimum AIC value
best_model_index <- which.min(aic_df$AIC)
# Extract the best model's identifier based on the index
best_model_id <- aic_df$model[best_model_index] #ARMA_7_GARCH_1_1_Distribution_sstd -4.908382
# Extract the best AIC value
best_aic <- aic_df$AIC[best_model_index]
cat("Best Model:", best_model_id, "with AIC:", best_aic, "\n")

model_1=garchFit(formula=~arma(1,0)+garch(1,1),data=ret_1,trace=F,cond.dist="sstd")
coef(model_1)



###similarly, select best model for model2

aic2 <- list()
dist2 <- c( "norm","snorm", "sstd", "std") # define distribution

# Iterate over ARMA and GARCH orders
for(d in dist2){
  for(arma_order in c(8)){
    for(q in 1:3){
      for(p in 1:3){
        # Dynamically construct the formula string
        formula_string <- paste("~arma(", arma_order, ",0)+garch(", p, ",", q, ")", sep = "")
        # Convert the string to a formula object
        formula_obj <- as.formula(formula_string)
        # fit the GARCH model 
        model <- tryCatch({
          garchFit(formula = formula_obj, data = ret_2, trace = F, cond.dist = d)
        }, error = function(e) {
          return(NA)  # Properly return NA on error
        })
        # define model id 
        model_id <- paste("ARMA", arma_order, "GARCH", p, q, "Distribution", d, sep="_")
        # Store the aic2 value if the model was successfully fitted
        if(!is.na(model)){
          aic2[[model_id]] <- unname(model@fit$ics["AIC"])  # Corrected to use 'model' object
        } else {
          aic2[[model_id]] <- NA
        }
      }
    }
  }
}

# select optimum AIC
aic2_df <- data.frame(model = names(aic2),AIC = unlist(aic2),row.names=NULL)
dim(aic2_df)
# Find the index of the model with the minimum AIC value
best_model_index <- which.min(aic2_df$AIC)
# Extract the best model's identifier based on the index
best_model_identifier2 <- aic2_df$model[best_model_index] #Best Model2: ARMA_8_GARCH_1_1_Distribution_sstd with AIC: -4.277397 
# Extract the best AIC value
best_aic2 <- aic2_df$AIC[best_model_index]
cat("Best Model2:", best_model_identifier2, "with AIC:", best_aic2, "\n")

model_2=garchFit(formula=~arma(8,0)+garch(1,1),data=ret_2,trace=F,cond.dist="sstd")

#################################################
# Step 3: Model checking
# extract residuals for return1 
res_1 <- residuals(model_1, standardize=TRUE )
par(mfrow=c(1,1))
pacf(res_1, col="green", lwd=2)
acf(res_1^2, col="red", lwd=2)
#box test to check residuals 
Box.test(res_1, lag = 10, type = c("Ljung-Box"), fitdf = 1)   ## fitdf = p, where p is from AR(p) 
Box.test(res_1^2, lag = 10, type = c("Ljung-Box"), fitdf = 1)
model_1@fit$ics

nu <- unname(coef(model_1)['shape']) # degree of freedom
xi <- unname(coef(model_1)['skew'])  # skewness


u_1<-psstd(res_1,  mean=0,sd=1,nu=nu,xi=xi)[9:length(ret_1)] 
##new edit (Alex): explain why we throw away the first 8 
#start with 9th, as we use AR(8) in ret_2 , the first 8 is 0  (Tutorial 2 : 55min). keep both u_1 and u_2 same length
#res_1

hist(u_1)



# Further distributional checks using formal tests, in additional to visual checks
#Kolmogorov-Smirnov test
KStest1<-LcKS(u_1, cdf = "punif")
KStest1$p.value
#Anderson-Darling test: this is more conservative test with higher power ( we don't need to run both)
ADtest1<-ad.test(u_1)  
ADtest1$p.value


# returns 2
res_2 <- residuals(model_2, standardize=TRUE)
par(mfrow=c(1,1))
acf(res_2, col="green", lwd=2)
acf(res_2^2, col="red", lwd=2)
Box.test(res_2, lag = 10, type = c("Ljung-Box"), fitdf = 8) # fitdf = p +1
Box.test(res_2^2, lag = 10, type = c("Ljung-Box"), fitdf = 8)#????
model_1@fit$ics


#coef(model_2)
##new edit: although the model is only on AR8, the R package includes AR1, AR2... AR8 (the in between lags). this limiation of R package is not ideal, but we don't lose marks on this. 
#when comes to simulating, we do need to simulate thie structure, we only need to simulate AR8, for example.( Alex)

# degree of freedom nu
# skewness xi

u_2<-psstd(res_2,  mean=0,sd=1,nu=coef(model_2)['shape'],xi=coef(model_2)['skew'])[9:length(ret_2)]
hist(u_2)

# Further distributional checks
#Kolmogorov-Smirnov test
KStest2<-LcKS(u_2, cdf = "punif")
KStest2$p.value
#Anderson-Darling test
ADtest2<-ad.test(u_2)
ADtest2$p.value






# We pass the test for uniformity for both transformed log-returns, so we can proceed 
# to copula modelling.
#======================================================================================================
# Using BiCopSelect function fit various copulas to the dataset and select
# the copula that provides the best fit based on the AIC criterion.
Copulamodel=BiCopSelect(u_1, u_2, familyset=NA, selectioncrit="AIC", indeptest=TRUE, level=0.05,se = TRUE)
Copulamodel
# The model of best fit is Bivariate copula:BB1 (par = 0.2, par2 = 1.06, tau = 0.14)  
#======================================================================================================





# Value-at-Risk uisng Monte Carlo simulation
N=1000
u_simulated=BiCopSim(N, family=Copulamodel$family, Copulamodel$par,  Copulamodel$par2)
# Next we apply the component-wise Inverse Probability Integral Transform (IPIT) to both ret1 and ret2. 

###CHECK THE PARAMETERS!!
coef(model_1)
res1_sim = qsstd(p = u_simulated[,1], nu = coef(model_1)['shape'], mean = coef(model_1)['mu'], sd = coef(model_1)['omega'], xi = coef(model_1)['skew'])
#i think this is the match for sstd, but not 100% sure
coef(model_2)
res2_sim = qsstd(p = u_simulated[,2], nu = coef(model_2)['shape'], mean = coef(model_2)['mu'], sd = coef(model_2)['omega'], xi = coef(model_2)['skew'])


#####################################################################################
# However, "res1_sim" and "res2_sim" are i.i.d.
# So, the next step is to re-introduce  autocorrelation and GARCH effects observed in data

# This section has been omitted as this will be part of your ICA group assignment
#####################################################################################

#first,simulate epsilon-t's from the sstd distribution(the conditional distribution from model 1)
n<-1000
#set the parameter values as in the summary of model_1
mu_1<-coef(model_1)['mu']
ar1_1<-coef(model_1)['ar1']
omega_1<-coef(model_1)['omega']
alpha1_1<-coef(model_1)['alpha1'] 
beta1_1<-coef(model_1)['beta1']
skew_1<-coef(model_1)['skew']
shape_1<-coef(model_1)['shape']  
#simulate the time-varying variance
sigma2s_1<-numeric(n)
y_1<-numeric(n)
sigma2s_1[1]<-omega_1#do we put omega here?
y_1[1]<-sigma2s_1[1]*res1_sim[1]
for (i in 2:n){
  sigma2s_1[i]<-omega_1+alpha1_1*y_1[i-1]^2+beta1_1*sigma2s_1[i-1]
  y_1[i]<-sigma2s_1[i]*res1_sim[i]
}
#simulate the SSE log return
sse_sim<-numeric(n)
sse_sim[1]<-mu_1+ar1_1*port$SSE_log[1161]+y_1[1]
for (i in 2:n){
  sse_sim[i]<-mu_1+ar1_1*sse_sim[i-1]+y_1[i]
}

#plot
plot(sse_sim,type='l',ylab='log return',xlab='time',main='Simulated log return for the SSE',ylim=c(-0.001,0.001))

#second, do the same for Nikkei
#set the parameter values as in the summary of model_1
mu_2<-coef(model_2)['mu']
ar1_2<-coef(model_2)['ar1']
ar2_2<-coef(model_2)['ar2']
ar3_2<-coef(model_2)['ar3']
ar4_2<-coef(model_2)['ar4']
ar5_2<--coef(model_2)['ar5']
ar6_2<-coef(model_2)['ar6']
ar7_2<-coef(model_2)['ar7']
ar8_2<-coef(model_2)['ar8']
omega_2<- coef(model_2)['omega']
alpha1_2<-coef(model_2)['alpha'] 
beta1_2<-coef(model_2)['beta']
skew_2<-coef(model_2)['skew']
shape_2<-coef(model_2)['shape'] 
ar<-c(ar1_2,ar2_2,ar3_2,ar4_2,ar5_2,ar6_2,ar7_2,ar8_2)
#simulate the time-varying variance
sigma2s_2<-numeric(n)
y_2<-numeric(n)
y_2[1]<-sigma2s_2[1]*res2_sim[1]
sigma2s_2[1]<-omega_2
for (i in 2:n){
  sigma2s_2[i]<-omega_2+alpha1_2*y_2[i-1]^2+beta1_2*sigma2s_2[i-1]
  y_2[i]<-sigma2s_1[i]*res2_sim[i]
}
#simulate the Nikkei log return
nk_sim<-numeric(n)
for (i in 1:8){
  nk_sim_i<-0
  for (j in 1:8){
    nk_sim_i<-nk_sim_i+ar[j]*port$NK_log[1161+i-j]
  }
  nk_sim[i]<-mu_2+y_2[i]
}
for (i in 9:1000){
  nk_sim_i<-0
  for (j in 1:8){
    nk_sim_i<-nk_sim_i+ar[j]*nk_sim[i-j]
  }
  nk_sim[i]<-nk_sim_i+mu_2+y_2[i]
}

#plot
plot(nk_sim,type='l',ylab='log return',xlab='time',main='Simulated log return for the Nikkei',ylim=c(-0.15,0.15))









#####################################################################################
portsim <- matrix(0, nrow = N, ncol = 1)
varsim <- matrix(0, nrow = 1, ncol = 2)

portsim=log(1+((exp(res1_sim)-1)+(exp(res2_sim)-1))*(1/2))
varsim=quantile(portsim,c(0.01,0.05))
varsim
