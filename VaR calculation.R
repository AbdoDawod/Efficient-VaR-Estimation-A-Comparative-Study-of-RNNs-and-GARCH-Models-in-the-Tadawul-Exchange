pkges <- c('scales','rugarch', 'xts', 'ggplot2','xtable', 'psych','quantmod', 'PerformanceAnalytics','cowplot','urca','FinTS','forecast' )[-c(5, 8:11)]
sapply(pkges, require, character.only = TRUE) 

getSymbols("^TASI.SR", src = "yahoo", from=as.Date("2007-03-19"),to=as.Date("2024-07-27"))

r_t <- na.omit(CalculateReturns(na.omit(TASI.SR[,6]), method="log"))
rm(TASI.SR)
#ts.plot(r_t$TASI.SR.Adjusted)

cl <- c(0.9, 0.95, 0.99, 0.995)

HS_VaR <- sapply(cl,function(alpha) VaR(as.vector(r_t), p= alpha, method="historical")[1])
HS_ES <- sapply(cl,function(alpha) ES(as.vector(r_t), p= alpha, method="historical")[1])

G_VaR <-  sapply(cl,function(alpha) VaR(as.vector(r_t), p= alpha, method="gaussian")[1])
G_ES <-  sapply(cl,function(alpha) ES(as.vector(r_t), p= alpha, method="gaussian")[1])

Cornish_FisherVaR <- sapply(cl, function(alpha) VaR(as.vector(r_t), p= alpha, method="modified")[1])
Cornish_FisherES <- sapply(cl, function(alpha) ES(as.vector(r_t), p= alpha, method="modified")[1])

print(percent( mean(r_t)+qnorm((1-cl))*sd(r_t) , accuracy = 0.001)) #VaR
print(percent( mean(r_t) + dnorm(qnorm((1-cl)))*sd(r_t)/(1-cl) , accuracy = 0.001)) # ES

##### MC
source("H:/MegSync/Unideb/My Research/Publications/IEEE conference/VaR Monte Carlo.R")
MC_VaR <- VaR.MC(r_t,Wo=1,n=1000,alpha=1-cl,k=1)
MC_ES <- CVaR.MC2(r_t,Wo=1,n=1000,alpha=1-cl,k=1)


##### GARCH VaR
#GARCH Spec - (Change Distribution here)
gspec11 <- ugarchspec(variance.model = list(model = "eGARCH", 
                                            garchOrder = c(1, 1)),
                      mean.model=list(armaOrder=c(0,0)), 
                      distribution="sstd")

#Rolling Estimation
roll11 <- ugarchroll(gspec11, r_t, solver = "hybrid", keep.coef = TRUE,n.ahead = 1, n.start=1,#
                     refit.every = 5, refit.window = "moving",
                     VaR.alpha = 1-cl )# c(0.025, 0.05))
percent(VaR <- quantile(roll11@forecast$VaR[,1] , probs = 1-cl), accuracy = 0.001)




#install.packages("ufRisk")
library(ufRisk)
tadawul <- as.vector(na.omit(TASI.SR)[1:3751,6])
results = varcast(tadawul, model = 'eGARCH', n.out = 750,   garchOrder = c(3, 3), a.e = cl[3], 
                  a.v = cl[4],
                  distr = c("norm", "std")[1])
cbind(results$VaR.e, results$VaR.v)[250,]*100
results$ES[250]*100




df <- rbind(HS=HS_VaR,gaussian=G_VaR, Cornish_Fisher=Cornish_FisherVaR, MC = MC_VaR)
rownames(df)=NULL
library(xtable)
xtable(apply(df,2,percent, accuracy = 0.001), type = "latex")

df2 <- rbind(HS=HS_ES,gaussian=G_ES, Cornish_Fisher=Cornish_FisherES, MC = MC_ES)
xtable(apply(df2,2,percent, accuracy = 0.001), type = "latex")

xtable(garch_fit@fit$matcoef, digits = 4, type = "latex")

############ Backtest
#install.packages("GAS")
library(GAS)

#BacktestVaR(data, VaR, alpha, Lags = 4)
mat <- as.vector(r_t[3001:3751])
BacktestVaR(mat, HS_VaR[1], 1-cl[1])

#install.packages("segMGarch")
library(segMGarch)
for(p in cl[4]){ print(kupiec(r_t[3001:3751], HS_VaR, 1-p, verbose = TRUE, test = "POF"))}
for(p in cl){ print(kupiec(r_t[3001:3751], Cornish_FisherVaR[4], 1-p, verbose = TRUE, test = "PoF"))}


for(i in 1:4){
  print(paste("alpha = ",1-cl[i]))
print(VaRTest(1-cl[i], as.numeric(r_t[3001:3751]), 
              rep(MC_VaR[,i],751), conf.level = cl[i] ))
}


rt0 <- as.numeric(r_t[3001:3751])
VaR0 <- rep(HS_VaR[4], length(rt0))

VaRLR(as.numeric(r_t[2001:3751]), rep(HS_VaR[1],1751), 0.1, "short")

  
################################################################
# Initialize an empty data frame
Result <- data.frame(Model = character(), AIC = numeric(), BIC = numeric(), SIC=numeric(), Hannan = numeric(), stringsAsFactors = FALSE)
dists <- c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "snig", "QMLE")[1:4]
var.modl <- c("sGARCH","fGARCH","eGARCH","gjrGARCH","apARCH","iGARCH")[-c(2,5)]
p = 1; q = 1
for(v.m in var.modl){
  # Loop over all distributions
  for (dis in dists) {
#    for (p in 1:3) {
#      for (q in 1:3) {
        spec.model <- tryCatch({
          ugarchspec(
            variance.model = list(model = v.m, garchOrder = c(p, q)
                                 # ,submodel="APARCH"
                                  ),
            mean.model = list(armaOrder = c(0, 0),
                              include.mean=TRUE ),
            distribution.model = dis
          )
        }, error = function(e) {
          message("Error in ugarchspec with distribution ", dis, " and GARCH order (", p, ", ", q, "): ", e$message)
          return(NULL)
        })
        
        # If spec.model is NULL, skip to the next iteration
        if (is.null(spec.model)) next
        
        fit <- tryCatch({
          ugarchfit(spec.model, data = r_t, trace = FALSE,solver = "hybrid")
        }, error = function(e) {
          message("Error in ugarchfit with distribution ", dis, " and GARCH order (", p, ", ", q, "): ", e$message)
          return(NULL)
        })
        
        # If fit is NULL, skip to the next iteration
        if (is.null(fit)) next
        # length(infocriteria(fit)[1])==0 | infocriteria(fit)[2]==0
        Result <- rbind(Result, data.frame(Model = paste( dis,"-",v.m, "(", p,",", q,")"), AIC = infocriteria(fit)[1], 
                                           BIC = infocriteria(fit)[2],SIC = infocriteria(fit)[3],Hannan=infocriteria(fit)[4]))
      }
    }
#  }
#}
print(Result)

Result[which.min(Result$AIC),]
Result[which.min(Result$BIC),]

xtable(Result, type = "latex")

plot(Result,col=c(1,2,3,4),pch=20)
#help("accuracy", package = "forecast")


model.arima = auto.arima(r_t[1:3000] , max.order = c(3 , 0 ,3) , stationary = TRUE , trace = T , ic = 'aicc')
model.arima$residuals %>% ggtsdisplay(plot.type = 'hist' , lag.max = 14)
ar.res = model.arima$residuals
Box.test(model.arima$residuals , lag = 14 , fitdf = 2 , type = 'Ljung-Box')
tsdisplay(ar.res^2 , main = 'Squared Residuals')

#################################################################################################################################
mod_spec <- list()
VaR_cal  <-  list()
ES_cal <- list()
# Define model and distribution combinations
models <- c("sGARCH", "eGARCH")[2]
distributions <- c("norm", "snorm", "std", "sstd", "ged", "sged")[-6]

res_names = expand.grid(models, distributions)
q=1
for( mod in models){
  for(dis in distributions){
    mod_spec[[q]] = ugarchspec(variance.model=list(model=mod, garchOrder=c(1,1)), 
                               mean.model=list(armaOrder=c(0,0)),  
                               distribution.model=distributions[q])
    #var.t = ugarchroll(mod_spec[[q]], data = r_t, n.ahead = 1, n.start =2000 , 
    #               refit.every = 5, refit.window = "moving", solver = "solnp",
    #               calculate.VaR = TRUE, VaR.alpha = 1-cl)
    #VaR_cal[[q]] = sapply(seq_along(cl), function(i) quantile(var.t@forecast$VaR[,i],1-cl[i]) )
    #ES_cal[[q]] <- sapply(sapply(seq_along(cl), function(i) var.t@forecast$VaR[, i][var.t@forecast$VaR[, i] <= VaR_cal[[q]][i]]), mean)
    setfixed(mod_spec[[q]]) <- list(mu = 0.01, ma1 = 0.2, ar1 = 0.5, omega = 1e-05,
                            alpha1 = 0.03, beta1 = 0.9, gamma1 = 0.01, delta = 1, shape = 5)
    filt = ugarchfilter(mod_spec[[q]],as.vector(r_t[3001:3751]),solver = "solnp")
    sig = 0.1 
    VaR = fitted(filt) + sigma(filt)*qdist(dis, p=sig, mu = 0, sigma = 1, 
                                           skew  = coef(fit)["skew"], shape=coef(fit)["shape"])

      print(VaRTest(sig, as.numeric(r_t[3001:3751]), 
                    VaR, conf.level = 1-sig ))
  
    q = q+1
  }
}
xtable(apply(matrix(unlist(VaR_cal), ncol = 4,byrow = T),2,percent, accuracy = 0.001), type = "latex")
xtable(apply(matrix(unlist(ES_cal), ncol = 4,byrow = T),2,percent, accuracy = 0.001), type = "latex")


mod_spec <- lapply(distributions, function(dis) ugarchspec(variance.model=
                                                             list(model="eGARCH", garchOrder=c(1,1)), 
                       mean.model=list(armaOrder=c(0,0)),  
                       distribution.model=dis) )

var_function <- function(Spec, data, sig){
  var.t = ugarchroll(Spec, data = data, n.ahead = 1, n.start = 1000, 
                     refit.every = 5, refit.window = "rolling", solver = "solnp",
                     calculate.VaR = TRUE, VaR.alpha = sig)
  report(var.t, type = "VaR", VaR.alpha = sig, conf.level = 1-sig)
}
lapply( 1-cl, function(s) var_function(mod_spec[[5]], r_t, s) )


lapply(seq_along(distributions), function(c) var_function(mod_spec[[c]], r_t, 1-cl[1]) )
