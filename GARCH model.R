rm(list=ls(all=TRUE))

pkges <- c('rugarch', 'xts', 'ggplot2','xtable', 'psych','quantmod', 'PerformanceAnalytics','cowplot','urca','FinTS','forecast' )[-c(5, 8:11)]
sapply(pkges, require, character.only = TRUE) 

tadwaul <- na.omit(getSymbols("^TASI.SR", src = "yahoo", from=as.Date("2007-03-19"),to=as.Date("2024-07-27"),, auto.assign = FALSE))
#df <- tadwaul;df$Date <- time(tadwaul)
r_t <- CalculateReturns(tadwaul[,6], method="log")[-1,]
n <- length(r_t)
testing <-  n*0.2


descrip_Results <- describe( data.frame(real_data = tadwaul$TASI.SR.Adjusted[-1,],returns = r_t))[-c(1,2,5, 6,7,10,13)]

xtable(descrip_Results, type = "latex")

ggplot() +
  geom_histogram(aes(x = r_t),binwidth = 0.01) +
  labs(x = "Returns", y = "Frequency", title = "Histogram of Returns")


p1 <- ggplot() +
  geom_line(aes(x = time(tadwaul), y = tadwaul$TASI.SR.Adjusted), color = "blue")+
   geom_hline(yintercept = mean(tadwaul$TASI.SR.Adjusted) , color = 'red') +
labs(x = "Time", y = "Price", title = "")

p2 <- ggplot() +
  geom_line(aes(x = time(r_t), y = r_t), color = "blue")+
  labs(x = "Time", y = "Return", title = "")


plot_grid(p1, p2, labels = NULL)

############################## Plotting

chart.Histogram(r_t,breaks = "fd", xlab = "Returns", ylab = "Frequency",
                methods = c('add.density', 'add.normal',"add.risk"),
                colorset = c('blue', 'red', 'black'),main = "")
legend("topright", legend = c("return", "kernel", "Normal dist"), 
       col = c("blue", "red","black"), fill= c("blue", "red","black"), lty = 1)
############################### Testing
#Check the Stationary:
#apply ADF test with drift
ADF_Returns = ur.df(r_t, type = "drift",selectlags = "AIC" )
#summary of he test
summary(ADF_Returns)

#Jarque Bera test to check normality
jarque.bera.test(r_t)

#Check for ARCH effect
# use Ljung Box.test from stats package to check auto correlation in squre retruns
Box.test(coredata(r_t^2), type="Ljung-Box", lag = 12)
#ARCH LM Test
ArchTest(r_t)

####################### model

spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                  mean.model = list(armaOrder = c(0,0), include.mean=TRUE),
                  distribution.model = "sstd" ) #the empty function specifies the default model. 
def.fit = ugarchfit(spec = spec, data = r_t, out.sample = 751)
forca <- ugarchforecast(def.fit,n.ahead = 250)
my_fitted <- fitted(forca)
plot(def.fit, which = "all")
head(as.vector(r_t))
head(my_fitted)


forca@forecast$sigmaFor*(250)^(1/2)
# plot of log-returns
plot(cbind(#"fitted"   = def.fit@fit$fitted.values,
           "forecast" = forca@forecast$sigmaFor*(250)^(1/2), #def.fit@fit$sigma,
           "original" = r_t[(n-249):n]), 
     col = c("blue", "red", "black"), lwd = c(0.5, 0.5, 2),
     main = "Forecast of synthetic log-returns", legend.loc = "topleft")

par(mfrow=c(1,2))
plot(forca,which=3)
plot(forca,which=4)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
validation_set <- r_t[c((n-testing):n)]
garch_fore_logreturns <- def.fit@forecast$seriesFor[1, ]
Rt_forecasted <- validation_set[-c(1:(n_val-250))] - garch_fore_logreturns[-1]

plot(validation_set[-c(1:3000)], type = "l", col = "blue", lwd = 2,lty=1, xlab="Time", ylab="Volatility", main = "Forecasting the volatility of MSFT with GARCH(1,2) model vs to validation data",cex.main=0.9)
lines(Rt_forecasted[-1], col = "red", lwd = 2,lty=2)
legend("topright", legend = c("Actual Returns", "Fitted Returns"), col = c("blue", "red"), lty = c(2,3), lwd = 3)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# Example of model specification and estimation

garch_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)), 
                         mean.model = list(armaOrder = c(0,0)),
                         distribution.model = "sstd")

garch_fit <- ugarchfit(spec = garch_spec, data = r_t, out.sample =0 )
head(as.vector(r_t[-1,]))
head(garch_fit@fit$fitted.values)

plot(garch_fit, which = 'all')
plot(garch_fit, which = 3)
par(mfrow=c(1,2))
for(i in c(1,2))
  plot(garch_fit,which=i)


summary(garch_fit)       # ARCH effects are filtered. However, 

rhat <- garch_fit@fit$fitted.values
#rhat2 <- best_model@fit$fitted.values
ggplot() +
  geom_line(aes(x = time(r_t)[1:3000], y = r_t[1:3000]), color = "blue")+
  geom_line(aes(x = time(r_t)[1:3000], y = rhat), color = "red")+
  labs(x = "Time", y = "Returns", title = "MSFT stock prices from 1986-03-14 to 2024-04-10")

hhat <- ts(garch_fit@fit$sigma^2)
plot.ts(hhat)


plot(garch_fit)          # conditional normality seems to be violated

# Forecasting one year ahead
forecast_dates <- Sys.Date() + 1:251
garch_fore <- ugarchforecast(garch_fit, n.ahead = 251, n.roll=50)
par(mfrow=c(1,2))
plot(garch_fore,which=1)
plot(garch_fore,which=2)

par(mfrow=c(1,2))
plot(garch_fore,which=3)
plot(garch_fore,which=4)



# Subsetting forecasted log returns and volatility based on the length of forecast_dates
forecast_log_returns <- xts(garch_fore@forecast$seriesFor[1:length(forecast_dates), 1], order.by = forecast_dates)
forecast_volatility <- xts(garch_fore@forecast$sigmaFor[1:length(forecast_dates), 1], order.by = forecast_dates)


plot(forecast_log_returns,validation_set)
plot(forecast_volatility)


plot(forecast_log_returns, type = "l", col = "blue", lwd = 2, main = "Forecasted vs Actual Log Returns")

# Adding the actual log returns to the plot
lines(validation_set, col = "red", lwd = 2)

# Adding a legend
legend("topright", legend = c("Forecasted Log Returns", "Actual Log Returns"), col = c("blue", "red"), lty = 1, lwd = 2)



# Obtain residuals from the fitted model
residuals <- residuals(garch_fit)

# Plot ACF and PACF of residuals
acf(residuals)
pacf(residuals)

# Check for autocorrelation using Ljung-Box test
Box.test(residuals, lag = 20, type = "Ljung-Box")
# Normality test using Shapiro-Wilk test
shapiro.test(residuals)

# Normal Q-Q plot
qqnorm(residuals)
qqline(residuals)
# Plot squared residuals to visualize volatility clustering
plot(residuals^2, type = "l", main = "Squared Residuals")

# Obtain fitted values from the model
fitted_values <- fitted(garch_fit)

# Plot actual returns against fitted values
plot(r_t, type = "l", col = "blue", lwd = 2, main = "Actual vs Fitted Returns")
lines(fitted_values, col = "red", lwd = 2)
legend("topright", legend = c("Actual Returns", "Fitted Returns"), col = c("blue", "red"), lty = 1, lwd = 2)



garch_fore_logreturns <- xts(ugarchforecast(garch_fit, n.ahead = 1, n.roll = 250)@forecast$seriesFor[1, ],  r_t)


g_spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 3)), 
                         mean.model = list(armaOrder = c(0, 0,1)),distribution.model = "std")

g_fit <- ugarchfit(spec = g_spec, data = r_t)
plot(g_fit)
garchvol = sigma(g_fit)
plot(garchvol, main = "Garch(1,1), alpha = 0.1, beta = 0.8")
sqrt(uncvariance(g_fit))

garchforecast = ugarchforecast(fitORspec = g_fit, 
                               n.ahead = 5)
annualvol = sqrt(252) * sigma(g_fit)
vt_weights = 0.05 / annualvol
plot(merge(annualvol, vt_weights), multi.panel = TRUE, main ='Annualized Portfolio Volatility vs Target Portfolio Weights with 5% annualized volatility')



g_spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 3)), 
                     mean.model = list(armaOrder = c(0, 0,3)),distribution.model = "std")

g_fit <- ugarchfit(spec = g_spec, data = r_t)
plot(g_fit, which = 'all')
garchforecast = ugarchforecast(fitORspec = g_fit, n.ahead = 250)

plot(garchforecast, which = 'all')


############################## model 2


out_of_sample <- round(n_val)
dates_out_of_sample <- tail(index(diff_returns), out_of_sample)

garch_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
                         variance.model = list(model = "sGARCH", garchOrder = c(1,2)),
                          dist='std')
garch_fit <- ugarchfit(spec = garch_spec, data = r_t, out.sample=750)
plot(garch_fit, which = 'all')
coef(garch_fit)
esiduals <- residuals(garch_fit)

# Calculate RMSE
rmse <- sqrt(mean(residuals^2))

dates_out_of_sample <- tail(index(diff_returns), 750)

garch_fore_logreturns <- (ugarchforecast(garch_fit, n.ahead = 1, n.roll = 750)@forecast$seriesFor[1, ])
Rt_forecasted <- r_t[3001:3751] - garch_fore_logreturns

Accual <- r_t[1:3001]
pridected <- as.xts(garch_fit@fit$residuals,order.by = time(r_t[1:3001]))
p <- plot(Accual, type = "l", col = "blue", lwd = 1,lty=1, 
     xlab="Time", ylab="Volatility", 
     main = "",# "Forecasting the volatility of MSFT with EGARCH(1,1) model vs to validation data",
     cex.main=0.9)
lines(pridected, col = "red", lwd = 1,lty=2)
legend("topright", legend = c("Actual Returns", "Fitted Returns"), 
       col = c("blue", "red"), lty = c(1,2), lwd = 1)
ggsave(filename ="ActualvsFitted.png", plot = p, width = 4, height = 4, units = "in", dpi = 300)


actual_values <- validation_set
forecast_values <- as.numeric(garch_fore_logreturns)

# Compute errors
errors <- actual_values - forecast_values[-1]

# Compute RMSE
rmse <- sqrt(mean(errors^2))

# Compute MAPE
mape <- mean(abs(errors/actual_values)) * 100

# Print results
print(paste("RMSE:", rmse))
print(paste("MAPE:", mape))

# Set a small constant value
epsilon <- 1e-8

# Compute MAPE
mape <- mean(abs(errors / pmax(abs(actual_values), epsilon))) * 100


library(rugarch)

# Example data
set.seed(123)
your_data <- rnorm(1000)


# Extract and standardize AIC and BIC values
info_criteria <- infocriteria(garch_fit)
AIC <- info_criteria[1]
BIC <- info_criteria[2]
n_obs <- n
AIC_standard <- AIC * n_obs
BIC_standard <- BIC * n_obs

# Output the standardized values
print(paste("Standardized AIC:", AIC_standard))
print(paste("Standardized BIC:", BIC_standard))

par(mfrow=c(1,2))
for(i in 10:11)
plot(garch_fit,which=i)


df <- data.frame(
  Time = time(r_t[3000:3751]),
  Actual_Returns = as.vector(r_t[3000:3751]),
  Fitted_Returns = as.vector(garch_fit@fit$residuals[3000:3751])
)

df_melted <- reshape2::melt(df, id = "Time", variable.name = "Type", value.name = "Returns")
df_melted$Type <- as.factor(df_melted$Type)
# Plot using ggplot2
ggplot(df_melted, aes(x = Time, y = Returns, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Actual_Returns" = "blue", "Fitted_Returns" = "red")
                     #,labels = c("Actual Returns", "Fitted Returns")
                     ) +
  #scale_linetype_manual(values = c("Actual_Returns" = "dashed", "Fitted_Returns" = "dashed"),
  #                      labels = c("Actual Returns", "Fitted Returns")) +
  labs(
    x = "Time", 
    y = "Volatility", 
   # title = "Forecasted the Volatility of the Tadawul Exchange Using a GARCH(1,1) Model with a Skewed Student-t Distribution Compared to Validation Data",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10)
  )




##################################
# Install and load the necessary packages

install.packages("patchwork")
library(ggplot2)
library(patchwork)

# Assuming r_date, MSFT_NASDAQ, and r_t are already defined and n is specified
# Example data (replace with your actual data)
# r_date <- seq(as.Date("2020-01-01"), by = "day", length.out = 100)
# MSFT_NASDAQ <- rnorm(100, mean = 200, sd = 10)
# r_t <- rnorm(100, mean = 0, sd = 1)
# n <- 100

# Create the first plot
p1 <- ggplot() +
  geom_line(aes(x = r_date[1:n], y = MSFT_NASDAQ[1:n]), color = "blue") +
  geom_hline(yintercept = mean(MSFT_NASDAQ[1:n]), color = 'red') +
  labs(x = "Time", y = "Price", title = "")

# Create the second plot
p2 <- ggplot() +
  geom_line(aes(x = r_date[1:n], y = r_t[1:n]), color = "blue") +
  labs(x = "Time", y = "Return", title = "")

# Combine the plots using patchwork
combined_plot <- p1 | p2

# Display the combined plot
print(combined_plot)


# Calculate VaR using PerformanceAnalytics
library(PerformanceAnalytics)
VaR(r_t, p = 0.95, method = "historical")
VaR(r_t, p = 0.95, method = "gaussian")
VaR(r_t, p = 0.95, method = "modified")
VaR(r_t, p = 0.95, method = "kernel")

chart.VaRSensitivity(r_t,
                     methods=c("HistoricalVaR", "ModifiedVaR", "GaussianVaR"),
                     colorset=bluefocus, lwd=2)

sapply(paste("package", pkges, sep = ":"), function(pk) detach(pk, unload = TRUE, character.only = TRUE, force = TRUE))
