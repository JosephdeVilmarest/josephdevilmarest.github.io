data <- readRDS('data.RDS')
head(data)

data <- data[which(data$Instant == 36),]
head(data)

##### Analyse graphique #####
sel <- which(data$Year < 2020)
plot(data$Toy[sel], data$Consommation[sel], pch=20, cex=0.3)
plot(data$Temperature[sel], data$Consommation[sel], pch=20, cex=0.3)
plot(data$TempLissLent[sel], data$Consommation[sel], pch=20, cex=0.3)
plot(data$TempLissRapide[sel], data$Consommation[sel], pch=20, cex=0.3)
plot(data$VitesseVent[sel], data$Consommation[sel], pch=20, cex=0.3)
plot(data$CPWind[sel], data$Consommation[sel], pch=20, cex=0.3)
plot(data$Nebulosite[sel], data$Consommation[sel], pch=20, cex=0.3)
plot(data$Lag1J[sel], data$Consommation[sel], pch=20, cex=0.3)


##### Partie 1

##### 1. Creation d'un GAM #####
train <- which(data$Year < 2019)
test <- which(data$Year == 2019)
# Avec termes AR
eq <- Consommation ~ s(Toy, k=20, bs='cc') + as.factor(DayType) + s(TempLissLent, k=5) +
  s(TempLissRapide, k=5) + s(CPWind) + Lag1J:as.factor(DayType) + Lag1S + te(Tend, Temperature, k=c(3,5)) +
  s(Tend)
# Sans termes AR
# eq <- Consommation ~ s(Toy, k=20, bs='cc') + as.factor(DayType) + s(TempLissLent, k=5) +
#   s(TempLissRapide, k=5) + s(CPWind) + te(Tend, Temperature, k=c(3,5)) + s(Tend)

g <- mgcv::gam(eq, data=data[train,])
yhat_gam <- predict(g, newdata=data)

plot(data$Date[test], data$Consommation[test], type='l')
lines(data$Date[test], yhat_gam[test], col='blue')

for (year in 2019:2022)
  print(c(year, sqrt(mean((data$Consommation-yhat_gam)[which(data$Year == year)]^2))))


##### 2. Definition de la matrice X #####
X <- predict(g, newdata=data, type='terms')
for (j in 1:ncol(X))
  X[,j] <- X[,j] / sd(X[train,j])
X <- cbind(X,1)


##### 3. Optimisation des hyper-parametres #####
sel_lent <- which(data$Year < 2019)
sel_rapide <- which(data$Year == 2020)
sel_moyen <- which(data$Year <= 2020)
params <- viking::iterative_grid_search(X[sel_lent,], data$Consommation[sel_lent], 
                                        2^(-30:0), p1=1, ncores = 4)

ssm <- viking::statespace(X, data$Consommation, kalman_params = params)
plot(ssm, date = data$Date, sel = which(data$Year > 2018), window_size = 28)

yhat_igs <- ssm$pred_mean
for (year in 2019:2022)
  print(c(year, sqrt(mean((data$Consommation-yhat_igs)[which(data$Year == year)]^2))))


##### 4. Introduction d'un delai #####
# Attention : si on a un delai on ne peut plus utiliser de termes AR dans le modele !
del <- 7
th_arr <- ssm$kf$theta_arr[sapply(1:nrow(X), function(t) {max(1,t-del+1)}),]
yhat_igs <- rowSums(X * th_arr)
for (year in 2019:2022) {
  print(c(year, sqrt(mean((data$Consommation-yhat_igs)[which(data$Year == year)]^2))))
}

params <- viking::iterative_grid_search(X[sel_lent,], data$Consommation[sel_lent], 
                                        2^(-30:0), p1=1, ncores = 4, delay = del)
ssm <- viking::statespace(X, data$Consommation, kalman_params = params)
th_arr <- ssm$kf$theta_arr[sapply(1:nrow(X), function(t) {max(1,t-del+1)}),]
yhat_igs <- rowSums(X * th_arr)
for (year in 2019:2022)
  print(c(year, sqrt(mean((data$Consommation-yhat_igs)[which(data$Year == year)]^2))))






##### Partie 2

sel_precovid <- which(data$Date < lubridate::as_date('2020-03-16'))
ssm0 <- viking::statespace(X[sel_precovid,], data$Consommation[sel_precovid],
                          kalman_params = params)
theta_final <- ssm0$kf$theta
P_final <- ssm0$kf$P

sel_covid <- which(data$Date >= lubridate::as_date('2020-03-16'))
Q_big <- diag(ncol(X), x=params$sig^2)
ssm1 <- viking::statespace(X[sel_covid,], data$Consommation[sel_covid],
                           kalman_params = list(theta=theta_final,
                                                P=P_final-params$Q+Q_big,
                                                Q=params$Q, sig=params$sig))

theta_arr <- rbind(ssm0$kf$theta_arr, ssm1$kf$theta_arr)
yhat_rupture <- rowSums(theta_arr * X)
for (year in 2019:2022)
  print(c(year, sqrt(mean((data$Consommation-yhat_rupture)[which(data$Year == year)]^2))))

sel <- which(data$Year == 2020)
w <- 28
plot(data$Date[sel], sapply(sel, function(t) {mean((data$Consommation-yhat_igs)[(t-w+1):t])}), type='l', ylim=c(-4000,2000), lwd=2)
lines(data$Date[sel], sapply(sel, function(t) {mean((data$Consommation-yhat_rupture)[(t-w+1):t])}), col='darkgreen', lwd=2)

plot(data$Date[sel], sapply(sel, function(t) {sqrt(mean((data$Consommation-yhat_igs)[(t-w+1):t]^2))}), type='l', lwd=2)
lines(data$Date[sel], sapply(sel, function(t) {sqrt(mean((data$Consommation-yhat_rupture)[(t-w+1):t]^2))}), col='darkgreen', lwd=2)


