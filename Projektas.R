library(fda.usc)
library(lattice)
library(fda)
library(plotly)
library(lubridate)
library(refund)

library("knitr")
knit2html("~/AG/BANDAU/Projektas.R")

plot_ly(x=1:365,  z = ~return_mat5) %>% add_surface()
plot_ly(x=1:26, y= men,z = returnfd5j$coefs, type = "contour")
LT<-read.delim("~/AG/BANDAU/LTload2.txt", header = T)
attach(LT)
price<-read.delim("~/AG/BANDAU/price.txt", header = T)
attach(price)
load<-as.matrix(load, ncol=1)
plot(load, type="l")
naujas<-matrix(load, nrow=24)

tvec<-seq(1:24)
matplot(tvec, naujas, type = "l")


price1<-as.matrix(price, ncol=1)
kaina<-matrix(price1, nrow=24)
colnames(kaina)<-time
matplot(tvec, kaina, type = "l")


date<-read.delim("~/AG/BANDAU/Date.txt", header = T)
format(as.POSIXct(date,format="%Y/%m/%d %H:%M:%S"),format="%m/%d/%Y")

time<-seq(as.Date("2017-04-30"), as.Date("2018-04-29"), by="days")
laikas<-as.character(time)

savaites<-weekdays(time, abbreviate=TRUE)
menesis<-month(time, label = TRUE)
savaites1<-as.numeric(format(time, format = "%u"))
t<-seq(0,1, len=24)
colnames(naujas)<-laikas

library(lubridate)
men<-month(as.POSIXlt(time, format="%d/%m/%Y"))

#smoothing
#fourier.basis
gcvs <- rep(NA, 24)

for(i in 1:23){
  fbasis <- create.fourier.basis(c(0,1),nbasis=i)
  gcvs[i] <- smooth.basis(t,naujas[,i],fbasis)$gcv
}
t<-seq(0,1, len=24)
t2<-seq(0,1, len=365)
f<-fdata(naujas, t)
a<-min.basis(f, type.CV = GCV.S, numbasis = 24)
a$gcv
matplot(t(a$fdata.est$data)[,2], type = "l", main="GCV with bspline basis")


plot(gcvs,xlab="bases",ylab="gcv", las=1, lwd=1)

i1 <- 15
gcvs[i1]#222.7191

loglam         = seq(-6, 0, 0.25)
Gcvsave        = rep(NA, length(loglam))
names(Gcvsave) = loglam
Dfsave         = Gcvsave
sinRMSE = Gcvsave
for(i in 1:length(loglam)){
  monthbasis<-create.fourier.basis(c(0,1), i1)
  sinefdPar  <- fdPar(monthbasis, Lfdobj=2, 10^loglam[i])
  sine.i    <- smooth.basis(t, naujas, sinefdPar)
  Gcvsave[i] <- sum(sine.i$gcv)
  Dfsave[i]  <-  sine.i$df
  sinemat.i  <- eval.fd(t, sine.i$fd )
  sineres.i  <- naujas - sinemat.i
  sinRMSE[i] <- mean(sqrt(sineres.i^2/(length(rownames(sine.i$fd$coefs))-round(sine.i$df,0))))
  
}

plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(GCV(lambda)), lwd=2 )

plot(loglam, sinRMSE, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(RMSE(lambda)), lwd=2 )

i <- which.min(Gcvsave)
lambda <- 10^loglam[i]
lambda#[1] 5.623413e-06
monthbasis<-create.fourier.basis(c(0,1), 15)
sinefdPar2  <- fdPar(monthbasis, Lfdobj=2, lambda)
sinesmooth2   <- smooth.basis(t, naujas, sinefdPar2)
plotfit.fd(naujas,t,sinesmooth2$fd)
sum(sinesmooth2$gcv)#[1] 172335.4
matplot(sinesmooth2$fd$coefs, type="l")
plot(sinesmooth2$fd)


nw<-eval.fd(t, sinesmooth2$fd)
matplot(naujas-nw, type="l", xlab = "t", ylab="Smoothing error")
skirt<-naujas-nw
max(abs(skirt))
min(abs(skirt))

#bspline basis
t2<-seq(0,1, len=13)
loglam4         <- seq(-6, 0, 0.25)
Gcvsave4        <- rep(NA, length(loglam4))
names(Gcvsave4) <- loglam4
Dfsave4         <- Gcvsave4
sinRMSE4 <- Gcvsave4
for(i in 1:length(loglam4)){
  monthbasis4<-create.bspline.basis(c(0,1), nbasis = 4+length(t)-2, norder = 4, breaks = t)
  sinefdPar4  <- fdPar(monthbasis4, Lfdobj=2, 10^loglam4[i])
  sine.i4    <- smooth.basis(t, naujas, sinefdPar4)
  Gcvsave4[i] <- sum(sine.i4$gcv)
  Dfsave4[i]  <-  sine.i4$df
  sinemat.i4  <- eval.fd(t, sine.i4$fd )
  sineres.i4  <- naujas - sinemat.i4
  sinRMSE4[i] <- mean(sqrt(sineres.i2^2/(length(rownames(sine.i2$fd$coefs))-round(sine.i2$df,0))))
}

plot(loglam4, Gcvsave4, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(GCV(lambda)), lwd=2 )
plot(loglam4, sinRMSE4, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(RMSE(lambda)), lwd=2 )

i4 <- which.min(Gcvsave4)
lambda4i <- 10^loglam4[i4]
lambda4i

t2<-seq(0,1, len=16)
monthbasis4<-create.bspline.basis(c(0,1), nbasis = 4+length(t2)-2, norder = 4, breaks = t2)
sinefdPar4i<- fdPar(monthbasis4, Lfdobj=2, lambda4i)
sinesmooth4i<- smooth.basis(t, naujas, sinefdPar4i)
returnfd4i<-sinesmooth4i$fd
plotfit.fd(naujas,t,returnfd4i)
plot(returnfd4i)
matplot(seq(0,1, length.out = 24), naujas, type = "l", xlab = "Hours", ylab = "Energy consumption")
sum(sinesmooth4i$gcv) #[1] 97628.13

nw2<-eval.fd(t, sinesmooth4i$fd)
matplot(naujas-nw2, type="l", xlab = "t", ylab="Smoothing error")
skirt2<-naujas-nw2
max(abs(skirt2))
min(abs(skirt2))


j <- which.min(sinRMSE4)
lambda2j <- 10^loglam2[j]
lambda2j
sinefdPar2j<- fdPar(monthbasis4, Lfdobj=2, lambda2j)
sinesmooth2j<- smooth.basis(t, naujas, sinefdPar2j)
returnfd2j<-sinesmooth2j$fd
plotfit.fd(naujas,t,returnfd2j)
plot(returnfd2j, col="grey")
matplot(seq(0,1, length.out = 24), naujas, type = "l", xlab = "Hours", ylab = "Energy consumption")
sum(sinesmooth2j$gcv) #[1] 97628.13

returnfd2j$fdnames$day<-savaites1
returnfd2j$fdnames$month<-men

out<-outliers.depth.trim(returnfd2j)#"rep247" "rep304" "rep305" "rep306" "rep270" "rep302"
out$outliers
outliers.depth.pond(returnfd2j)$outliers

matplot(seq(0,1, length.out = 24), naujas, type = "l", xlab = "Hours", ylab = "Energy consumption", col="grey")
lines(returnfd2j$coefs[,247])
lines(returnfd2j[,270])
lines(returnfd2j[,304])
lines(returnfd2j[,305])
lines(returnfd2j[,306])
lines(returnfd2j[,303])
lines(returnfd2j[,215])
lines(returnfd2j[,302])


#remove outlier

a<-naujas[,-306]
b<-a[,-305]
c<-b[,-304]
d<-c[,-303]
e<-d[,-270]
naujas2<-e[,-247]

a<-kaina[,-306]
b<-a[,-305]
c<-b[,-304]
d<-c[,-303]
e<-d[,-270]
kaina2<-e[,-247]


n1<-savaites1[-306]
n2<-n1[-305]
n3<-n2[-304]
n4<-n3[-303]
n5<-n4[-270]
savaites2<-n5[-247]


men<-men[-306]
men<-men[-305]
men<-men[-304]
men<-men[-303]
men<-men[-270]
men<-men[-247]


matplot(t, naujas2, type="l")
matplot(t, kaina2, type="l")
axis(1, at=t, labels = 1:24, padj=2, col="blue", col.axis = "blue")

#bspline basis

loglam5         <- seq(-6, 0, 0.25)
Gcvsave5        <- rep(NA, length(loglam5))
names(Gcvsave5) <- loglam5
Dfsave5         <- Gcvsave5
sinRMSE5 <- Gcvsave5
for(i in 1:length(loglam5)){
  monthbasis5<-create.bspline.basis(c(0,1), nbasis = 4+length(t)-2, norder = 4, breaks = t)
  sinefdPar5  <- fdPar(monthbasis5, Lfdobj=2, 10^loglam5[i])
  sine.i5    <- smooth.basis(t, naujas, sinefdPar5)
  Gcvsave5[i] <- sum(sine.i5$gcv)
  Dfsave5[i]  <-  sine.i5$df
  sinemat.i5  <- eval.fd(t, sine.i5$fd )
  sineres.i5  <- naujas - sinemat.i5
  sinRMSE5[i] <- mean(sqrt(sineres.i5^2/(length(rownames(sine.i5$fd$coefs))-round(sine.i5$df,0))))
}

plot(loglam5, Gcvsave5, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(GCV(lambda)), lwd=2 )
plot(loglam5, sinRMSE5, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(RMSE(lambda)), lwd=2 )

j <- which.min(sinRMSE5)
lambda5j <- 10^loglam5[j]
lambda5j
sinefdPar5j<- fdPar(monthbasis5, Lfdobj=2, lambda5j)
sinesmooth5j<- smooth.basis(t, naujas, sinefdPar5j)
returnfd5j<-sinesmooth5j$fd
plotfit.fd(naujas2,t,returnfd5j)
plot(returnfd5j, xlab="Hours (time)", ylab="Energy Consumption, MW")
matplot(seq(0,1, length.out = 24), naujas2, type = "l", xlab="Hours (time)", ylab="Energy use MW")
sum(sinesmooth5j$gcv) #[1] 94696.99


returnfd5j$fdnames$day<-savaites1
returnfd2j$fdnames$month<-men

outliers.depth.trim(returnfd2j)$outliers#"rep215" "rep302"
outliers.depth.pond(returnfd5j)$outliers#"character(0)"


#vidurkis ir standartinis nuokrypis
vidurkis<-mean(returnfd2j)
plot(returnfd5j, col="grey", lty=1)
lines(vidurkis, lwd=2, col="red")

sfd5 <- std.fd(returnfd5j)
lines(vidurkis+2*sfd5,lwd=2,lty=2, col="blue")
lines(vidurkis-2*sfd5,lwd=2,lty=2, col="blue")

centr<-matrix(NA, ncol=365, nrow = 26)
for (a in 1:26) {
  centr[a,]<-returnfd2j$coefs[a,]-vidurkis$coefs[a]
}
colnames(centr)<-laikas

fd<-returnfd2j
returnfd2j$coefs<-centr
matplot(returnfd2j$coefs, type = "l")

loglam5         <- seq(-6, 0, 0.25)
Gcvsave5        <- rep(NA, length(loglam5))
names(Gcvsave5) <- loglam5
Dfsave5         <- Gcvsave5
sinRMSE5 <- Gcvsave5
for(i in 1:length(loglam5)){
  monthbasis5<-create.bspline.basis(c(0,1), nbasis = 4+length(t)-2, norder = 4, breaks = t)
  sinefdPar5  <- fdPar(monthbasis5, Lfdobj=2, 10^loglam5[i])
  sine.i5    <- smooth.basis(t, centruotas, sinefdPar5)
  Gcvsave5[i] <- sum(sine.i5$gcv)
  Dfsave5[i]  <-  sine.i5$df
  sinemat.i5  <- eval.fd(t, sine.i5$fd )
  sineres.i5  <- naujas - sinemat.i5
  sinRMSE5[i] <- mean(sqrt(sineres.i5^2/(length(rownames(sine.i5$fd$coefs))-round(sine.i5$df,0))))
}

plot(loglam5, Gcvsave5, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(GCV(lambda)), lwd=2 )
plot(loglam5, sinRMSE5, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(RMSE(lambda)), lwd=2 )

j <- which.min(sinRMSE5)
lambda5j <- 10^loglam5[j]
lambda5j
sinefdPar5j<- fdPar(monthbasis5, Lfdobj=2, lambda5j)
sinesmooth5j<- smooth.basis(t, naujas, sinefdPar5j)
returnfd5j<-sinesmooth5j$fd
plotfit.fd(naujas2,t,returnfd5j)
plot(c,returnfd5j, xlab="Hours (time)", ylab="Energy use, MW")
matplot(seq(0,1, length.out = 24), naujas2, type = "l", xlab="Hours (time)", ylab="Energy use MW")
sum(sinesmooth5j$gcv) #[1] 94696.99

#isvestines
plot(deriv.fd(returnfd2j))
plot(deriv.fd(returnfd2j, 2))


#kovariaciju matrica
return.var5<-var.fd(centr)
return_mat5<-eval.bifd(t, t, return.var5)

levelplot(return_mat5,xlab="Hour",ylab="Hour", col.regions = terrain.colors(100))
levelplot((cor.fd(t, returnfd2j)),xlab="Hour",ylab="Hour", col.regions = terrain.colors(100))
persp(t, t, return_mat5, theta = 60, phi = 25, expand = 0.3,ltheta = 120, shade = 0.55,  ticktype = "detailed",xlab = "X", ylab = "Y", zlab = "Z", col = "lightblue")
plot_ly(x=t, y=t, z=return_mat5)

#boxplot/depth
par(mfrow=c(1,2))
fbplot(returnfd2j$coefs, method = "MBD")#101
fbplot(returnfd2j$coefs, method = "BD2")#122
par(mfrow=c(1,1))

out.FM <- depth.FM(returnfd5j, draw=T)
out.RT<-depth.RT(returnfd5j, draw = T)
out.RP<- depth.RP(returnfd5j, draw = T)
print(cur<-c(out.FM$lmed, out.RT$lmed, out.RP$lmed))
plot(returnfd5j,col="grey", xlab = "Time (hours)", ylab = "MW")
lines(returnfd5j[cur], lwd = 2, lty = 1:3, col = 1:3)
legend("topleft", c("FMD", "RPD", "RTD"), lwd = 2, lty = 1:3,col = 1:3, cex = 0.8)

#pca

return.pcalist5 <-pca.fd(returnfd2j, nharm=7)
plot(cumsum(return.pcalist5$varprop), type = "b") #iki 3-ju

#three components
return.pcalist5 <- pca.fd(returnfd2j, nharm=2)
plot.pca.fd(return.pcalist5, expand=.2)

return.rotpcalist5 <- varmx.pca.fd(return.pcalist5)
plot.pca.fd(return.rotpcalist5) 


plot(return.rotpcalist5$scores[,1],return.rotpcalist5$scores[,2], type="p", pch="o",
     xlab="Rotated Harmonic I", ylab="Rotated Harmonic II")
text(return.rotpcalist5$scores[,1]~return.rotpcalist5$scores[,2], labels = colnames(naujas), pos = 4)
returnmat5 <- eval.fd(t, returnfd2j)
returnres5 <-  centr[-(1:2),] - returnmat5

matplot((naujas-(return.pcalist5$scores[,1]+return.pcalist5$scores[,2])), type="l")

#issiaiskinti
returnres.fd5 <- smooth.basis(t, returnres5, sinefdPar2j)$fd
plot(returnres.fd5, lwd=2, col=4, lty=1)
return.pca2 <- pca.fd(returnres.fd5, 5)
plot(return.pca2, expand=0.01)



# grafiskai
returnfd2j$fdnames$day.week<-savaites1
dayw <- ifelse(returnfd2j$fdnames$day.week<5, 4, 3)
matplot(returnfd2j$coefs,col=dayw, type = "l")
legend(1, 1900, legend=c("Working days", "Non-working days"), col = 3:4, cex = 0.7)

m.y <- ifelse(returnfd5j$fdnames$month<4, 28, 8)
matplot(returnfd5j$coefs,col=m.y, type = "l")
par(mfrow=c(1,1))
library(rainbow)
plot(returnfd2j$coefs, plot.type="depth")
data<-fts(x=seq(0,1, len=26),y=returnfd2j$coefs, xname="t", yname="graza")
plot(data, plotlegend=TRUE, plot.type="depth")
plot(data, plotlegend=TRUE, plot.type="function")
legend("topleft", c("FMD", "RPD", "RTD"), lwd = 2, lty = 1:3,col = 1:3, cex = 0.8)
centr<-center.fd(returnfd2j)
par(mfrow=c(1,1))
data5<-fts(x=seq(0,1, len=14),y=returnfd2j$coefs, xname="Hours (time)", yname="Energy consumption")
graf<-fdepth(data, type = c("FM"), trim = 0.25)
plot(data5, plotlegend=TRUE, plot.type="depth")
plot(returnfd2j, plot.type="depth", plotlegend=TRUE, col=rainbow(7), lty=1)
plot.fdepth(graf)

plot.fds(data5, plot.type = c("depth"), colorchoice = c("rainbow"), plotlegend = TRUE)
col1<- ifelse(returnfd5j$fdnames$time<4,3,4)
plot(returnfd5j,col=col1,lwd=1)

par(mfrow=c(1,1))
returnfd2j$fdnames$day.week<-savaites1
dayw <- ifelse(returnfd2j$fdnames$day.week<5, 4, 3)
matplot(returnfd2j$coefs,col=dayw, type = "l")

returnfd2j$fdnames$season<-sezonai

energy.data.frame2<-as.data.frame(returnfd2j$coefs)
x.energy2 <- vector("list", 1) 
x.energy2[[1]] <- as.matrix(energy.data.frame2[, 1:365]) 
group.label.energy2 <- ifelse(returnfd2j$fdnames$day.week<5, 1, 2)

library(fdANOVA)
plotFANOVA(x = x.energy2[[1]], group.label = as.character(group.label.energy2), int = c(0.025, 0.975), col=3:4, separately = TRUE)
fanova.tests(x.energy2[[1]], group.label = group.label.energy2, test = c("L2b"))

energy.data.frame<-as.data.frame(returnfd2j$coefs)
x.energy <- vector("list", 1) 
x.energy[[1]] <- as.matrix(energy.data.frame[, 1:365]) 
group.label.energy <- returnfd2j$fdnames$season

plotFANOVA(x = x.energy[[1]], group.label = as.character(group.label.energy), int = c(0.025, 0.975))
fanova.tests(x.energy[[1]], group.label = group.label.energy)
plotFANOVA(x = x.energy[[1]], group.label = as.character(group.label.energy), int = c(0.025, 0.975), separately = TRUE) 
plotFANOVA(x = x.energy[[1]], group.label = as.character(group.label.energy), int = c(0.025, 0.975), means = TRUE)

returnfd2j$fdnames$month<-men

menesiai <- unique(returnfd2j$fdnames$week) 
p <- length(menesiai) + 1 
menList <- vector("list", p) 
menList[[1]] <- c(rep(1,365),0) 
for (j in 2:p) { 
  xj <- returnfd2j$fdnames$week == menesiai[j-1] 
  menList[[j]] = c(xj,1) 
} 

koef<-returnfd2j$coefs
koef360<-cbind(koef, matrix(0,26,1))
n360fd<-fd(koef360, monthbasis4, returnfd2j$fdnames)
betab<-create.fourier.basis(c(0,1),15)
betafdPar2<-fdPar(betab)
btalist<-vector("list",p)
for (j in 1:p) btalist[[j]] <- betafdPar2
fRegressList <- fRegress(n360fd, menList, btalist) 
btaestList <- fRegressList$betaestlist 
menFit <- fRegressList$yhatfd 
menesiai <- c("1", menesiai) 
par(mfrow=c(2,2),cex=1) 
names(menesiai)<-c("Intercept", "Weekends","Working days", "Working days")
for (j in 1:p) plot(btaestList[[j]]$fd, lwd=2, xlab="Hours", ylab="MW", main=names(menesiai)[j]) 
par(mfrow=c(1,1))
plot(menFit, lwd=2, col=1, lty=1, xlab="Hours", ylab="", main="Prediction")



29+25+30+31+30+31
saltasis<-matrix(NA, ncol = 176, nrow = 26)
saltasisC<-1
30+31+30+31+31+30
siltasis<-matrix(NA, ncol=183, nrow = 26)
siltasisC<-1

for (k in 1:359) {
  if ( returnfd5j$fdnames$month[k] > 3 && returnfd5j$fdnames$month[k] < 10 ){
    siltasis[,siltasisC] = returnfd5j$coefs[,k]
    siltasisC<-siltasisC+1
  } else {
    saltasis[,saltasisC] = returnfd5j$coefs[,k]
    saltasisC<-saltasisC+1
  }
}


ziema<-matrix(NA, ncol=200, nrow = 26)#90
pavasaris<-matrix(NA, ncol = 200, nrow = 26)#92
vasara<-matrix(NA, ncol = 200, nrow = 26)#92
ruduo<-matrix(NA, ncol=200, nrow = 26)#91
ziemaC<-1
pavasarisC<-1
vasaraC<-1
ruduoC<-1

for (k in 1:365) {
  if ( returnfd2j$fdnames$month[k] < 3 && returnfd2j$fdnames$month[k] > 11 ){
    ziema[,ziemaC] = returnfd2j$coefs[,k]
    ziemaC<-ziemaC+1
  } else if ( returnfd2j$fdnames$month[k] > 2 && returnfd2j$fdnames$month[k] < 6 ) {
      pavasaris[,pavasarisC] = returnfd2j$coefs[,k]
      pavasarisC<-pavasarisC+1
  } else if ( returnfd2j$fdnames$month[k] > 5 && returnfd2j$fdnames$month[k] < 9 ) {
        vasara[,vasaraC] = returnfd2j$coefs[,k]
        vasaraC<-vasaraC+1
  } else {
    ruduo[,ruduoC] = returnfd2j$coefs[,k]
    ruduoC<-ruduoC+1
  }
} for (k in 1:364) {
  
}
ifelse(returnfd2j$fdnames$month > 2 && returnfd2j$fdnames$month < 6,1, 0)
matplot(pavasaris, type = "l")
matplot(vasara, type="l")

sezonai<-rep(NA, 365)
for (a in 1:365) {
  if (returnfd2j$fdnames$month[a] == 1) {
    sezonai[a]=1
  }
   else if (returnfd2j$fdnames$month[a] == 2) {
    sezonai[a]=1
  }
  else if (returnfd2j$fdnames$month[a] == 12) {
    sezonai[a]=1
  }
  
  else if (returnfd2j$fdnames$month[a] == 3) {
    sezonai[a]=2
  }
  
  else if (returnfd2j$fdnames$month[a] == 4) {
    sezonai[a]=2
  }
  else if (returnfd2j$fdnames$month[a] == 5) {
    sezonai[a]=2
  }
  else if (returnfd2j$fdnames$month[a] == 6) {
    sezonai[a]=3
  }
  else if (returnfd2j$fdnames$month[a] == 7) {
    sezonai[a]=3
  }
  else if (returnfd2j$fdnames$month[a] == 8) {
    sezonai[a]=3
  }
  else if (returnfd2j$fdnames$month[a] == 9) {
    sezonai[a]=4
  }
  else if (returnfd2j$fdnames$month[a] == 10) {
    sezonai[a]=4
  }
  else if (returnfd2j$fdnames$month[a] == 11) {
    sezonai[a]=4
  }
}
W<-rep("Winter", 90)
Sp<-rep("Spring",92)
S<-rep("Summer", 92)
A<-rep("Autumn",91)
92-32
names(sezonai)<-c(Sp[1:32], S, A, W, Sp[1:60])

returnfd2j$fdnames$season<-sezonai
sum(ifelse(returnfd2j$fdnames$season==1, 1,0))

savd<-rep(NA, 365)
for (a in 1:365) {
  if (returnfd2j$fdnames$day.week[a] < 6) {
    savd[a]=1
  }
  else {
    savd[a]=0
  }
}

returnfd2j$fdnames$week<-savd

yhatmat <- eval.fd(t, menFit$fd) 
rmatb <- vidurkis$coefs - yhatmat
SigmaEb <-var(t(rmatb)) 
y2cMap <- menFit$y2cMap
stderrList <- fRegress.stderr(fRegressList, y2cMap, SigmaEb) 
betastderrlist <- stderrList$betastderrlist 


returnmat5 <- eval.fd(t, returnfd5j)
returnres5 <-  centr[-(1:2),] - returnmat5
returnres.fd5 <- smooth.basis(t, returnres5, sinefdPar5j)$fd
plot(returnres.fd5, lwd=2, col=4)
return.pca2 <- pca.fd(returnres.fd5, 5)
plot(return.pca2, expand=0.01)
constraints <- matrix(c(0,1,1,1,1), 1)

m<-cbind(1, model.matrix(~ factor(returnfd2j$fdnames$) - 1))
olsm<-fosr(fdobj = returnfd2j, X = m, method = "OLS", con=constraints)
plot(olsm,1)

glmm<-fosr(fdobj = returnfd2j, X = m, method = "GLS", con=constraints)
par(mfrow=c(2,3), mar=c(5,2,4,1))
plot(olsm, split=1, set.mfrow=FALSE,titles=c("OLS: Intercept ", levels(factor(returnfd2j$fdnames$season))),ylab="", xlab=" Hours")
plot(glmm, split=1, set.mfrow=FALSE,titles=c("GLS: Intercept ", levels(factor(returnfd2j$fdnames$season))),ylab="", xlab=" Hours")
par(mfrow=c(1,1))
matplot(seq(0,1, len=26), olsm$yhat$coefs, type="l", xlab="Hours", main="Prediction")
apply(returnfd2j$coefs,2, mean)

constraints2 <- matrix(c(0,1,1), 1)
a<-cbind(1, model.matrix(~ factor(returnfd2j$fdnames$week) - 1))
olsma<-fosr(fdobj = returnfd2j, X = a, method = "OLS", con=constraints2)
plot(olsma, split=1, set.mfrow=FALSE,titles=c("OLS: Intercept ", levels(factor(returnfd2j$fdnames$week))),ylab="", xlab=" Hours")
plot(olsma, split=1)
matplot(seq(0,1, len=26), olsma$yhat$coefs, type="l", xlab="Hours", main="Prediction")

library(refund)
constraints2 <- matrix(c(0,1,1,1,1,1,1,1), 1)
a<-cbind(1, model.matrix(~ factor(returnfd5j$fdnames$day.week) - 1))
olsma<-fosr(fdobj = returnfd5j, X = a, method = "OLS", con=constraints2)
plot(olsma, split=1, set.mfrow=FALSE,titles=c("OLS: Intercept ", levels(factor(returnfd5j$fdnames$day.week))),ylab="", xlab=" Hours")
plot(olsma, split=1)
matplot(seq(0,1, len=26), olsma$yhat$coefs, type="l", xlab="Hours", main="Prediction")
