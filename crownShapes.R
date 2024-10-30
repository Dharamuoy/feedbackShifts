library(dplyr)

crownDat <- read.csv("crowns.csv")

w_Mod <- lm(w~Rank+C_H, crownDat)
base::summary(w_Mod)
crownDat$predicted_w <- predict(w_Mod, newdata = crownDat)

Hc_Mod <- lm(Hc~Rank, crownDat[crownDat$Hc != 0,])
base::summary(Hc_Mod)
crownDat$predicted_Hc <- predict(Hc_Mod, newdata = crownDat)

He_Mod <- lm(He~Rank, crownDat[crownDat$Hc != 0,])
base::summary(He_Mod)
crownDat$predicted_He <- predict(He_Mod, newdata = crownDat)

Ht_Mod <- lm(Ht~Rank, crownDat)
base::summary(Ht_Mod)
crownDat$predicted_Ht <- predict(Ht_Mod, newdata = crownDat)
crownDat$mean_Ht <- mean(crownDat$Ht)

# Ash tree height
Mokany <- read.csv("Mokany.csv")
control=nls.control(maxiter=10000, tol=1e-7, minFactor = 1/999999999)
init <- c(a=0.6)
x <- as.numeric(Mokany$X)
y <- as.numeric(Mokany$Y)

Mok<-nls(y~a*x,data=Mokany,start=init,trace=T, control = control)
BSum <- base::summary(Mok)
Ba <- BSum$coefficients[1]
BRSE <- BSum$sigma

#Negative exponential
init1<-c(r=0.5, K = 60)
NE<-nls(y~K * (1-exp(-r*x)),data=Mokany,start=init1,trace=T)
NESum <- base::summary(NE)
rGrowth <- NESum$coefficients[1]
KGrowth <- NESum$coefficients[2]
NERSE <- NESum$sigma
NERsq <- cor(predict(NE, newdata=x), y)**2

windows(3.8,3)
ggplot(data = Mokany,aes(x = as.numeric(X), y = as.numeric(Y))) +
  ggdensity::geom_hdr(probs = c(0.99,0.5,0.25), fill = "darkolivegreen4") +
  geom_point(shape = 21, alpha = 0.5, size = 1) +
  geom_line(aes(x = as.numeric(X), y = predict(NE)), linewidth = 0.75, linetype = 2)+
  labs(x = "Stand age (years)", y = "Tree height (m)") +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=11, colour = "black"),
        axis.text.y  = element_text(size=10, colour = "black"),
        axis.title.y = element_text(size=11, face="bold"),
        axis.title.x = element_text(size=11, face="bold"),
        plot.title = element_text(vjust=1.5, face="bold", size=11, colour = "black"))
