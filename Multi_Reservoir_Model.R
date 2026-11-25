#Lepto_Model_Simulation_LAXMAN2026
rm(list = ls())
# Date: 16/11/2026
# ctrl shift R for sections

## Set working directory and load packages ####
setwd("~/Leptospirosis_Model_Hons")

library(deSolve)

## Model Parameters ####

# Parameter Equation Components 
lambda1 <- 20000
lambda2 <- 500
c1 <- 0.03016
c2 <- 0.00028
h1 <- 10000
h2 <- 2500000
da <- 0.1071
dh <-  0.0012895
dL1 <- 0.15
dL2 <- 0.04
rH <- 0.4046
tau1 <- 0.06
tau2 <- 0.2
mub <- 0.05
alpha <- 3.3940
Tmin <- 7.1 #celcius - antima26ref
Tmax <- 34 # celcius - antima26ref
T <- 27
beta <- 0.3458 #antima26ref
beta2 <- 1.46*10^-6
kappa <- 10000

## Parameter list
parameters <- c(lambda1 = lambda1,
                lambda2 = lambda2,
                c1 = c1,
                c2 = c2,
                beta2 = beta2,
                h1 = h1,
                h2 = h2,
                kappa = kappa,
                da = da,
                dh = dh,
                dL1 = dL1,
                dL2 = dL2,
                rH = rH, 
                tau1 = tau1,
                tau2 = tau2,
                mub = mub,
                alpha = alpha,
                Tmin = Tmin,
                Tmax = Tmax,
                T = T,
                beta = beta)

## Model Variables + Initial Conditions  ####
state <- c(Sh = 50000000,
           Ih = 1,
           Sa = 15000,
           Ia = 1000,
           Bl = 100)

## Model Equations ####
LeptoModelAntimaw_Enviro<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    # force of infection - cant fit in params bcs of time depend. var. 
    En = ((Bl*beta2)/(kappa + Bl))
    Td = alpha*(T - Tmin)*(1 - exp(beta*(T-Tmax)))
    lambdah = c2*Td*(Sh/(Sh+h2))*Ia + En
    lambdaa = c1*Td*(Sa/(Sa+h1))*Ia + En
    
    # rate of change equns
    dSh = lambda2 - lambdah - dh*Sh + rH*Ih
    dIh = lambdah - da*Ih - dL2*Ih - rH*Ih
    dSa = lambda1 - lambdaa - da*Sa 
    dIa = lambdaa - da*Ia - dL1*Ia
    dBl = tau1*Ih + tau2*Ia - mub*Bl
    
    list(c(dSh, dIh, dSa, dIa, dBl)) # rate of change
  }) 
}


## Model Application - TimeScale ####
times <- seq(0, 300, by = 0.01)


## Model Integration ####
out <- ode(y = state, times = times, func = LeptoModelAntimaw_Enviro, parms = parameters)
summary(out)


df <- as.data.frame(out)

## Plotting Time - png ####

png("LeptoModelAntimaw_Enviro.png5", width = 1600, height = 1200, res = 150)

par(oma = c(0, 0, 2, 0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "Sh"], out[, "Bl"], pch = ".")
mtext(outer = TRUE, side = 3, "Leptospirosis Modelled Equations for Transmission", cex = 1.5)

dev.off() ## plots go into file not plots plane 

## Plots for each population ####

# make output into dataframe 
df <- as.data.frame(out)

# making/adding total pops from dataframe 
df$Na <- df$Sa + df$Ia
df$Nh <- df$Sh + df$Ih
df$N <- df$Na + df$Nh


png("LeptoModelAntimaw_Enviro_simulations5.png", width = 1600, height = 1200, res = 160)

# replicate the papers plots via overlapping variables againt time 
LeptoModelAntimaw_Enviro_plot <- function(vars, colrs, lty = 1, main, ylab = "Population", lwd = 2)
{matplot(df$time,df[,vars], type = "l", col = colrs, lty = 1, lwd = 2, xlab = "time", ylab = ylab, main = main)
  legend("right", legend = vars, col = colrs, lty = 1, lwd = 2, cex =0.7, bty = "n")}

## 2x2 grid also png otherwise code breaks 
par(mfrow = c(2,2), oma = c(0, 0, 3, 0), mar = c(4, 4, 3, 1))


## plot 1 -humans 
LeptoModelAntimaw_Enviro_plot(vars = c("Sh", "Ih", "Nh"),
                       colrs = c("steelblue", "navy", "blue"),
                       main = "Human Transmission")

## plot 2 - human exposed, infected, recovered 
LeptoModelAntimaw_Enviro_plot(vars = c("Ih", "Ia"),
                       colrs = c("navy", "darkred"),
                       main = "Infections")

## plot 3 - rodents 
LeptoModelAntimaw_Enviro_plot(vars = c("Sa", "Na", "Ia"),
                       colrs = c("orange", "red", "darkred"),
                       main = "Rodent Transmission")

## plot 4 - populations
LeptoModelAntimaw_Enviro_plot(vars = c("N", "Na", "Nh"),
                       colrs = c("green4", "red", "blue"),
                       main = "Populations")

mtext("LeptoModelAntimaw_Enviro_plot", 
      outer = TRUE, side = 3, cex = 1.4)

dev.off()


