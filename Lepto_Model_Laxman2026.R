#Lepto_Model_Simulation_LAXMAN2026

# Date: 16/11/2026
# ctrl shift R for sections

## Set working directory and load packages ####
setwd("~/Leptospirosis_Model_Hons")

library(deSolve)

## Model Parameters ####

# Parameter Equation Components 
Lambdah <- 37013.84
Lambdaa <- 140.8
delta <- 0.3
r <- 0.5351
rho <- 0.3
muh <- 0.001108
mua <- 0.00704
#beta1 <- 3.1*10^-7 # bifurcation parameter
zeta <- 1.3
alpha <- 0.2
gamma <- 0.623
sigma <- 0.1312
v <- 0.0409
d1 <- 0.0004
d2 <- 0.1434
xi <- 0.5039
beta2 <- 1.46*10^-6
beta3 <- 3.85*10^-6
  
## Parameter list
parameters <- c(Lambdah = Lambdah,
                Lambdaa = Lambdaa,
                delta = delta,
                #beta1 = beta1,
                beta2 = beta2,
                beta3 = beta3,
                gamma = gamma,
                muh = muh,
                mua = mua,
                xi = xi,
                alpha = alpha, 
                v = v,
                rho = rho,
                sigma = sigma,
                zeta = zeta,
                r = r,
                d1 = d1,
                d2 = d2,
                k1 = muh + gamma,
                k2 = d1 + sigma + muh + v,
                k3 = muh + xi,
                k4 = muh + r,
                k5 = d2 + mua)

## Model Variables + Initial Conditions  ####
state <- c(Shh = 5000000,
           Slh = 28400000,
           Eh = 3000,
           Ah = 235,
           Ih = 10,
           Rh = 0,
           Sa = 19000,
           Ia = 1000)

## Model Equations ####
LeptoModLaxman<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    # force of infection - cant fit in params bcs of time depend. var. 
    lambda1 = beta1*Ia
    lambda2 = beta3*Ia + beta2*Ih
    beta1 = (1 - Rv)/Ra
    
    # rate of change equns
    dShh = delta*Lambdah + rho*r*Rh - muh*Shh - lambda1*zeta*Shh
    dSlh = (1-delta)*Lambdah + (1-rho)*r*Rh - muh*Slh - lambda1*Slh
    dEh = (zeta*Shh + Slh)*lambda1 - k1*Eh
    dAh = (1-alpha)*gamma*Eh - k2*Ah
    dIh = sigma*Ah + alpha*gamma*Eh - k3*Ih
    dRh = v*Ah + xi*Ih - k4*Rh
    dSa = Lambdaa - lambda2*Sa - mua*Sa
    dIa = lambda2*Sa - k5*Ia
    
    
    list(c(dShh, dSlh, dEh, dAh, dIh, dRh, dSa, dIa)) # rate of change
  }) 
}


## Model Application - TimeScale ####
times <- seq(0, 300, by = 0.01)


## Model Integration ####
out <- ode(y = state, times = times, func = LeptoModLaxman, parms = parameters)
summary(out)


df <- as.data.frame(out)

## Plotting Time - png ####

png("lepto_laxman.png", width = 1600, height = 1200, res = 150)

par(oma = c(0, 0, 2, 0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "Shh"], out[, "Rh"], pch = ".")
mtext(outer = TRUE, side = 3, "Leptospirosis Modelled Equations for Transmission", cex = 1.5)

dev.off() ## plots go into file not plots plane 

## Plots for each population ####

# make output into dataframe 
df <- as.data.frame(out)

# making/adding total pops from dataframe 
df$Na <- df$Sa + df$Ia
df$Nh <- df$Shh + df$Eh + df$Ih + df$Rh + df$Slh + df$Ah
df$N <- df$Na + df$Nh


png("lepto_simulation_laxman.png", width = 1600, height = 1200, res = 160)

# replicate the papers plots via overlapping variables againt time 
paper_replicant_laxman <- function(vars, colrs, lty = 1, main, ylab = "Population", lwd = 2)
{matplot(df$time,df[,vars], type = "l", col = colrs, lty = 1, lwd = 2, xlab = "time", ylab = ylab, main = main)
  legend("right", legend = vars, col = colrs, lty = 1, lwd = 2, cex =0.7, bty = "n")}

## 2x2 grid also png otherwise code breaks 
par(mfrow = c(2,2), oma = c(0, 0, 3, 0), mar = c(4, 4, 3, 1))


## plot 1 -humans 
paper_replicant_laxman(vars = c("Shh", "Slh", "Eh", "Ah", "Ih", "Rh"),
                     colrs = c("steelblue","steelblue", "deepskyblue","blue", "blue", "navy"),
                     main = "Human Transmission")

## plot 2 - human exposed, infected, recovered 
paper_replicant_laxman(vars = c("Eh", "Ah", "Ih", "Rh"),
                     colrs = c("deepskyblue", "blue", "navy", "black"),
                     main = "Human Exposed, Acute, Chronic, Recovered")

## plot 3 - rodents 
paper_replicant_laxman(vars = c("Sa", "Ia", "Na"),
                     colrs = c("orange", "red", "darkred"),
                     main = "Rodent Transmission")

## plot 4 - populations
paper_replicant_laxman(vars = c("N", "Na", "Nh"),
                     colrs = c("green4", "darkred", "navy"),
                     main = "Populations")

mtext("Leptospirosis Model Plot - Laxman", 
      outer = TRUE, side = 3, cex = 1.4)

dev.off()


