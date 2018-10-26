
library(deSolve)

#### RUNTIME ###################################################################

# start timing script
start <- proc.time()

#### VARIABLES & PARAMETERS ####################################################

# Create parameter vector ------------------------------------------------------

# Temperature
t0   <- 292.15 # Reference temperature (K)
k    <- 8.617*10^-5 # Boltzmann's constant (eV * K^-1)

# Mortality
phi1 <- 0.0015 # Allometric scalar of background mortality (g^(1-phi2) * day^-1)
phi2 <- 0.75   # Allometric exponent of background mortality (no units)
E.mortality  <- 0.45   # Activation energy of background mortality (eV)

# Metabolism
rho1  <- 0.0123 # Allometric scalar of metabolism (g^(1-rho2) * day^-1)
rho2  <- 0.77 # Allometric exponent of metabolism at reference temperature (no units)
E.metabolism  <- 0.594 # Activation energy of metabolism (eV)

# Consumption
max.attack <- 3 # Maximum zooplankton (1mm) attack rate (100 m^3 * day^-1)
mass.opt <- 41 # Optimal forager size for 1mm zooplankton (g)
alpha <- 0.75 # Allometric exponent of attack rate at reference temperature (no units)
eps1 <- 0.25 # Allometric scalar of maximum ingestion rate (g^(1-eps2) * day^-1)
eps2 <- 0.77 # Allometric exponent of maximum ingestion rate at reference temperature (no units)
E.consume <- 1.21 # Activation energy of functional response parameters (eV)
sigma <- 0.5 # Assimilation efficiency (no units)

# Ontogeny and Maturation
length.j <- 12 # Length at onset of active feeding (mm)
length.a <- 140 # Length at onset of active feeding (mm)
lam1 <- 0.00794 # Constant in length-weight relationship (mm * g^(-lam2))
lam2 <- 3.15 # Exponent in length-weight relationship (no units)
size.j <- lam1 * length.j^lam2 # Size of newborn (g)
size.a <- lam1 * length.a^lam2 # Size of adult (g)
z <- size.j/size.a # Newborn to adult size ratio (no units)
rm(length.j, length.a, lam1, lam2, size.j, size.a)

# Resource
delta <- 0.1 # Turnover rate of shared and adult resource (day^-1)
R.max <- 1 # Maximum resource biomass density (g * (100 m^3)^-1)

parameters <- c(t0 = t0, k = k, phi1 = phi1, phi2 = phi2, 
                E.mortality = E.mortality, rho1 = rho1, rho2 = rho2, 
                E.metabolism = E.metabolism, max.attack = max.attack, 
                mass.opt = mass.opt, alpha = alpha, eps1 = eps1, eps2 = eps2, 
                E.consume = E.consume, sigma = sigma, z = z, delta = delta, 
                R.max = R.max)

# Create state variable vector -------------------------------------------------

j.initial <- 1 # juveniles
a.initial <- 1 # adults
r.initial <- 1 # resource

state.variables <- c(J = j.initial, A = a.initial, R = r.initial)


#### CHANGING PARAMETERS #######################################################

# Temeprature (K)
t.sequence <- seq(from = t0-10, to = t0+10, length = 3)
# Temperature effect on mass-scaling exponent of metabolism (K^-1)
cM.sequence <- seq(from = -0.2, to = 0.2, by = 0.1)

# create parameter matrix, one row represents parameters for one simulation
simulations <- length(t.sequence) * length(cM.sequence)
t <- rep(x = t.sequence, each = length(cM.sequence))
cM <- rep(x = cM.sequence, times = length(t.sequence))
rm(t.sequence, cM.sequence)

parameters.matrix <- matrix(data = parameters,
                            nrow = simulations,
                            ncol = length(parameters))
colnames(parameters.matrix) <- names(parameters)
parameters.matrix <- cbind(t, parameters.matrix, cM)


#### ODE FUNCTIONS #############################################################

Model <- function(t, y, p){
    with(as.list(c(y, p)),{
        
        # Mortality
        mortality.J <- exp((E.mortality*(t-t0))/(k*t*t0)) * phi1 * J^phi2
        mortality.A <- exp((E.mortality*(t-t0))/(k*t*t0)) * phi1 * A^phi2
        
        # Metabolism
        metabolism.J <- exp((E.metabolism*(t-t0))/(k*t*t0)) * rho1 * J^(rho2 + cM*(t-t0))
        metabolism.A <- exp((E.metabolism*(t-t0))/(k*t*t0)) * rho1 * A^(rho2 + cM*(t-t0))
        
        # Consumption
        attack.J <- exp((E.consume*(t-t0))/(k*t*t0)) * max.attack * ((J/mass.opt)*exp(1-(J/mass.opt)))^alpha
        attack.A <- exp((E.consume*(t-t0))/(k*t*t0)) * max.attack * ((A/mass.opt)*exp(1-(A/mass.opt)))^alpha

        encounter.J <- attack.J * R
        encounter.A <- attack.A * R
        
        max.ingestion.J <- exp((E.consume*(t-t0))/(k*t*t0)) * eps1 * J^(eps2)
        max.ingestion.A <- exp((E.consume*(t-t0))/(k*t*t0)) * eps1 * A^(eps2)

        ingestion.J <- encounter.J / (1 + encounter.J/max.ingestion.J)
        ingestion.A <- encounter.A / (1 + encounter.A/max.ingestion.A)

        net.biomass.J <- sigma * ingestion.J - metabolism.J
        net.biomass.A <- sigma * ingestion.A - metabolism.A

        if (net.biomass.J > 0)  {positive.biomass.J <- net.biomass.J}
        else                    {positive.biomass.J <- 0}
        if (net.biomass.A > 0)  {positive.biomass.A <- net.biomass.A}
        else                    {positive.biomass.A <- 0}
        
        # Maturation
        maturation.J <- (positive.biomass.J - mortality.J) / (1 - z^(1 - (mortality.J/positive.biomass.J)))
        
        # Rate of change
        dJ.dt <- positive.biomass.A*A + net.biomass.J*J - maturation.J*J - mortality.J*J
        dA.dt <- maturation.J*J + net.biomass.A*A - positive.biomass.A*A - mortality.A*A
        dR.dt <- delta*(R.max-R) - ingestion.J*J - ingestion.A*A
        
        # Return rate of change
        return(list(c(dJ.dt, dA.dt, dR.dt)))
    })
}

#### SIMULATION ################################################################


# Choose How Long to Run -------------------------------------------------------
time.steps <- 0:100
time.steps <- time.steps/100

output <- list()
for (ss in 1:simulations) {
    # print progres statement
    print(paste0("Starting simulation for simulation ", ss))
    # run model
    output[[ss]] <- ode(y = state.variables, 
                        times = time.steps, 
                        func = Model,
                        parms = parameters.matrix[ss, ])
    # remove time column
    output[[ss]] <- output[[ss]][,-1]
    # plot output
    matplot(output[[ss]], type = "l", lty = 1, pch = 0.5, col = 1:3,
            xlab = '', ylab = '')
    title(main = paste0("Temp = ", parameters.matrix[ss, "t"], 
                        "  |  cM = ", parameters.matrix[ss, "cM"]))
}


#### RUNTIME ###################################################################

# stop timing script
end <- proc.time()
# print elapsed time
print(end[3]-start[3])


