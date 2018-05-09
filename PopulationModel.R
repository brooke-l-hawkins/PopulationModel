
library(deSolve)

#### RUNTIME ###################################################################

# start timing script
start <- proc.time()

#### PARAMETERS ################################################################

# t subscript in variable name indicates temperature-sensitive
# base terms are now the scaling factor for temperature-sensitive terms

# z: adult to juvenile size ratio
z<-0.2
# p: adult modifier on production of juveniles
# can use to modify stage structure, seems like it needs to be <1
p<-0.4
# sig: conversion efficiency of resource to useful energy for fish
sig<-0.7
# M: maximum ingestion rate
# make mass specific?
M<-0.5
# MS: function breadth for max intake rate
# i.e. attack rate, opens downward
MS<-10
# H: handling time
# speed at which fish can eat resources (smaller is faster)
H<-1
# HS: function breadth for handling time
# temperature parabola opens upward
HS<-10
# uJ: juvenile mortality
# uJ = 0.05 when C = 20
uJ<-1.625
# uJe: activation energy of juvenile mortality
uJe<--0.006
# uA: adult mortality
# uA = 0.005 when C = 20
uA<-1.625
# uAe: activation energy og adult mortality
uAe<--0.006
# uR: resource mortality
# uR = 0.005 when C = 20
uR<-0.163
# uRe: activation energy of resource mortality
uRe<--0.006
# t: costs of maintaining somatic growth/turnover
# i.e. base level of resource intake you must exceed to mature/reproduce
# t = 0.1 when C = 20
t<-0.326
# te: activation energy of metabolic waste
te<--0.006
# r: resource growth rate
r<-1.5
# rS: function breadth for resource growth rate
# temperature parabola opens downward
rS<-15
# K: resource carrying capacity
K<-5
# B: adult reproductive rate
B<-0.5
# boltzmann's constant
kb<-8.617*10^-5
# t: temperature in degrees celsius
# 20 is optimal, 10 is cold, 30 is hot
C<-20

#### ODE FUNCTIONS #############################################################

# t subscript in variable name indicates temperature-sensitive

# Mt: maximum ingestion rate
# Ht: handling time
# tt: costs of maintaining somatic growth/turnover
# uJt: mortality rate for juveniles
# uAt: mortality rate for adults
# uRt: mortality rate for resources
# rt: resource growth rate
# ca: functional response for adults
# cj: functional response for juveniles
# mj: juvenile maturation rate
# ra: reproduction per adult

BaseStaget<-function(t,y,p){
    {
        J<-y[1]
        A<-y[2]
        R<-y[3]
    }
    with(as.list(p),{
        Mt<-M*exp(-(C-20)^2/(2*MS)^2)
        Ht<-H*exp((C-20)^2/(2*HS)^2)
        tt<-t*exp(te/(kb*C))
        uJt<-uJ*exp(uJe/(kb*C))
        uAt<-uA*exp(uAe/(kb*C))
        uRt<-uR*exp(uRe/(kb*C))
        rt<-r*exp(-(C-23)^2/(2*rS)^2)
        
        ca<-Mt*(R/(Ht+R))
        cj<-Mt*(R/(Ht+R))
        ifelse(((sig*cj)-tt)<0, mj<-0, mj<-(((sig*cj)-tt)-uJt)/(1-z^(1-(uJt/((sig*cj)-tt)))))
        ifelse(((sig*ca)-tt)<0, ra<-0, ra<-((sig*ca)-tt)*B)
        
        dJ.dt<- ra*A -mj*J - uJt *J
        
        dA.dt<- mj*J - uAt*A -ra*p*A
        
        dR.dt<- rt*R*(1-(R/K)) - cj*J - ca*A - uRt*R
        
        return(list(c(dJ.dt,dA.dt,dR.dt)))
    })
}

#### SIMULATION ################################################################

# Set Plot Preferences ---------------------------------------------------------

# plot juvenile, adult, and resource dynamics?
dynamics.plot = T
# plot bifurcation plot?
bifrucation.plot = T
# layout for plots
par(mfcol=c(4,5), mar=c(1,1,1,1))

iterations <- 0:10000
iterations <- iterations/10
if (bifrucation.plot) {
    # which state variable should be plotted?
    bifurcation.plot.variable <- "Adults"
    # which proportion of indices should be plotted?
    # ex. 0.5 will plot last half of iterations
    bifurcation.plot.portion <- 0.5
}

# Set State Variables ----------------------------------------------------------

juv.vec <- c(1,2,3,4,5) # juveniles
a <- 1 # adults
r <- 1 # resources

# Set Parameters ---------------------------------------------------------------

# choose parameter to change
parm.name <- "C"
# choose range of parameters
parm.seq <- seq(10,30,length = 3)

# ------------------------------------------------------------------------------

# create parameter vector
parms<-c(z=z, p=p, sig=sig, M=M, MS=MS, H=H,HS=HS, uJ=uJ, uJe=uJe,
         uA=uA, uAe=uAe, uR=uR, uRe=uRe, t=t, te=te,r=r, rS=rS, K=K, B=B, kb=kb, C=C)
# find index of parameter to change
parm.index <- which(names(parms)==parm.name)

# initialize list to store output
output <- list()

# loop to change initial value of juveniles
for (j in juv.vec) {
    y <- c(j,a,r)
    names(y) <- c("Juveniles", "Adults", "Resources")
    
    
    # loop to change parameter
    for (i in 1:length(parm.seq)) {
        # set parameters
        parms.loop <- parms
        parms.loop[parm.index] <- parm.seq[i]

        # run simulation
        output[[i]] <- ode(y=y, time=iterations, BaseStaget, parms=parms.loop)
        # remove time column
        output[[i]] <- output[[i]][,-1]
        
        # juvenile, adult, and resource dynamics plot
        if (dynamics.plot) {
            matplot(output[[i]],type="l",lty=1,pch=0.5,col=1:3,
                    xlab=paste0("J=",j," ","A=",a," ","R=",r," ", parm.name,"=", parm.seq[i]))
        }

    } # end parm.seq loop

    # bifurcation plot
    if (bifrucation.plot) {
        
        # choose rows to plot
        total.rows <- length(iterations)
        bifurcation.first.row <- round(bifurcation.plot.portion * total.rows)
        bifurcation.last.row <- total.rows
        bifurcation.rows <- bifurcation.first.row:bifurcation.last.row
        
        # set up plot axes and labels
        # TODO change range.lim to depend on bifurcation.plot.portion
        range.lim <- lapply(output, function(x) apply(x, 2, range))
        range.lim <- apply(do.call("rbind", range.lim), 2, range)
        plot(0, 0, pch = "", xlab=paste0("J=",j," ","A=",a," ","R=",r," "),
             ylab = bifurcation.plot.variable, xlim = range(parm.seq),
             ylim = range.lim[,bifurcation.plot.variable])
        
        # plot points
        for (i in 1:length(parm.seq)) {
            points(rep(parm.seq[i], length(bifurcation.rows)), output[[i]][bifurcation.rows,bifurcation.plot.variable])
        }
    }
    
} # end juv.vec loop

#### RUNTIME ###################################################################

# stop timing script
end <- proc.time()
# print elapsed time
print(end[3]-start[3])

