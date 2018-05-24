
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
# uJt = 0.05 when C = 293.15
uJ<-10000
# uJe: activation energy of juvenile mortality
uJe<--0.308
# uA: adult mortality
# uAt = 0.005 when C = 293.15
uA<-10000
# uAe: activation energy of adult mortality
uAe<--0.308
# uR: resource mortality
# uR = 0.005 when C = 293.15
uR<-0.163
# uRe: activation energy of resource mortality
uRe<--0.006
# t: costs of maintaining somatic growth/turnover
# i.e. base level of resource intake you must exceed to mature/reproduce
# t = 0.1 when C = 293.15
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
# C: temperature in kelvin
# 283.15 is cold, 293.15 is optimal, 303.15 is hot
C<-293.15
# Copt: optimal temperature in kelvin
Copt<-293.15

#### ODE FUNCTIONS #############################################################

# t subscript in variable name indicates temperature-sensitive

# Mt: maximum ingestion rate
# Ht: handling time
# tt: costs of maintaining somatic growth/turnover
# uJt: mortality rate for juveniles
# uAt: mortality rate for adults
# uRt: mortality rate for resources
# rt: resource growth rate
# qt: competitive difference between adults and juveniles
    # qt=1 when juveniles and adults are equal competitors
    # qt=0.5 when juveniles are 3x better than adults
    # qt=1.5 when adults 3x better than juveniles
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
        Mt<-M*exp(-(C-Copt)^2/(2*MS)^2)
        Ht<-H*exp((C-Copt)^2/(2*HS)^2)
        tt<-t*exp(te/(kb*C))
        uJt<-uJ*exp(uJe/(kb*C))
        uAt<-uA*exp(uAe/(kb*C))
        uRt<-uR*exp(uRe/(kb*C))
        rt<-r*exp(-(C-Copt)^2/(2*rS)^2)
        qt<-(1/20)*C
        
        ca<-qt*Mt*(R/(Ht+R))
        cj<-(2-qt)*Mt*(R/(Ht+R))

        ifelse((sig*cj-tt)<0, mj<-0, mj<-(sig*cj-tt-uJt)/(1-z^(1-(uJt/(sig*cj-tt)))))
        ifelse((sig*ca-tt)<0, ra<-0, ra<-(sig*ca-tt)*B)
        
        dJ.dt<- ra*A -mj*J - uJt *J
        
        dA.dt<- mj*J - uAt*A -ra*p*A
        
        dR.dt<- rt*R*(1-(R/K)) - cj*J - ca*A - uRt*R
        
        return(list(c(dJ.dt,dA.dt,dR.dt)))
    })
}

#### SIMULATION ################################################################

# Set State Variables ----------------------------------------------------------
j.initial <- 1 # juveniles
a.initial <- 1 # adults
r.initial <- 1 # resources

# Set Parameters ---------------------------------------------------------------

# enter name of temperature parameter
temp.name <- "C"
# choose range of temperatures
temp.seq <- seq(Copt-10,Copt+10,length = 21)

# Set Plot Preferences ---------------------------------------------------------

# plot juvenile, adult, and resource dynamics?
# TODO add dynamics.to.plot <- c(10,15,20,25,30)
dynamics.plot = T
# plot bifurcation plot?
b.plot = T

# layout for plots
#TODO generalize
if (dynamics.plot & b.plot) {
    # first rows have dynamics plots
    # last row has bifurcation plots
    par(mfcol=c(6,1),
        mar=c(2,3,2,1), mgp=c(1.5,0.5,0))
} else if (dynamics.plot) {
    # rows have dynamics plots for one value in temp.seq
    par(mfcol=c(length(temp.seq),1),
        mar=c(2,3,2,1), mgp=c(1.5,0.5,0))
} else if (b.plot) {
    # one bifurcation plot
    par(mfcol=c(1,1),
        mar=c(2,3,2,1), mgp=c(1.5,0.5,0))
}

iterations <- 0:50000
iterations <- iterations/10
if (b.plot) {
    # which state variable should be plotted in bifurcation plot?
    b.var <- "JARatio"
    # which proportion of iterations should be plotted in bifurcation plot?
    # ex. 0.25 will plot last quarter of iterations
    b.portion <- 0.25
}

# Run Simulation ---------------------------------------------------------------

# create parameter vector
parms<-c(z=z, p=p, sig=sig, M=M, MS=MS, H=H,HS=HS, uJ=uJ, uJe=uJe, uA=uA, 
         uAe=uAe, uR=uR, uRe=uRe, t=t, te=te,r=r, rS=rS, K=K, B=B, kb=kb, C=C)
# find index of temperature to change
temp.index <- which(names(parms)==temp.name)

# initialize list to store output
output <- list()

# initialize state variables
y <- c(j.initial, a.initial, r.initial)
names(y) <- c("Juveniles", "Adults", "Resources")

# initialize parameters to change within loop
parms.loop <- parms

# loop to change temperature
for (j in 1:length(temp.seq)) {
    print(paste0("Starting simulation for ",temp.name,"=",temp.seq[j]))
    
    # set temperature
    parms.loop[temp.index] <- temp.seq[j]

    # run simulation
    output[[j]] <- ode(y=y, time=iterations, BaseStaget, parms=parms.loop)
    # remove time column
    output[[j]] <- output[[j]][,-1]
    # add juvenile to adult ratio
    JARatio <- output[[j]][,"Juveniles"] / output[[j]][,"Adults"]
    output[[j]] <- cbind(output[[j]], JARatio)
    
    # juvenile, adult, and resource dynamics plot
    if (dynamics.plot & sum(temp.seq[j]==c(10,15,20,25,30))==1) {
        matplot(output[[j]][,1:3],type="l",lty=1,pch=0.5,col=1:3, xlab='', ylab='')
        title(ylab=paste0(temp.name,"=",temp.seq[j]), cex.lab=1, font.lab=2)
    }

} # end temp.seq loop

# bifurcation plot
if (b.plot) {
    
    # choose rows to plot
    total.rows <- length(iterations)
    b.first.row <- round((1-b.portion) * total.rows)
    b.last.row <- total.rows
    b.rows <- b.first.row:b.last.row
    
    # set up plot axes and labels
    range.lim <- lapply(output, function(x) apply(x, 2, range))
    range.lim <- apply(do.call("rbind", range.lim), 2, range)
    plot(0, 0, pch = "", xlim = range(temp.seq), ylim = range.lim[,b.var],
         xlab='', ylab=b.var)
    
    # plot points
    for (j in 1:length(temp.seq)) {
        points(rep(temp.seq[j], length(b.rows)), output[[j]][b.rows,b.var])
    }
}

#### RUNTIME ###################################################################

# stop timing script
end <- proc.time()
# print elapsed time
print(end[3]-start[3])

