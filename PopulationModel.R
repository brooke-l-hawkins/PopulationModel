
library(deSolve)

#### RUNTIME ###################################################################

# start timing script
start <- proc.time()

#### PARAMETERS ################################################################

# z: adult to juvenile size ratio
z<-0.2
# p: adult modifier on production of juveniles
# can use to modify stage structure, seems like it needs to be <1
p<-0.4
# sig: conversion efficiency of resource to useful energy for fish
sig<-0.7
# M: maximum ingestion rate
# make mass specific?
M<-1
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
# uAt = 0.05 when C = 293.15
uA<-10000
# uAe: activation energy of adult mortality
uAe<--0.308
# uR1: shared resource mortality
# uR2: adult-specific resource mortality
# uR = 0.005 when C = 293.15
uR1<-uR2<-1000
# uR1e: activation energy of shared resource mortality
# uR2e: activation energy of adult-specific resource mortality
uR1e<-uR2e<--0.308
# t: costs of maintaining somatic growth/turnover
# i.e. base level of resource intake you must exceed to mature/reproduce
# t = 0.1 when C = 293.15
t<-10000
# te: activation energy of metabolic waste
te<--0.291
# r1: shared resource growth rate
# r2: adult-specific resource growth rate
r1<-r2<-1.5
# r1S: function breadth for shared resource growth rate
# r2S: function breadth for adult-specific resource growth rate
# temperature parabola opens downward
r1S<-r2S<-15
# r2p: adult preference for adult-specific resource
# fraction between 0 and 1 (inclusive)
    # r2p=0 when adults ignore adult-specific resource
    # r2p=0.5 when adults prefer shared and adult-specific resources equally
    # r2p=1 when adults ignore shared resource
r2p<-0
#r1p: adult preference for shared resource
# fraction between 1 and 0 (inclusive)
r1p<- 1-r2p
# K1: shared resource carrying capacity
# K2: adult-specific resource carrying capacity
K1<-K2<-5
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
# uR1t: mortality rate for shared resources
# uR2t: mortality rate for adult-specific resources
# r1t: resource growth rate
# r2t: adult-specific resource growth rate
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
        R1<-y[3]
        R2<-y[4]
    }
    with(as.list(p),{
        Mt<-M*exp(-(C-Copt)^2/(2*MS)^2)
        Ht<-H*exp((C-Copt)^2/(2*HS)^2)
        tt<-t*exp(te/(kb*C))
        uJt<-uJ*exp(uJe/(kb*C))
        uAt<-uA*exp(uAe/(kb*C))
        uR1t<-uR1*exp(uR1e/(kb*C))
        uR2t<-uR2*exp(uR2e/(kb*C))
        r1t<-r1*exp(-(C-Copt)^2/(2*r1S)^2)
        r2t<-r2*exp(-(C-Copt)^2/(2*r2S)^2)
        qt<--0.01*(C-Copt)^2+1.5
        
        ca1<-qt*Mt*r1p*R1/(Ht+r1p*(R1+R2))
        ca2<-qt*Mt*r2p*R2/(Ht+r2p*(R1+R2))
        ca<-ca1+ca2
        cj<-(2-qt)*Mt*R1/(Ht+R1)

        ifelse((sig*cj-tt)<0, mj<-0, mj<-(sig*cj-tt-uJt)/(1-z^(1-(uJt/(sig*cj-tt)))))
        ifelse((sig*ca-tt)<0, ra<-0, ra<-(sig*ca-tt)*B)
        
        dJ.dt<- ra*A - mj*J - uJt*J
        dA.dt<- mj*J - uAt*A - ra*p*A
        dR1.dt<- r1t*R1*(1-(R1/K1)) - ca1*A - uR1t*R1 - cj*J
        dR2.dt<- r2t*R2*(1-(R2/K2)) - ca2*A - uR2t*R2
        
        return(list(c(dJ.dt, dA.dt, dR1.dt, dR2.dt)))
    })
}

#### SIMULATION ################################################################

# Choose How Long to Run -------------------------------------------------------
iterations <- 0:25000
iterations <- iterations/10

# Set State Variables ----------------------------------------------------------
j.initial <- 1 # juveniles
a.initial <- 1 # adults
r1.initial <- 1 # shared resources
r2.initial<- 1 # adult-specific resources

# Set Parameters ---------------------------------------------------------------

# enter name of temperature parameter
temp.name <- "C"
# choose range of temperatures
temp.seq <- seq(Copt-10,Copt+10,length = 21)

# Set Plot Preferences ---------------------------------------------------------

# plot juvenile, adult, and resource dynamics?
dynamics.plot = T
if (dynamics.plot){
    # which dynamics to plot?
    dynamics.to.plot <- c(Copt-10, Copt-5, Copt, Copt+5, Copt+10)  
}

# plot bifurcation plot?
b.plot = T
if (b.plot) {
    # which state variable should be plotted in bifurcation plot?
    b.vars <- c("JARatio", "Juveniles", "Adults")
    # which proportion of iterations should be plotted in bifurcation plot?
    # ex. 0.25 will plot last quarter of iterations
    b.portion <- 0.5
}

# layout for plots
if (dynamics.plot & b.plot) {
    # starting rows have dynamics plots
    # last rows have bifurcation plots
    par(mfcol=c(length(dynamics.to.plot)+length(b.vars),1),
        mar=c(2,3,2,1), mgp=c(1.5,0.5,0))
} else if (dynamics.plot) {
    # rows have dynamics plots for one value in temp.seq
    par(mfcol=c(length(dynamics.to.plot),1),
        mar=c(2,3,2,1), mgp=c(1.5,0.5,0))
} else if (b.plot) {
    # rows have bifurcation plots
    par(mfcol=c(length(b.vars),1),
        mar=c(2,3,2,1), mgp=c(1.5,0.5,0))
}

# Run Simulation ---------------------------------------------------------------

# create parameter vector
parms<-c(z=z, p=p, sig=sig, M=M, MS=MS, H=H, HS=HS, uJ=uJ, uJe=uJe, uA=uA, uAe=uAe,
         uR1=uR1, uR2=uR2, uR1e=uR1e, uR2e=uR2e, t=t, te=te, r1=r1, r2=r2,
         r1S=r1S, r2S=r2S, K1=K1, K2=K2, B=B, kb=kb, C=C, Copt=Copt)
# find index of temperature to change
temp.index <- which(names(parms)==temp.name)

# initialize list to store output
output <- list()

# initialize state variables
y <- c(j.initial, a.initial, r.initial, ar.initial)
names(y) <- c("Juveniles", "Adults", "Shared Resources", "Adult Resources")

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
    if (dynamics.plot) {
        if (sum(temp.seq[j]==dynamics.to.plot)==1) {
            matplot(output[[j]][,1:4],type="l",lty=1,pch=0.5,col=1:4, xlab='', ylab='')
            title(ylab=paste0(temp.name,"=",temp.seq[j]), cex.lab=1, font.lab=2)
        }
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
    
    for (b in b.vars) {
        # set up axes
        plot(0, 0, pch = "", xlim = range(temp.seq), ylim = range.lim[,b],
             xlab='', ylab=b)
        
        # plot points
        for (j in 1:length(temp.seq)) {
            points(rep(temp.seq[j], length(b.rows)), output[[j]][b.rows,b])
        }
    }
    
}

#### RUNTIME ###################################################################

# stop timing script
end <- proc.time()
# print elapsed time
print(end[3]-start[3])

