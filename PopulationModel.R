
library(deSolve)

#### RUNTIME ###################################################################

# start timing script
start <- proc.time()

#### PARAMETERS ################################################################

# z: adult to juvenile size ratio
z<-0.02
# p: adult modifier on production of juveniles
# can use to modify stage structure, seems like it needs to be <1
p<-0.06
# sig: conversion efficiency of resource to useful energy for fish
sig<-0.75
# M: maximum ingestion rate
# make mass specific?
M<-0.5
# H: andling time - speed at which fish can eat resources, smaller is faster
H<-3
# uJ: juvenile mortality
uJ<-0.005
# uA: adult mortality
uA<-0.005
# uR: resource mortality
# will be useful when we add temp dependence
uR<-0.005
# t: costs of maintaining somatic growth/turnover
# i.e. base level of resource intake you must exceed to mature/reproduce
t<-0.01
# r: resource growth rate
r<-5
# K: resource carrying capacity
K<-50
# B: adult reproductive rate
B<-0.5

# have to tell desolve which parameters to care about
parms<-c(z=z, p=p, sig=sig, M=M, H=H, uJ=uJ, uA=uA, uR=uR, t=t, r=r, K=K, B=B)

#### ODE FUNCTIONS #############################################################

# ca: functional response for adults
# cj: functional response for juveniles
# mj: juvenile maturation rate
# ra: reproduction per adult
# uA: mortality rate for adults
# uJ: mortality rate for uveniles

BaseStage<-function(t,y,p){
	{
		J<-y[1]
		A<-y[2]
		R<-y[3]
	}
	with(as.list(p),{
		ca<-M*(R/(H+R))
		cj<-M*(R/(H+R))
		ifelse(((sig*cj)-t)<0, mj<-0, mj<-(((sig*cj)-t)-uJ)/(1-z^(1-(uJ/((sig*cj)-T)))))
		ifelse(((sig*ca)-t)<0, ra<-0, ra<-((sig*ca)-t)*B)
		
		dJ.dt<- ra*A -mj*J - uJ
		
		dA.dt<- mj*J - uA -ra*p*A
		
		dR.dt<- r*R*(1-(R/K)) - cj*J - ca*A -uR
		
        return(list(c(dJ.dt,dA.dt,dR.dt)))
	})
}

#### SIMULATION ################################################################

# state variable initial conditions
J<-1
A<-5
R<-10
y <- c(J,A,R)
names(y) <- c("Juveniles", "Adults", "Resources")

# duration of simulation
end.time<-1000

# define how time works for simulation
days<-(seq(0,end.time,by=0.1))

# run desolve to simulate the model through time (days)
BS.out<-data.frame(ode(y=y,time=days,func=BaseStage, parms=parms))

matplot(BS.out[,2:4],type="l",lty=1,pch=0.5,col=1:3)
legend('right', names(y), lty=1,col=1:3, bty = "n")

#### RUNTIME ###################################################################

# stop timing script
end <- proc.time()
# print elapsed time
print(end[3]-start[3])
