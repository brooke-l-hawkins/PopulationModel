
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
# uAt = 0.05 when C = 293.15
uA<-10000
# uAe: activation energy of adult mortality
uAe<--0.308
# uR: resource mortality
# uR = 0.005 when C = 293.15
uR<-1000
# uRe: activation energy of resource mortality
uRe<--0.308
# t: costs of maintaining somatic growth/turnover
# i.e. base level of resource intake you must exceed to mature/reproduce
# t = 0.1 when C = 293.15
t<-10000
# te: activation energy of metabolic waste
te<--0.291
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
# Copt: optimal temperature in kelvin
# 283.15 is cold, 293.15 is optimal, 303.15 is hot
Copt<-293.15

# C: temperature in kelvin
# 283.15 is cold, 293.15 is optimal, 303.15 is hot
C<-seq(from=Copt-10,to=Copt+10,length = 21)

# initial state values
R <- J <- A <- 1

# plotting settings
par(mfcol=c(3,2), mar=c(4,3,3,1))

#### MODEL #####################################################################

# t subscript in variable name indicates temperature-sensitive
# base terms are now the scaling factor for temperature-sensitive terms

# Mt: maximum ingestion rate
# Ht: handling time
# tt: costs of maintaining somatic growth/turnover
# uJt: mortality rate for juveniles
# uAt: mortality rate for adults
# uRt: mortality rate for resources
# rt: resource growth rate
# qt: competitive difference between adults and juveniles
# ca: functional response for adults
# cj: functional response for juveniles
# mj: juvenile maturation rate
# ra: reproduction per adult

Mt<-M*exp(-(C-Copt)^2/(2*MS)^2)
Ht<-H*exp((C-Copt)^2/(2*HS)^2)
tt<-t*exp(te/(kb*C))
uJt<-uJ*exp(uJe/(kb*C))
uAt<-uA*exp(uAe/(kb*C))
uRt<-uR*exp(uRe/(kb*C))
rt<-r*exp(-(C-Copt)^2/(2*rS)^2)

qt<--0.01*(C-Copt)^2+1.5

ca<-qt*Mt*(R/(Ht+R))
cj<-(2-qt)*Mt*(R/(Ht+R))

mj.condition <- !(sig*cj-tt)<0
mj <- (sig*cj-tt-uJt)/(1-z^(1-(uJt/(sig*cj-tt))))
mj <- mj.condition * mj

ra.condition <- !(sig*ca-tt)<0
ra <- (sig*ca-tt)*B
ra <- ra.condition * ra

dJ.dt<- ra*A -mj*J - uJt *J
dA.dt<- mj*J - uAt*A -ra*p*A
dR.dt<- rt*R*(1-(R/K)) - cj*J - ca*A - uRt*R

#### PLOTS #####################################################################

plot(qt~C, type='l', col=1, main="competitive difference", ylim = range(qt))

plot(Mt~C, type='l', col=4, main="physiology rates",
     ylim=range(Mt, Ht, tt))
lines(C, Ht, col=5)
lines(C, tt, col=6)
legend("topleft", legend=c("Mt", "Ht", "tt"),
       col=4:6, lty=1)

plot(uJt~C, type='l', col=1, main="mortality rates",
     ylim=range(uJt, uAt, uRt))
lines(C, uAt, col=2)
lines(C, uRt, col=3)
legend("topleft", legend=c("uJt", "uAt", "uRt"),
       col=1:3, lty=1)

plot(dJ.dt~C, type='l', col=1, main="growth rates", 
     ylim=range(dJ.dt, dA.dt, dR.dt))
lines(C, dA.dt, col=2)
lines(C, dR.dt, col=3)
legend("topleft", legend=c("dJ.dt", "dA.dt", "dR.dt"),
       col=1:3, lty=1)

plot(cj~C, type='l', col=1, main="functional response",
     ylim=range(cj,ca))
lines(C, ca, col=2)
legend("topleft", legend=c("cj", "ca"),
       col=1:2, lty=1)

plot(mj~C, type='l', col=1, main="stage specific rates",
     ylim=range(mj,ra))
lines(C, ra, col=2)
legend("topleft", legend=c("mj", "ra"),
       col=1:2, lty=1)




