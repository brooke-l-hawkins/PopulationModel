
library(deSolve)

#### RUNTIME ###################################################################

# start timing script
start <- proc.time()

#### VARIABLES #################################################################

# Parameters -------------------------------------------------------------------

z    <- 0.2         # adult to juvenile size ratio
p    <- 0.4         # adult modifier on production of juveniles
                    # can use to modify stage structure
                    # seems like it needs to be <1
sig  <- 0.7         # conversion efficiency of resource to useful energy for fish
M    <- 1           # maximum ingestion rate
                    # make mass-specific?
                    # Mt vs. temperature forms a negative parabola
MS   <- 10          # function breadth for max intake rate
                    # attack rate
H    <- 1           # handling time
                    # speed at which fish can eat resources (smaller is faster)
HS   <- 10          # function breadth for handling time
                    # HSt vs. temperature forms a negative parabola
uJ   <- 10000       # juvenile mortality
                    # uJt = 0.05 when C = 293.15
uA   <- 10000       # adult mortality
                    # uAt = 0.05 when C = 293.15
uJe  <- -0.308      # activation energy of juvenile mortality
uAe  <- -0.308      # activation energy of adult mortality
uR1  <- 1000        # shared resource mortality
uR2  <- 1000        # adult-specific resource mortality
                    # uRt = 0.005 when C = 293.15
uR1e <- -0.308      # activation energy of shared resource mortality
uR2e <- -0.308      # activation energy of adult-specific resource mortality
t    <- 10000       # costs of maintaining somatic growth/turnover
                    # base level of resource intake you must exceed to mature/reproduce
                    # tt = 0.1 when C = 293.15
te   <- -0.291      # activation energy of metabolic waste
r1   <- 1.5         # shared resource growth rate
r2   <- 1.5         # adult-specific resource growth rate
                    # r1t and r2t vs. temperature forms a negative parabola
r1S  <- 15          # function breadth for shared resource growth rate
r2S  <- 15          # function breadth for adult-specific resource growth rate
r1p  <- 1           # proportion of shared resource available to adults
r2p  <- 0.5         # proportion of adult-specific resource available to adults
                    # r2p is a fraction between 0 and 1 (inclusive)
                    # r2p=0   when adults access none of adult-specific resource
                    # r2p=0.5 when adults access half of adult-specific resource
                    # r2p=1   when adults access all of adult-specific resource
K1   <- 5           # shared resource carrying capacity
K2   <- 5           # adult-specific resource carrying capacity
B    <- 0.5         # adult reproductive rate
kb   <- 8.617*10^-5 # boltzmann's constant
C    <- 293.15      # temperature in kelvin
                    # 283.15 is cold, 293.15 is optimal, 303.15 is hot
Copt <- 293.15      # optimal temperature in kelvin

# Initial State Variables ------------------------------------------------------

j.initial  <- 1     # juveniles
a.initial  <- 1     # adults
r1.initial <- 1     # shared resources
r2.initial <- 1     # adult-specific resources

# Create Variable Vector -------------------------------------------------------

variable.vec <- c(z=z, p=p, sig=sig, M=M, MS=MS, H=H, HS=HS, uJ=uJ, uA=uA,
                  uJe=uJe, uAe=uAe, uR1=uR1, uR2=uR2, uR1e=uR1e, uR2e=uR2e, t=t,
                  te=te, r1=r1, r2=r2, r1S=r1S, r2S=r2S, r1p=r1p, r2p=r2p,
                  K1=K1, K2=K2, B=B, kb=kb, C=C, Copt=Copt, j.initial=j.initial,
                  a.initial=a.initial, r1.initial=r1.initial,
                  r2.initial=r2.initial)

# enter indices of parameters in variable.vec
parameters.indices <- 1:29

# enter indices of initial state variables in variable.vec
state.indices <- 30:33

#### CHANGING VARIABLES ########################################################

# enter name of temperature parameter
C.name <- "C"
# choose range of temperatures
C.seq <- seq(from=Copt-10, to=Copt+10, length=5)

# enter name of proportion of adult-specific resource available to adults
r2p.name <- "r2p"
# choose range of proportions
r2p.seq <- seq(from=r2p-0.25, to=r2p+0.25, length=3)

# determine repetition length
repetitions.length <- length(C.seq)*length(r2p.seq)

# create parameter matrix
variable.matrix <- matrix(data=variable.vec, nrow=length(variable.vec),
                          ncol=repetitions.length)
row.names(variable.matrix) <- names(variable.vec)

# update temperature
variable.matrix[C.name,] <- rep(C.seq, times=length(r2p.seq))

# update proportions
variable.matrix[r2p.name,] <- rep(r2p.seq, each=length(C.seq))

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
        qt<-1
        
        ca1<-qt*Mt*r1p*R1/(Ht+r1p*(R1+r2p*R2))
        ca2<-qt*Mt*r2p*R2/(Ht+r2p*(R1+r2p*R2))
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

# Dynamics Plot Preferences ----------------------------------------------------

# plot juvenile, adult, and resource dynamics?
dynamics.plot = T

if (dynamics.plot){
    # which dynamics to plot?
    C.dynamics.to.plot <- C.seq[c(-2,-4)]
    r2p.dynamics.to.plot <- r2p.seq
    # plot layout
    par(mfcol=c(length(C.dynamics.to.plot),length(r2p.dynamics.to.plot)),
        mar=c(2,3,2,1), mgp=c(1.5,0.5,0))
}

# Run Simulation ---------------------------------------------------------------

# initialize list to store output
output <- list()

for (r in 1:repetitions.length) {
    
    print(paste0("Starting simulation for repetition ",r))
    
    # set parameters
    parms <- variable.matrix[parameters.indices, r]

    # set initial state variables
    y <- variable.matrix[state.indices, r]
    names(y) <- c("Juveniles", "Adults", "Shared Resources", "Adult Resources")
    
    # run simulation
    output[[r]] <- ode(y=y, time=iterations, BaseStaget, parms=parms)
    # remove time column
    output[[r]] <- output[[r]][,-1]
    # add juvenile to adult ratio
    JARatio <- output[[r]][,"Juveniles"] / output[[r]][,"Adults"]
    output[[r]] <- cbind(output[[r]], JARatio)
    
    # juvenile, adult, and resources dynamics plot
    if (dynamics.plot) {
        # if parameter is included in parameters for dynamics to plot
        if (sum(variable.matrix[C.name,r]==C.dynamics.to.plot)==1 &
            sum(variable.matrix[r2p.name,r]==r2p.dynamics.to.plot)==1) {
            # plot dynamics
            matplot(output[[r]][,1:4],type="l",lty=1,pch=0.5,col=1:4, xlab='', ylab='')
            title(ylab=paste0(C.name,"=",variable.matrix[C.name,r]), 
                  main=paste0(r2p.name,"=",variable.matrix[r2p.name,r]),
                  cex.lab=1, font.lab=2)
        }
    }

} # end repetitions loop

# Bifurcation Plots ------------------------------------------------------------

# plot bifurcation plot?
b.plot = T

# bifurcation plot
if (b.plot) {

    # which state variable should be plotted in bifurcation plot?
    b.vars <- c("JARatio", "Juveniles", "Adults")
    # which proportion of iterations should be plotted in bifurcation plot?
    # ex. 0.25 will plot last quarter of iterations
    b.portion <- 0.5
    
    # plot layout
    par(mfcol=c(length(b.vars),1), mar=c(2,3,2,1), mgp=c(1.5,0.5,0))
    
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
        plot(0, 0, pch = "", xlim = range(variable.matrix[C.name,]), ylim = range.lim[,b],
             xlab='', ylab=b)
        
        # plot points
        for (r in 1:repetitions.length) {
            max.val <- max(output[[r]][b.rows,b])
            min.val <- min(output[[r]][b.rows,b])
            # which.r2p used to assign colors and offset lines in plot
            which.r2p <- which(r2p.seq==variable.matrix[r2p.name,r])
            lines(x=rep(variable.matrix[C.name,r], 2)+which.r2p*0.1,
                  y=c(min.val, max.val), col=which.r2p)
        }
    }
    
}

#### RUNTIME ###################################################################

# stop timing script
end <- proc.time()
# print elapsed time
print(end[3]-start[3])

