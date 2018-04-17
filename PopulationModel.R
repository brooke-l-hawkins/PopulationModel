library(deSolve)

#parameters to define

#adult to juvenile size ratio
z<-0.02
#adult modifier on production of juveniles - can use to modify stage structure, seems like it needs to be less than 1
p<-0.06
#conversion efficiency of resource to useful energy for fish
sig<-0.75
#Maximum ingestion rate ~ make mass specific?
M<-0.5
#Handling time - speed at which fish can eat resources, smaller is faster
H<-3
#Mortality for juveniles
uJ<-0.005
#Mortality for adults
uA<-0.005
#Mortality for resources (will be useful when we add temp dependence)
uR<-0.005
#costs of maintaining somatic growth/turnover i.e. the base level of resource intake you have to exceed to mature or reproduce
t<-0.01
#Resource growth rate
r<-5
#Resource carrying capacity
K<-50
#Adult reproductive rate parameter
B<-0.5

#ODE functions for this model are defined below

#ra is reproduction per adult
#mj is the juvenile maturation rate
#u with a stage following it is the mortality rate for that stage
#ca and cj are the functional responses for adults and juveniles respectively (though they are identical right now)
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

#have to tell desolve which parameters to care about
parms<-c(M=M, H=H, sig=sig, t=t, uJ=uJ, B=B, uA=uA,p=p ,r=r, K=K, uR=uR)


#Initial conditions for state variable
J<-1
A<-5
R<-10



#duration of simulation
end.time<-1000

#defining how time works for simulation
days<-(seq(0,end.time,by=0.1))

#running desolve to simulate the model through time (days)
BS.out<-data.frame(ode(y=c(J,A,R),time=days,func=BaseStage, parms=parms))

matplot(BS.out[,2:4],type="l",lty=1,pch=0.5,col=1:3)
legend('right', c("Resource","Adults","Juveniles"), lty=1,col=3:1, bty = "n")