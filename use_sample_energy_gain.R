#TempSwimSpeed_v1_02
#   Predator efficiency 

library("Rcpp")
library("BH")
sourceCpp("TempSwimSpeedOptim.cpp")

pi = acos(-1)

tnum = 24
t = 0:23+0.5

vmin = 1.5
vmax = 2.5
vmaxt = 18

umin = 0.5
umax = 2.0
umaxt = 15

#Define the speed of predator (V) and prey (U) at each time step
V = vmin + (vmax - vmin)*(cos(2*pi*(t-vmaxt)/tnum)+1)/2
U = umin + (umax - umin)*(cos(2*pi*(t-umaxt)/tnum)+1)/2

K = rep(0.02, length=tnum)		#amount of food availability for prey
C = rep(0.02, length=tnum)		#metaboric cost for predators when they go out for predation
L = 0.2+0.8*exp(-2.0*cos(2*pi*t/tnum))/exp(2.0)	#influence of light on the predation rate



#following three parameters determine the prey traits
e = 0.3		#relative risk of predation for resting prey
d = 0.001	#relative density of predator/prey
#following three parameters determine the predation rate: a*(v-u)^b / {1 + h*a*(v-u)^b} 
b = 3.0		#non-linear influence of speed difference
h = 2.0	#handling time for predation a prey

#Variation of model
#if we set d=0, the prey decide their behaviour only based on the food availability for them regardless of the predation.
# d = 0.0
#
#if we set b=0, the swim speed (i.e., temperature) has no influence on the predation rate.
# b = 0.0
#
#if fix the V and U as constant values, the influence of temperature is ignored in this model.
# V=seq(1.5, length=tnum)
# U=seq(1.2, length=tnum)

#Optimization
#Because prey have discrete choice (f = 0, f* or 1 where f* is the threshold prey frequency for predators activity),
#there are two potential optimal behaviour of prey (before and after the peak of fitness)
#Return values include the following members
#	PreyL		Optimal prey behaviour at lower point of the fitness peak
#	PreyH		Optimal prey behaviour at upper point of the fitness peak
#	PreyWL	Fitness of prey at PreyL
#	PreyWH	Fitness of prey at PreyH
#	PredatorL	Optimal predator behaviour at PreyL
#	PredatorH	Optimal predator behaviour at PreyH
#	PredatorWL	Fitness of predator at PreyL
#	PredatorWH	Fitness of predator at PreyH
#	ThresholdPreyFreq	The threshold prey frequency for predators activity f* 
#	PreyReward	Reward of prey (independent from the strategy of prey/predators)
#	PreyMortality0	Mortality rate of prey caused by predation when f=0
#	PreyMortalityF	ortality rate of prey caused by predation when f=f*
#	PreyMortality1	ortality rate of prey caused by predation f=1
Ans = tss_probforage_energygain_optimize(V, U, K, C, L, d, e, b, h)


#Summarized figure
if(Ans$PreyWL>Ans$PreyWH){
	Prey= Ans$PreyL
	Predator = Ans$PredatorL
}else{
	Prey= Ans$PreyH
	Predator = Ans$PredatorH
}
plot(0,0,type="n",xlim=c(0,tnum),ylim=c(0,max(c(vmax,umax))))
lines(t,Predator,col="red",lwd=3)
lines(t,Prey,col="blue",lwd=3)

points(t,V,col="red")
lines(t,V,col="red")
points(t,U,col="blue")
lines(t,U,col="blue")
points(t,L,col="black")
lines(t,L,col="black")

lines(t,C,col="red",lty=3)
lines(t,K,col="blue",lty=3)


#Optimal prey behaviour f
#	grey color is the difference of upper and lower optimal points.
barplot(rbind(Ans$PreyL,Ans$PreyH-Ans$PreyL),xlab="time")
points(1:length(Ans$PreyL)+0.5,Ans$ThresholdPreyFreq)

#Optimal predatpr behaviour p
#	grey color is the difference of upper and lower optimal points.
barplot(rbind(Ans$PredatorL,Ans$PredatorH-Ans$PredatorL),xlab="time")

#Mortality rate at each time step
plot(Ans$PreyMortality0,type="b",pch=19,xlab="time",ylab= "predation risk")
plot(Ans$PreyMortalityF,type="b",pch=19,xlab="time",ylab= "predation risk")
plot(Ans$PreyMortality1,type="b",pch=19,xlab="time",ylab= "predation risk")

