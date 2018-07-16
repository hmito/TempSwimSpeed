#TempSwimSpeed_v1_03
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

K = rep(0.01, length=tnum)		#amount of food availability for prey
C = rep(0.02, length=tnum)		#metaboric cost for predators when they go out for predation
L = 0.2+0.8*exp(-2.0*cos(2*pi*t/tnum))/exp(2.0)	#influence of light on the predation rate



#following three parameters determine the prey traits
e = 0.0		#relative risk of predation for resting prey
d = 0.001		#relative density of predator/prey
#[NEW 18/07/16] influence of swim speed on foraging efficiency
omega = 0	#obtained reward is alpha*((1-omega)+omega*U)
#following two parameters determine the predation rate: a*(v-u)^b / {1 + h*a*(v-u)^b} 
b = 3.0		#non-linear influence of speed difference
h = 2.0		#handling time for predation a prey

#following twp parameters determine the cost of prey
cb = 0.0001	#metabolic cost for prey (should pay both for resting and foraging)
cf	= 0.0001	#foraging cost for prey (should pay only for foraging)


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
#	Prey		Optimal prey behaviour
#	PreyW		Fitness of prey
#	Predator	Optimal predator behaviour
#	PredatorW	Fitness of predator
#	ThresholdPreyFreq	The threshold prey frequency for predators activity f* 
#	PreyReward	Reward of prey (independent from the strategy of prey/predators)
#	PreyCost0	Cost of prey caused by predation when f=0
#	PreyCostF	Cost of prey caused by predation when f=f*
#	PreyCost1	Cost of prey caused by predation f=1
Ans = tss_probforage_energygain_optimize(V, U, K, C, L, d, e, omega, b, h,cb,cf)


plot(0,0,type="n",xlim=c(0,tnum),ylim=c(0,max(c(vmax,umax))))
lines(t,Ans$Predator,col="red",lwd=3)
lines(t,Ans$Prey,col="blue",lwd=3)

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
barplot(rbind(Ans$Prey,Ans$Prey-Ans$Prey),xlab="time")
points(1:length(Ans$Prey)+0.5,Ans$ThresholdPreyFreq)

#Optimal predatpr behaviour p
#	grey color is the difference of upper and lower optimal points.
barplot(rbind(Ans$Predator),xlab="time")

#Mortality rate at each time step
plot(rep(1:length(Ans$PreyCost),times=2),c(Ans$PreyCost,Ans$PredatorCost),type="n",pch=19,xlab="time",ylab= "predation risk", col="red")
#lines(Ans$PreyCostF,type="b",pch=19,xlab="time",ylab= "predation risk", col="green")
lines(Ans$PreyCost,type="b",pch=19,xlab="time",ylab= "predation risk",col="red")
lines(Ans$PreyReward,type="b",pch=19,xlab="time",ylab= "predation risk",col="blue")

plot(Ans$PredatorReward,type="b",pch=19,xlab="time",ylab= "predation risk",col="blue")
lines(Ans$PredatorCost,type="b",pch=19,xlab="time",ylab= "predation risk",col="red")

