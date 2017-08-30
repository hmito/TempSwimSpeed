#TempSwimSpeed_v1_01
#   Predator 

library("Rcpp")
library("BH")
sourceCpp("TempSwimSpeedOptim.cpp")

pi = acos(-1)

vmin = 0.5
vmax = 2.5
umin = 0.8
umax = 2.0
n = 0.5

t = 0:24

#Define the speed of predator (V) and prey (U) at each time step
V = vmin + (vmax - vmin)*sin(pi*t/24)^n
U = umin + (umax - umin)*sin(pi*t/24)

#following three parameters determine the predation rate: a*(v-u)^b / {1 + h*a*(v-u)^b} 
a = 1		#inverse of searching time	
b = 2.0		#non-linear influence of speed difference
h = 2.0		#handling time for predation a prey

#following three parameters determine the prey traits
k = 0.01		#coefficient of foraging reward for prey (k*u is the reward)
e = 0.3		#relative risk of predation for resting prey
	
#following two parameters determine the predator traits
C = rep(0.01,length=length(t))	#relative density of predator/prey
d = 0.001		#predation cost for predators

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
Ans = tss_probforage_energygain_optimize(V, U, C, a, b, h, k, d, e)

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
