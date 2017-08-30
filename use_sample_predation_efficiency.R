#TempSwimSpeed_v1_01
#   Predator efficiency 

library("Rcpp")
library("BH")
sourceCpp("TempSwimSpeedOptim.cpp")

pi = acos(-1)

vmin = 1.0
vmax = 2.5
umin = 1.0
umax = 2.0
n = 2

t = 0:12

#Define the speed of predator (V) and prey (U) at each time step
V = vmin + (vmax - vmin)*sin(pi*t/24)^n
U = umin + (umax - umin)*sin(pi*t/24)

#following three parameters determine the predation rate: a*(v-u)^b / {1 + h*a*(v-u)^b} 
a = 1		#inverse of searching time	
b = 2.0		#non-linear influence of speed difference
h = 2.0		#handling time for predation a prey

#following three parameters determine the prey traits
k = 0.0011		#coefficient of foraging reward for prey (k*u is the reward)
e = 0.01		#relative risk of predation for resting prey

#following two parameters determine the predator traits
C = rep(0.01,length=length(t))	#metaboric cost for predators when they go out for predation
base_c =0.1								#basic metaboric cost for predators
d = 0.01									#relative density of predator/prey

#Optimization in prey-predator game model
#Return values include the following members
#	Prey			Optimal prey behaviour
#	PreyW			Fitness of prey 
#	Predator		Optimal predator behaviour at PreyL
#	PredatorW	Fitness of predator
#	PreyR			Reward of prey
#	PreyM			Mortality rate of prey
#	PredatorR	Reward (Predation benefit) of predators
#	PredatorC	Metaboric cost of predators
Ans = tss_probforage_predeff_optimize(V, U, C, base_c, a, b, h, k, d, e)

#Optimal prey behaviour f
barplot(Ans$Prey,xlab="time")

#Optimal predator behaviour f
barplot(Ans$Predator,xlab="time")

#Performance (i.e., reward/mortality) at each time step (blue:f < f*, red: f >= f*)
plot(Ans$PredatorR,type="b",ylim=c(0,max(Ans$PredatorR)),pch=19,col="blue",xlab="time",ylab= "dr/dm")
par(new=TRUE)
plot(Ans$PreyR,type="b",ylim=c(0,max(Ans$PredatorR)),pch=19,col="red",xlab="time",ylab= "dr/dm")
par(new=TRUE)
plot(Ans$PreyM,type="b",ylim=c(0,max(Ans$PredatorR)),pch=19,col="green",xlab="time",ylab= "dr/dm")


