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
C = rep(0.01,length=length(t))	#relative density of predator/prey
base_c =0.1
d = 0.01		#predation cost for predators


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
#	ThresholdPreyFreq	The threshold prey frequency for predators activity	f* 
#	PreyReward	Reward of prey (independent from the strategy of prey/predators)
#	PredationMortality_1	Mortality rate of prey caused by predation when f=1
#	drdm_f	Performance of foraging for prey when f < f*
#	drdm_1	Performance of foraging for prey when f = 1
Ans = tss_ppgame_optimize(V, U, C, base_c, a, b, h, k, d, e)
Mutant = tss_ppgame_stability(Ans$Prey,V, U, C, base_c, a, b, h, k, d, e, 10000000)

Ans2 = tss_ppgame_optimize_hill_climb(V, U, C, base_c, a, b, h, k, d, e, 10000000)

Ans$PreyW
Ans2$PreyW

Prey = Ans$Prey
Prey[7]=1

Ans2 = tss_ppgame_fitness(Prey,V, U, C, base_c, a, b, h, k, d, e)

#Optimal prey behaviour f (optimal predator behaviour is 1 if f=1, or zero otherwise)
#	grey color is the difference of upper and lower optimal points.
barplot(Ans$Prey,xlab="time")
barplot(Mutant,xlab="time")
barplot(Ans$Predator,xlab="time")

barplot(rbind(Ans$PredatorR,Ans$PredatorH-Ans$PredatorL),xlab="time")
Ans$PreyWH
Ans$PreyWL
Ans$PredatorWH
Ans$PredatorWL

#Mortality rate at each time step
plot(U/Ans$PreyMortality0,type="b",pch=19,xlab="time",ylab= "predation risk")
plot(U/Ans$PreyMortalityF,type="b",pch=19,xlab="time",ylab= "predation risk")
plot(U/Ans$PreyMortality1,type="b",pch=19,xlab="time",ylab= "predation risk")

#Performance (i.e., reward/mortality) at each time step (blue:f < f*, red: f >= f*)
plot(Ans$PredatorR,type="b",ylim=c(0,max(Ans$PredatorR)),pch=19,col="blue",xlab="time",ylab= "dr/dm")
par(new=TRUE)
plot(Ans$PreyR,type="b",ylim=c(0,max(Ans$PredatorR)),pch=19,col="red",xlab="time",ylab= "dr/dm")
par(new=TRUE)
plot(Ans$PreyM,type="b",ylim=c(0,max(Ans$PredatorR)),pch=19,col="green",xlab="time",ylab= "dr/dm")


