library("Rcpp")
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
a = 0.5		#inverse of searching time	
b = 2.0		#non-linear influence of speed difference
h = 2.0		#handling time for predation a prey

#following three parameters determine the prey traits
k = 0.1			#coefficient of foraging reward for prey (k*u is the reward)
mu_a = 0.0005	#basic mortality rate of active prey (not include the mortality by predation)
mu_r = 0.0005	#basic mortality rate of resting prey (mu_a >= mu_r)

#following two parameters determine the predator traits
c = 0.0001	#relative density of predator/prey
d = 0.01		#predation cost for predators

#Optimization
#Because prey have discrete choice (f = 0, f* or 1 where f* is the threshold prey frequency for predators activity),
#there are two potential optimal behaviour of prey (before and after the peak of fitness) 
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
Ans = tss_probforage_optimize(V, U,a, b, h, k, mu_a, mu_r, c, d)

#Optimal prey behaviour f (optimal predator behaviour is 1 if f=1, or zero otherwise)
#	grey color is the difference of upper and lower optimal points.
barplot(rbind(Ans$PreyL,Ans$PreyH-Ans$PreyL),xlab="time")

#Fitness of prey
Ans$PreyWL #lower point
Ans$PreyWH #upper point

#Fitness of predator
Ans$PredatorWL #lower point
Ans$PredatorWH #upper point

#Mortality rate at each time step
plot(Ans$PredationMortality_1,type="b",pch=19,xlab="time",ylab= "predation risk")

#Performance (i.e., reward/mortality) at each time step (blue:f < f*, red: f >= f*)
plot(Ans$drdm_f,type="b",ylim=c(0,max(Ans$drdm_f)),pch=19,col="blue",xlab="time",ylab= "dr/dm")
par(new=TRUE)
plot(Ans$drdm_1,type="b",ylim=c(0,max(Ans$drdm_f)),pch=19,col="red",xlab="time",ylab= "dr/dm")

