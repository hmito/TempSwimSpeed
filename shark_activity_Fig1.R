source("shark_activity_functions.R")
#=== constant parameters ===
#set time
pi = acos(-1)
tnum = 24 # time of day
t = 1:tnum-0.5

# temperature of the water
t_w=15
wmin = 25
wmax = 30

# amount of food availability for prey
alpha = rep(1.0, length=tnum)
# baseline mortality for prey (should pay both for resting and foraging)
mb = 0.2	
# average prey swim speed
u0 = 1.0
# metabolic cost for predators when they go out for predation
predcost=0.15
predCostMass=0.0 # due to endothermy

#=== default parameter values ===
omega = 1.0	#foraging efficiency increment by speed 
beta = 1.0 	#predation efficiency
phi = 0.1  	#probability of failing to hide in safe place
h = 1.0    	#handling time
r = 0.3  	#body radius (meter)
uk = 0.2 #influence of bodytemp
v0 = 1.5  #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp

#light effects
l0 = 1.0				#REQUIRE non-zero value! predation efficiency at noon
l1 = 1.0					#predation efficiency in perfect dark
light_influence = 5.0			#determines the sensitivity for small light
l_min=1;
twilight_coef=0.3;
#if light_influence is small, predation efficiency under twilight is mostly same with that in perfect dark.
#if light_influence is large, predation efficiency under twilight is mostly same with that at noon

#mortality rate of prey by predation
mx = 1.0	#predation by other predators
my = 0.5 #predation by sharks

watertemp=calc.watertemp(t,t_w,wmin,wmax)
L = calc.light_effect(t,l_min, light_influence, twilight_coef)

# assumption figure
par(mfrow=c(1,2))
plot(t,watertemp,col="blue",pch=15, xlab="time (t)",ylab="temperature")
lines(t,watertemp,col="blue",lty=1)
#	text(7,25.0,bquote('w'['t']))
points(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="red",pch=16)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="red",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,0.1),col="purple",pch=17)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.1),col="purple",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,0.6),col="orange",pch=18)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.6),col="orange",lty=1)
#	text(4,28,bquote('s'['t']))
#	par(new =T)
#	plot(t,L,col="orange",pch=17, type = "p", axes = FALSE, ylab = "",ylim=c(0,max(L)),main="temperature")
lines(t,L,col="orange",lty=1)
#	text(5,0.25,bquote(lambda['t']))
axis(4)

sharktemp=calc.sharktemp(t,t_w,wmin,wmax,r) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)

plot(t,V,type="n",xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),
	  xlab="time (t)",ylab="burst speed")
segments(-100,0,100,0)
points(t,V,col="red",pch=16)
lines(t,V,col="red",lty=1)
#	text(15,2.3,bquote('v'['t']))
points(t,U,col="blue",pch=15)
lines(t,U,col="blue",lty=1)
#	text(18,1.9,bquote('u'['t']))
points(t,V-U,col="purple",pch=17)
lines(t,V-U,col="purple",lty=1)