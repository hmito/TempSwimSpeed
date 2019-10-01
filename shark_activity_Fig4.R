#load library
source("shark_activity_functions.R")


#=== constant parameters ===
#set time
pi = acos(-1)
tnum = 24 # time of day
t = 1:tnum-0.5

# temperature of the water
t_w = 15
wmin = 25
wmax = 30
watertemp = wmin+(wmax-wmin)*(cos(2*pi*(t-t_w)/length(t))+1)/2

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

#prey and predator speed
r = 0.3 #effective body size in the context of heat balance
uk = 0.2 #influence of bodytemp
v0 = 1.5  #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp

#light effects
l0 = 1.0				#REQUIRE non-zero value! predation efficiency at noon
l1 = 1.0					#predation efficiency in perfect dark
l_shape = 5.0			#determines the sensitivity for small light
#if l_shape is small, predation efficiency under twilight is mostly same with that in perfect dark.
#if l_shape is large, predation efficiency under twilight is mostly same with that at noon

#mortality rate of prey by predation
mx = 1.0	#predation by other predators
my = 0.5 #predation by sharks




# Fig4a  -------------------------
v0 = 1.5  #average swim speed (prey is always 1.0)
my = 0.5  #predation by sharks

sharktemp=calc.sharktemp(t,t_w,wmax,wmin,r) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)
L = calc.light_effect(t,l0,l_shape,0.3)

plot_and_save.sim_result_with_wave("Fig4a",V,U,L,alpha,beta,mx,my,mb,phi,omega,h)


# Fig4b  -------------------------
v0 = 1.5  #average swim speed (prey is always 1.0)
my = 1.0  #predation by sharks

sharktemp=calc.sharktemp(t,t_w,wmax,wmin,r) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)
L = calc.light_effect(t,l0,l_shape,0.3)

plot_and_save.sim_result_with_wave("Fig4b",V,U,L,alpha,beta,mx,my,mb,phi,omega,h)


# Fig4c  -------------------------
v0 = 1.5  #average swim speed (prey is always 1.0)
my = 0.0  #predation by sharks

sharktemp=calc.sharktemp(t,t_w,wmax,wmin,r)  
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)
L = calc.light_effect(t,l0,l_shape,0.3)

plot_and_save.sim_result_with_wave("Fig4c",V,U,L,alpha,beta,mx,my,mb,phi,omega,h)

# Fig4d  -------------------------
v0 = 1.75  #average swim speed (prey is always 1.0)
my = 0.5  #predation by sharks

sharktemp=calc.sharktemp(t,t_w,wmax,wmin,r) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)
L = calc.light_effect(t,l0,l_shape,0.3)

plot_and_save.sim_result_with_wave("Fig4d",V,U,L,alpha,beta,mx,my,mb,phi,omega,h)



