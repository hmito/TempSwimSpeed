#load library
source("shark_activity_functions.R")

#=== constant parameters ===
#set time
pi = acos(-1)
tnum = 24 # time of day
t = 1:tnum-0.5

# temperature of the water
tw=15
wmin = 25
wmax = 30

#prey and predator speed
ub = 1.0 #[can be fixed] average swim speed
uk = 0.2 #influence of bodytemp
vb = 1.5 #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp

#light effects
lm = 0.5		#predation efficiency at noon
lk = 0.2		#determines the sensitivity for small light
ld = 0.3		#duration of twilight

#mortality rate of prey by predation
alpha = 1.0 #amount of food for prey
omega = 1.0	#foraging efficiency increment by speed 
phi = 0.1  	#probability of failing to hide in safe place
mb = 0.2	#baseline mortality rate
mx = 1.0	#predation by other predators
my = 0.5 #predation by sharks

#=== default parameter values ===
r = 0.3  	#predator's body radius (meter)
cost=0.15	#predation cost
beta = 1.0 	#predation efficiency
h = 1.0    	#handling time

# Fig4a  -------------------------
vb = 1.5  #average swim speed (prey is always 1.0)
my = 0.0  #predation by sharks

plot_and_save.sim_result_with_wave("Fig4a", t, tw, wmin, wmax, ub, uk, vb, vk, lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h)

# Fig4b  -------------------------
vb = 1.5  #average swim speed (prey is always 1.0)
my = 1.25  #predation by sharks

plot_and_save.sim_result_with_wave("Fig4b", t, tw, wmin, wmax, ub, uk, vb, vk, lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h)

# Fig4c  -------------------------
vb = 1.5  #average swim speed (prey is always 1.0)
my = 0.5  #predation by sharks

plot_and_save.sim_result_with_wave("Fig4c", t, tw, wmin, wmax, ub, uk, vb, vk, lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h)

# Fig4d  -------------------------
vb = 1.7  #average swim speed (prey is always 1.0)
my = 1.0  #predation by sharks

plot_and_save.sim_result_with_wave("Fig4d", t, tw, wmin, wmax, ub, uk, vb, vk, lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h)

# Fig4e  -------------------------
vb = 1.8  #average swim speed (prey is always 1.0)
my = 1.5  #predation by sharks

plot_and_save.sim_result_with_wave("Fig4e", t, tw, wmin, wmax, ub, uk, vb, vk, lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h)
