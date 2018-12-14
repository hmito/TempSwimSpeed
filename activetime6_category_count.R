#load library
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
watertemp = wmin+(wmax-wmin)*(cos(2*pi*(t-t_w)/length(t))+1)/2

# amount of food availability for prey
alpha = rep(1.0, length=tnum)
# baseline mortality for prey (should pay both for resting and foraging)
mb = 0.1	
# average prey swim speed
u0 = 1.0
# metabolic cost for predators when they go out for predation
C = rep(0.1, length=tnum)


#=== default parameter values ===
omega = 1.0	#foraging efficiency increment by speed 
beta = 1.0 	#predation efficiency
phi = 0.1  	#probability of failing to hide in safe place
h = 1.0    	#handling time

#prey and predator speed
mass = 10 #effective bpdy size in the context of heat balance
uk = 0.2 #influence of bodytemp
v0 = 1.4  #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp

#light effects
l0 = 0.2					#REQUIRE non-zero value! predation efficiency at noon
l1 = 1.0					#predation efficiency in perfect dark
l_shape = 1.0			#determines the sensitivity for small light
#if l_shape is small, predation efficiency under twilight is mostly same with that in perfect dark.
#if l_shape is large, predation efficiency under twilight is mostly same with that at noon

#mortality rate of prey by predation
mx = 0.5	#predation by other predators
my = 0.5 #predation by sharks


grid = 51
x.ax = seq(0.0,0.4, length=grid)
y.ax = seq(1.0,3.0,length=grid)
z.ax = seq(1,30.0, length=grid)

#plot categories of simulation results with changing v0 and phi
nolist = numeric(64)
for(omega in c(0,1,2)){
for(beta in c(0.5,1,2)){
for(l0 in c(0.1,0.2,0.4)){
print(sprintf("%.1f,%.1f,%.1f",omega,beta,l0))
for(z in 1:length(z.ax)){
mass = z.ax[z]
sharktemp=calc.bodytemp(watertemp,mass) 
for(y in 1:length(y.ax)){
	v0 = y.ax[y]
	for(x in 1:length(x.ax)){
		phi = x.ax[x]
		
		#reclaculation of U,V,L (not necessary for U and L)
		U = u0 + uk*(watertemp-(wmax+wmin)/2)
		V = v0 + vk*(sharktemp-(wmax+wmin)/2)
		L = calc.light_effect(t,l0,l1,l_shape)
		
		#run simulation
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		#calculate category
		no = activetime6.get_category(Ans)+1
		nolist[no] = nolist[no]+1
	}
}
}
}
}
}

for(mx in c(0.0,0.5,1.0)){		
for(my in c(0.0,0.5,1.0)){		
print(sprintf("%.1f,%.1f",mx,my))
for(z in 1:length(z.ax)){
	mass = z.ax[z]
	sharktemp=calc.bodytemp(watertemp,mass) 
	for(y in 1:length(y.ax)){
		v0 = y.ax[y]
		for(x in 1:length(x.ax)){
			phi = x.ax[x]
			
			#reclaculation of U,V,L (not necessary for U and L)
			U = u0 + uk*(watertemp-(wmax+wmin)/2)
			V = v0 + vk*(sharktemp-(wmax+wmin)/2)
			L = calc.light_effect(t,l0,l1,l_shape)
			
			#run simulation
			Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
			#calculate category
			no = activetime6.get_category(Ans)+1
			nolist[no] = nolist[no]+1
		}
	}
}
}
}

write.csv(nolist,"nolist3.csv")
