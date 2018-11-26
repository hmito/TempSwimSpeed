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


#=== draw phi-v0 figures === 
sharktemp=calc.bodytemp(watertemp,mass) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)
L = calc.light_effect(t,l0,l1,l_shape)
#draw assumptions (temperature and swim speed)
plot.assumption(t, watertemp, sharktemp, V, U, L)

x.seq = seq(0.0,0.25,length=7)
y.seq = seq(1.0,2.0,length=7)

grid = 101
x.ax = seq(0.0,0.4, length=grid)
y.ax = seq(1.0,3.0,length=grid)

name = "phi_v0_default"

#plot multiple results of simulations with changing v0 and phi
png(paste(name,"_map.png",sep=""),height=2000,width=2000)
par(mfrow=c(length(x.seq),length(y.seq)),cex=2.0,mex=0.3)
for(v0 in rev(y.seq)){
	for(phi in x.seq){
		#reclaculation of U,V,L (not necessary for U and L)
		U = u0 + uk*(watertemp-(wmax+wmin)/2)
		V = v0 + vk*(sharktemp-(wmax+wmin)/2)
		L = calc.light_effect(t,l0,l1,l_shape)

		#run simulation
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		#plot simulation results
		plot.sim_result(Ans,bquote(list(varphi==.(phi),'v'['0']==.(v0))))
	}
}
dev.off()

#plot categories of simulation results with changing v0 and phi
no = matrix(0,grid,grid)
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
		no[x,y] = majortime5.get_category(Ans,0.20)
	}
}
#reset to default values
phi = 0.1
v0 = 1.4
#plot category image
plotmode = majortime5.get_plotmode(no)
png(paste(name,"_image.png",sep=""),height=1600,width=1600)
par(mfrow=c(1,1),cex=4.0,bg=rgb(0,0,0,0))
image.plotmode(x.ax,y.ax,plotmode, xlab="phi",ylab="v0")
dev.off()

#list of categorization error (grey colors) 
plotmode$err_category

#list of categorization error (grey colors) by ignoring the number of peaks
sort(unique((plotmode$err_category)%%100))



#=== draw mass-v0 figures === 
sharktemp=calc.bodytemp(watertemp,mass) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)
L = calc.light_effect(t,l0,l1,l_shape)
plot.assumption(t, watertemp, sharktemp, V, U, L)

x.seq = seq(1,20.0,length=7)
y.seq = seq(1.0,2.0,length=7)

grid = 101
x.ax = seq(1,30.0, length=grid)
y.ax = seq(1.0,3.0,length=grid)

name = "mass_v0_default"

#plot multiple results of simulations with changing v0 and mass
png(paste(name,"_map.png",sep=""),height=2000,width=2000)
par(mfrow=c(length(x.seq),length(y.seq)),cex=2.0,mex=0.3)
for(v0 in rev(y.seq)){
	for(mass in x.seq){
		#reclaculation of U,V,L (not necessary for U and L)
		sharktemp=calc.bodytemp(watertemp,mass) 
		U = u0 + uk*(watertemp-(wmax+wmin)/2)
		V = v0 + vk*(sharktemp-(wmax+wmin)/2)
		L = calc.light_effect(t,l0,l1,l_shape)

		#run simulation
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		#plot simulation results
		plot.sim_result(Ans,bquote(list(varphi==.(phi),'v'['0']==.(v0))))
	}
}
dev.off()

#plot categories of simulation results with changing v0 and mass
no = matrix(0,grid,grid)
for(y in 1:length(y.ax)){
	v0 = y.ax[y]
	for(x in 1:length(x.ax)){
		mass = x.ax[x]
		
		#reclaculation of U,V,L (not necessary for U and L)
		sharktemp=calc.bodytemp(watertemp,mass) 
		U = u0 + uk*(watertemp-(wmax+wmin)/2)
		V = v0 + vk*(sharktemp-(wmax+wmin)/2)
		L = calc.light_effect(t,l0,l1,l_shape)
		
		#run simulation
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		#calculate category
		no[x,y] = majortime5.get_category(Ans,0.20)
	}
}
#reset to default values
mass = 10
v0 = 1.4
#plot category image
png(paste(name,"_image.png",sep=""),height=1600,width=1600)
par(mfrow=c(1,1),cex=4.0,bg=rgb(0,0,0,0))
image.plotmode(x.ax,y.ax,majortime5.get_plotmode(no),xlab="mass",ylab="v0")
dev.off()
