#load library
source("shark_activity_functions.R")


#=== constant parameters ===
#set time
pi = acos(-1)
tnum = 24 # time of day
t = 1:tnum-0.5

# temperature of the water
t_w=15
wmax = 25
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


#=== plot with changing phi and v0
omega = 1.0	#foraging efficiency increment by speed 
beta = 1.0 	#predation efficiency
phi = 0.3  	#probability of failing to hide in safe place
h = 1.0    	#handling time

#prey and predator speed
mass = 10 #mass of shark
uk = 0.2 #influence of bodytemp
v0 = 1.4  #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp

#light effects
l0 = 1.0
l1 = 1.0

#mortality rate of prey by predation
mx = 0.5
my = 0.5

name = "l10"

sharktemp=calc.bodytemp(watertemp,mass) 
U = u0 + uk*(watertemp-(tmax+tmin)/2)
V = v0 + vk*(sharktemp-(tmax+tmin)/2)
L = calc.light_effect(t,l0,l1)

plot.assumption(t, watertemp, sharktemp, V, U, L)

x.seq = seq(0.0,0.4,length=5)
y.seq = seq(1.0,4.0,length=5)

grid = 51
x.ax = seq(0.0,0.4, length=grid)
y.ax = seq(1.0,4.0,length=grid)

png(paste(name,"_map.png",sep=""),height=1600,width=1600)
par(mfrow=c(length(x.seq),length(y.seq)),cex=2.0,mex=0.3)
for(v0 in rev(y.seq)){
	V = v0 + vk*(sharktemp-(tmax+tmin)/2)
	for(phi in x.seq){
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		plot.sim_result(Ans,bquote(list(varphi==.(phi),'v'['0']==.(v0))))
	}
}
dev.off()

no = matrix(0,grid,grid)
for(y in 1:length(y.ax)){
	v0 = y.ax[y]
	V = v0 + vk*(sharktemp-(tmax+tmin)/2)
	for(x in 1:length(x.ax)){
		phi = x.ax[x]
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		#calculate peak based category number
		no[x,y] = allpeak_no.sim_result(Ans,0.0)
		#	plot.sim_result(Ans,"")
	}
}
#transform categorized no to mode
mode = allpeak_no.mode(no)
png(paste(name,"_image.png",sep=""),height=1600,width=1600)
par(mfrow=c(1,1),cex=4.0,bg=rgb(0,0,0,0))
allpeak_no.image(x.ax,y.ax,mode,xlab=bquote(varphi),ylab=bquote('v'['0']))
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)
dev.off()




#=== plot with changing mass and v0
omega = 1.0	#foraging efficiency increment by speed 
beta = 1.0 	#predation efficiency
phi = 0.1  	#probability of failing to hide in safe place
h = 1.0    	#handling time

#prey and predator speed
mass = 10 
uk = 0.2 #influence of bodytemp
v0 = 1.4  #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp

#light effects
l0 = 1.0
l1 = 1.0

#mortality rate of prey by predation
mx = 0.5
my = 0.5

name = "mv0_default"

sharktemp=calc.bodytemp(watertemp,mass) 
U = u0 + uk*(watertemp-(tmax+tmin)/2)
V = v0 + vk*(sharktemp-(tmax+tmin)/2)
L = calc.light_effect(t,l0,l1)

#	plot.assumption(t, watertemp, sharktemp, V, U, L)

x.seq = seq(1, 20,length=5)
y.seq = seq(1.0, 2.0,length=5)

grid = 51
x.ax = seq(1, 20, length=grid)
y.ax = seq(1.0, 2.0,length=grid)

png(paste(name,"_map.png",sep=""),height=1600,width=1600)
par(mfrow=c(length(x.seq),length(y.seq)),cex=2.0,mex=0.3)
for(v0 in rev(y.seq)){
	for(mass in x.seq){
		sharktemp=calc.bodytemp(watertemp,mass) 
		V = v0 + vk*(sharktemp-(tmax+tmin)/2)
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		plot.sim_result(Ans,bquote(list('mass'==.(mass),'v'['0']==.(v0))))
	}
}
dev.off()

no = matrix(0,grid,grid)
for(y in 1:length(y.ax)){
	v0 = y.ax[y]
	for(x in 1:length(x.ax)){
		mass = x.ax[x]
		sharktemp=calc.bodytemp(watertemp,mass) 
		V = v0 + vk*(sharktemp-(tmax+tmin)/2)
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		no[x,y] = allpeak_no.sim_result(Ans,0.0)
		#	plot.sim_result(Ans,"")
	}
}

mode = allpeak_no.mode(no)
png(paste(name,"_image.png",sep=""),height=1600,width=1600)
par(mfrow=c(1,1),cex=4.0,bg=rgb(0,0,0,0))
allpeak_no.image(x.ax,y.ax,mode,xlab=bquote('mass'),ylab=bquote('v'['0']))
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)
dev.off()