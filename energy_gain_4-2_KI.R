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
mx = 1.0
my = 1.0

sharktemp=calc.bodytemp(watertemp,mass) 
U = u0 + uk*(watertemp-(tmax+tmin)/2)
V = v0 + vk*(sharktemp-(tmax+tmin)/2)
L = calc.light_effect(t,l0,l1)

plot.assumption(t, watertemp, sharktemp, V, U, L)

x.seq = seq(0.0,0.4,length=5)
y.seq = seq(1.0,2.0,length=5)

grid = 51
x.ax = seq(0.0,0.4, length=grid)
y.ax = seq(1.0,2.0,length=grid)

name = "mx10my10"

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
		no[x,y] = allpeak_no.sim_result(Ans,0.0)
		#	plot.sim_result(Ans,"")
	}
}

mode = allpeak_no.mode(no)
png(paste(name,"_image.png",sep=""),height=1600,width=1600)
par(mfrow=c(1,1),cex=4.0,bg=rgb(0,0,0,0))
allpeak_no.image(x.ax,y.ax,mode,xlab=bquote(varphi),ylab=bquote('v'['0']))
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)
dev.off()

#=== plot with changing mx and my
omega = 1.0
beta = 1.0
phi = 0.1
h = 1.0

# work out the shark's body temperature
#    mass is kind of the size of shark (like sharkradius)
#    I used this just for simplifying parameters
mass = 10 
sharktemp=calc.bodytemp(watertemp,mass) 
#bodytemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)

# prey immediately track temperature
utemp = 0.2 #influence of bodytemp
U = uave + utemp*(watertemp-(tmax+tmin)/2)

# sharks track their own temperature
vave = 1.4  #average swim speed (prey is always 1.0)
vtemp = 0.2 #influence of bodytemp
V = vave + vtemp*(sharktemp-(tmax+tmin)/2)

# calc light effect, or predation efficiency
#    e.g., in muddy (inclear) water this value will increase?
lmin = 1.0
lmax = 1.0
L = calc.light_effect(t,lmin,lmax)

plot.assumption(t, watertemp, sharktemp, V, U, L)

x.seq = seq(0.0,1.0,length=5)
y.seq = seq(0.0,1.0,length=5)

grid = 51
x.ax = seq(0,1,length=grid)
y.ax = seq(0,1,length=grid)

name = "case1"

png(paste(name,"_map.png",sep=""),height=1600,width=1600)
par(mfrow=c(length(x.seq),length(y.seq)),cex=2.0,mex=0.3)
for(my in rev(y.seq)){
	for(mx in x.seq){
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		plot.sim_result(Ans,bquote(list('m'['x']==.(mx),'m'['y']==.(my))))
	}
}
dev.off()

no = matrix(0,grid,grid)
for(y in 1:length(y.ax)){
	for(x in 1:length(x.ax)){
		mx = x.ax[x]
		my = y.ax[y]
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		no[x,y] = allpeak_no.sim_result(Ans,0.0)
		#	plot.sim_result(Ans,"")
	}
}

mode = allpeak_no.mode(no)
png(paste(name,"_image.png",sep=""),height=1600,width=1600)
par(mfrow=c(1,1),cex=4.0,bg=rgb(0,0,0,0))
allpeak_no.image(x.ax,y.ax,mode,xlab=bquote('m'['x']),ylab=bquote('m'['y']))
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)
dev.off()

#=== plot with changing mass (sharkradius) and phi ===
omega = 0.5
beta = 1.0
phi = 0.0
h = 0.1

mx = 0.5
my = 0.5
vave=1.2

# calc light effect, or predation efficiency
#    e.g., in muddy (inclear) water this value will increase?
lmin = 1.0
lmax = 1.0
L = calc.light_effect(t,lmin,lmax)

par(mfrow=c(6,4),mex=0.3)
#plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
for(mass in c(1,4,16,64)){
	# work out the shark's body temperature
	sharktemp=calc.bodytemp(watertemp,mass) #just for simplifying parameters
	
	plot(t,watertemp,col="blue",pch=15, xlab="",ylab="")
	lines(t,watertemp,col="blue",lty=1)
	points(t,sharktemp,col="red",pch=16)
	lines(t,sharktemp,col="red",lty=1)
}
for(phi in seq(0.0,1.0,length=5)){
	for(mass in c(1,4,16,64)){
		# work out the shark's body temperature
		sharktemp=calc.bodytemp(watertemp,mass) #just for simplifying parameters
		V = vave + vtemp*(sharktemp-(tmax+tmin)/2)
		
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		plot.sim_result(Ans,bquote(list('mass'==.(mass),varphi==.(phi))))
	}
}




#=== plot with changing mass (sharkradius) and vave (average shark swim speed) ===
omega = 0.1
beta = 1
phi = 0.0
h = 0.1

mx = 0.4
my = 0.4
par(mfrow=c(6,5),mex=0.3)
#plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
for(mass in c(1,4,16,64)){
	# work out the shark's body temperature
	sharktemp=calc.bodytemp(watertemp,mass) #just for simplifying parameters
	
	plot(t,watertemp,col="blue",pch=15, xlab="",ylab="")
	lines(t,watertemp,col="blue",lty=1)
	points(t,sharktemp,col="red",pch=16)
	lines(t,sharktemp,col="red",lty=1)
}
for(vave in seq(0.4,length=4,by=0.4)){
	for(mass in c(1,4,16,64)){
		# work out the shark's body temperature
		sharktemp=calc.bodytemp(watertemp,mass) #just for simplifying parameters
		V = vave + vtemp*(sharktemp-(tmax+tmin)/2)
		
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		plot.sim_result(Ans,bquote(list('mass'==.(mass),'v'['ave']==.(vave))))
	}
}




#=== plot with changing mass (sharkradius) and phi ===
omega = 0.0
beta = 0.25
phi = 0.1
h = 0.1

# work out the shark's body temperature
#    mass is kind of the size of shark (like sharkradius)
#    I used this just for simplifying parameters
mass = 10
sharktemp=calc.bodytemp(watertemp,mass) 
#bodytemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)

# prey immediately track temperature
utemp = 0.1 #influence of bodytemp
U = uave + utemp*(watertemp-(tmax+tmin)/2)

# sharks track their own temperature
vave = 1.2 #average swim speed (prey is always 1.0)
vtemp = 0.1 #influence of bodytemp
V = vave + vtemp*(sharktemp-(tmax+tmin)/2)

# calc light effect, or predation efficiency
#    e.g., in muddy (inclear) water this value will increase?
lmin = 1
lmax = 1
L = calc.light_effect(t,lmin,lmax)

#plot.assumption(t,watertemp,sharktemp, V,U,L)

grid = 51
no = matrix(0,grid,grid)
x.ax = seq(1,20,length=grid)
y.ax = seq(0,0.4,length=grid)

mx = 0.5
my = 0.5
for(y in 1:length(y.ax)){
	for(x in 1:length(x.ax)){
		mass = x.ax[x]
		phi = y.ax[y]
		# work out the shark's body temperature
		sharktemp=calc.bodytemp(watertemp,mass) #just for simplifying parameters
		V = vave + vtemp*(sharktemp-(tmax+tmin)/2)
		
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		no[x,y] = allpeak_no.sim_result(Ans,0.0)
	#	plot.sim_result(Ans,"")
	}
}

mode = allpeak_no.mode(no)

allpeak_no.image(x.ax,y.ax,mode,xlab="mass",ylab="phi")
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)
#image(x.ax,y.ax,mode,zlim=c(0,15),col=mode.colors())



#=== plot with changing mass (sharkradius) and phi ===
omega = 1.0
beta = 1.0
phi = 0.1
h = 0.1

# work out the shark's body temperature
#    mass is kind of the size of shark (like sharkradius)
#    I used this just for simplifying parameters
mass = 10
sharktemp=calc.bodytemp(watertemp,mass) 
#bodytemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)

# prey immediately track temperature
utemp = 0.1 #influence of bodytemp
U = uave + utemp*(watertemp-(tmax+tmin)/2)

# sharks track their own temperature
vave = 1.2 #average swim speed (prey is always 1.0)
vtemp = 0.1 #influence of bodytemp
V = vave + vtemp*(sharktemp-(tmax+tmin)/2)

# calc light effect, or predation efficiency
#    e.g., in muddy (inclear) water this value will increase?
lmin = 1
lmax = 1
L = calc.light_effect(t,lmin,lmax)

#plot.assumption(t,watertemp,sharktemp, V,U,L)

grid = 51
no = matrix(0,grid,grid)
x.ax = seq(0,2,length=grid)
y.ax = seq(0,2,length=grid)

mx = 1.0
my = 1.0
for(y in 1:length(y.ax)){
	for(x in 1:length(x.ax)){
		mx = x.ax[x]
		my = y.ax[y]
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		no[x,y] = allpeak_no.sim_result(Ans,0.0)
		#	plot.sim_result(Ans,"")
	}
}

mode = allpeak_no.mode(no)

allpeak_no.image(x.ax,y.ax,mode,xlab="mx",ylab="my")
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)



#=== plot with changing mass (sharkradius) and phi ===
omega = 1.0
beta = 1.0
phi = 0.2
h = 1.0

# work out the shark's body temperature
#    mass is kind of the size of shark (like sharkradius)
#    I used this just for simplifying parameters
mass = 10
sharktemp=calc.bodytemp(watertemp,mass) 
#bodytemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)

# prey immediately track temperature
utemp = 0.2 #influence of bodytemp
U = uave + utemp*(watertemp-(tmax+tmin)/2)

# sharks track their own temperature
vave = 1.4 #average swim speed (prey is always 1.0)
vtemp = 0.2 #influence of bodytemp
V = vave + vtemp*(sharktemp-(tmax+tmin)/2)

# calc light effect, or predation efficiency
#    e.g., in muddy (inclear) water this value will increase?
lmin = 1
lmax = 1
L = calc.light_effect(t,lmin,lmax)

#plot.assumption(t,watertemp,sharktemp, V,U,L)

grid = 51
no = matrix(0,grid,grid)
x.ax = seq(0,2,length=grid)
y.ax = seq(0,2,length=grid)

mx = 0.5
my = 0.5
for(y in 1:length(y.ax)){
	for(x in 1:length(x.ax)){
		beta  = x.ax[x]
		omega = y.ax[y]
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		no[x,y] = allpeak_no.sim_result(Ans,0.0)
		#	plot.sim_result(Ans,"")
	}
}

mode = allpeak_no.mode(no)
allpeak_no.image(x.ax,y.ax,mode,xlab="beta",ylab="omega")
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)


unique(no[mode$main==0])

#=== plot with changing mass (sharkradius) and phi ===
omega = 1.0
beta = 0.0
phi = 0.2
h = 1.0

# work out the shark's body temperature
#    mass is kind of the size of shark (like sharkradius)
#    I used this just for simplifying parameters
mass = 10
sharktemp=calc.bodytemp(watertemp,mass) 
#bodytemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)

# prey immediately track temperature
uk = 0.2 #influence of bodytemp
U = u0 + uk*(watertemp-(tmax+tmin)/2)

# sharks track their own temperature
v0 = 1.4 #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp
V = v0 + vk*(sharktemp-(tmax+tmin)/2)

# calc light effect, or predation efficiency
#    e.g., in muddy (inclear) water this value will increase?
l0 = 1
l1 = 1
L = calc.light_effect(t,l0,l1)

mx = 0.5
my = 0.5

#plot.assumption(t,watertemp,sharktemp, V,U,L)

#=== for mass and v0
grid = 51
no = matrix(0,grid,grid)
x.ax = seq(1,20,length=grid)
y.ax = seq(1.0,2.0,length=grid)

for(y in 1:length(y.ax)){
	for(x in 1:length(x.ax)){
		mass  = x.ax[x]
		v0 = y.ax[y]
 
		sharktemp=calc.bodytemp(watertemp,mass) 
		V = v0 + vk*(sharktemp-(tmax+tmin)/2)
		
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		no[x,y] = allpeak_no.sim_result(Ans,0.0)
		#	plot.sim_result(Ans,"")
	}
}
par(mfrow=c(1,1))
mode = allpeak_no.mode(no)
allpeak_no.image(x.ax,y.ax,mode,xlab="mass",ylab="v0")
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)


#=== plot with changing mass (sharkradius) and phi ===
omega = 1.0
beta = 1.0
phi = 0.2
h = 1.0

# work out the shark's body temperature
#    mass is kind of the size of shark (like sharkradius)
#    I used this just for simplifying parameters
mass = 1
sharktemp=calc.bodytemp(watertemp,mass) 
#bodytemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)

# prey immediately track temperature
uk = 0.2 #influence of bodytemp
U = u0 + uk*(watertemp-(tmax+tmin)/2)

# sharks track their own temperature
v0 = 1.4 #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp
V = v0 + vk*(sharktemp-(tmax+tmin)/2)

# calc light effect, or predation efficiency
#    e.g., in muddy (inclear) water this value will increase?
l0 = 1
l1 = 1
L = calc.light_effect(t,l0,l1)

mx = 0.5
my = 0.5

grid = 51
no = matrix(0,grid,grid)
x.ax = seq(0,0.4,length=grid)
y.ax = seq(1.0,2.0,length=grid)

for(y in 1:length(y.ax)){
	for(x in 1:length(x.ax)){
		phi  = x.ax[x]
		v0 = y.ax[y]
		
		sharktemp=calc.bodytemp(watertemp,mass) 
		V = v0 + vk*(sharktemp-(tmax+tmin)/2)
		
		Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
		no[x,y] = allpeak_no.sim_result(Ans,0.0)
		#	plot.sim_result(Ans,"")
	}
}
par(mfrow=c(1,1))
mode = allpeak_no.mode(no)
allpeak_no.image(x.ax,y.ax,mode,xlab="phi",ylab="v0")
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)



#=== plot with changing mass (sharkradius) and phi ===
omega = 1.0
beta = 1.0
phi = 0.0
h = 1.0

# work out the shark's body temperature
#    mass is kind of the size of shark (like sharkradius)
#    I used this just for simplifying parameters
mass = 10
sharktemp=calc.bodytemp(watertemp,mass) 
#bodytemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)

# prey immediately track temperature
uk = 0.2 #influence of bodytemp
U = u0 + uk*(watertemp-(tmax+tmin)/2)

# sharks track their own temperature
v0 = 1.4 #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp
V = v0 + vk*(sharktemp-(tmax+tmin)/2)

# calc light effect, or predation efficiency
#    e.g., in muddy (inclear) water this value will increase?
l0 = 1
l1 = 1
L = calc.light_effect(t,l0,l1)

mx = 0.5
my = 0.5

#plot.assumption(t,watertemp,sharktemp, V,U,L)

grid = 101
active = matrix(0,length(watertemp),grid)
y.ax = seq(1.0,20.0,length=grid)

for(y in 1:length(y.ax)){
	mass = y.ax[y]
	sharktemp=calc.bodytemp(watertemp,mass) 
	V = v0 + vk*(sharktemp-(tmax+tmin)/2)
	
	Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
	active[,y] = Ans$Predator
	#	plot.sim_result(Ans,"")
}

image(t,y.ax,active)


par(mfrow=c(1,1))
mode = allpeak_no.mode(no)
allpeak_no.image(x.ax,y.ax,mode,xlab="mass",ylab="v0")
allpeak_no.point(x.ax,y.ax,mode,cex=0.5)

#===base color===
#black	: no use
#white	: all use
#cyan		: from afternoon  to beforenoon (rest around noon)
#blue		: from afternoon  to mid-night
#red		: from beforenoon to mid-night (rest early morning)
#darkgreen	: beforenoon
#darkred	: afternoon
#darkblue: midnight
#purple	: late beforenoon to early beforenoon (rest befornoon)
#green: midnight to beforenoon (rest befornoon)
#yellowgreen	: late afternoon to early afternoon (rest afternoon)
#yellow	: beforenoon to afternoon (daytime)
#===spot color===
#darkgreen-circle : sub peak in beforenoon
#green-circle : sub peak in between midnight and beforenoon

mode[no==14]=13 #orange : late midnight to early midnight (rest midnight)