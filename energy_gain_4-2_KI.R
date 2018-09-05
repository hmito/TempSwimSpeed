#load library
library("Rcpp")
library("BH")
sourceCpp("TempSwimSpeedOptim.cpp")

#change of body size based on the mass of individual
#	mathematically same with calc.shark.temp function if we set
#  	mass = 1/(2*0.6*/(sharkradius*log(sharkradius/(sharkradius-skinthickness)))*0.031593/60 )
calc.bodytemp = function(watertemp, mass, error = 1e-10){
	bodytemp=rep(mean(watertemp),length(watertemp))
	for (i in c(1:1000)) {
		prev_bodytemp = bodytemp
		for (j in 1:length(watertemp)) {
			bodytemp[(j%%length(watertemp))+1] = bodytemp[j]+(watertemp[j]-bodytemp[j])/mass
		}
		
		#return if the difference from previous iteration is smaller than error
		if(sum(abs(bodytemp - prev_bodytemp))<error){
			return(bodytemp)
		}
	}
	#fail to calculate stable bodytemp
	return(rep(NA,length(watertemp)))
}

#light effect
calc.light_effect = function(t,lmin,lmax){
	lwave=exp(-1.0*cos(2*pi*(length(t)/2-t)/length(t)))/exp(1.0)
	return(lmin+(lmax-lmin)*lwave)
}

#Summarized figure
plot.sim_result = function(Ans,title){
	Prey= Ans$Prey
	Predator = Ans$Predator
	plot(0,0,type="n",
		  #		  xlab="time (t)",ylab="foraging predator (red), prey (blue)",
		  xlab="",ylab="",
		  xlim=c(0,tnum),ylim=c(-0.02,1.02),       
		  main=title
	)
	lines(t,Prey,col="blue",lwd=3)
	lines(t,Predator*0.99,col="red",lwd=3,lty="dashed")
}

#plot pair of temp and V,U,L
plot.assumption=function(t,watertemp,sharktemp,V,U,L){
	# assumption figure
	par(mfrow=c(1,2))
	plot(t,watertemp,col="blue",pch=15, xlab="time (t)",ylab="temperature")
	lines(t,watertemp,col="blue",lty=1)
	#	text(7,25.0,bquote('w'['t']))
	points(t,sharktemp,col="red",pch=16)
	lines(t,sharktemp,col="red",lty=1)
	#points(t,calc.bodytemp(watertemp,4),col="purple",pch=16)
	#lines(t,calc.bodytemp(watertemp,4),col="purple",lty=1)
	#points(t,calc.bodytemp(watertemp,25),col="orange",pch=16)
	#lines(t,calc.bodytemp(watertemp,25),col="orange",lty=1)
	#	text(4,28,bquote('s'['t']))
	par(new =T)
	plot(t,L,col="orange",pch=17, type = "p", axes = FALSE, ylab = "",ylim=c(0,max(L)),main="temperature")
	lines(t,L,col="orange",lty=1)
	#	text(5,0.25,bquote(lambda['t']))
	axis(4)
	
	plot(t,V,col="red",pch=16,xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),
		  xlab="time (t)",ylab="burst speed",main="swim speed")
	segments(-100,0,100,0)
	lines(t,V,col="red",lty=1)
	#	text(15,2.3,bquote('v'['t']))
	points(t,U,col="blue",pch=15)
	lines(t,U,col="blue",lty=1)
	#	text(18,1.9,bquote('u'['t']))
	points(t,V-U,col="purple",pch=2)
	lines(t,V-U,col="purple",lty=1)
	#	text(17,0.9,bquote(list('v'['t'],'- u'['t'])))
}

#=== constant parameters ===
#set time
pi = acos(-1)
tnum = 24 # time of day
t = 1:tnum-0.5

# temperature of the water
tempmaxt=15
tmin = 25
tmax = 30
watertemp = tmin+(tmax-tmin)*(cos(2*pi*(t-tempmaxt)/tnum)+1)/2

#amount of food availability for prey
alpha = rep(1.0, length=tnum)
# baseline mortality for prey (should pay both for resting and foraging)
mb = 0.1	
#average prey swim speed
uave = 1.0
#metabolic cost for predators when they go out for predation
C = rep(0.1, length=tnum)

#=== plot with changing mx (mortality by other predator) and my (mortality by shark)
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
vave = 1.2  #average swim speed (prey is always 1.0)
vtemp = 0.1 #influence of bodytemp
V = vave + vtemp*(sharktemp-(tmax+tmin)/2)

# calc light effect, or predation efficiency
#    e.g., in muddy (inclear) water this value will increase?
lmin = 1.0
lmax = 1.0
L = calc.light_effect(t,lmin,lmax)

plot.assumption(t, watertemp, sharktemp, V, U, L)

par(mfrow=c(5,5),mex=0.3)
for(my in seq(0.0,0.6,length=5)){
for(mx in seq(0.0,0.6,length=5)){
	Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
	plot.sim_result(Ans,bquote(list('m'['x']==.(mx),'m'['y']==.(my))))
}
}


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
par(mfrow=c(5,4),mex=0.3)
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
