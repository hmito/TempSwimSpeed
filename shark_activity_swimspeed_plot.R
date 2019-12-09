source("shark_activity_functions.R")

light_mode = FALSE

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
lm = 0.5		#predation efficiency at noon
lk = 2		#determines the sensitivity for small light
ld = 0.3		#duration of twilight

#if light_influence is small, predation efficiency under twilight is mostly same with that in perfect dark.
#if light_influence is large, predation efficiency under twilight is mostly same with that at noon

#mortality rate of prey by predation
mx = 1.0	#predation by other predators
my = 0.5 #predation by sharks

watertemp=calc.watertemp(t,t_w,wmin,wmax)

png("Fig1a.png",width=1200,height = 1200)
# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
#plot(t,watertemp,col="blue",pch=15, xlab="time (t)",ylab="temperature")
plot(t,watertemp,col="blue",pch=15,type ="n", xlab="",ylab="",xaxt = "n")
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}
points(t,watertemp,col="blue",pch=15, xlab="",ylab="",xaxt = "n")
lines(t,watertemp,col="blue",lty=1)
#	text(7,25.0,bquote('w'['t']))
points(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="red",pch=16)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="red",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,0.1),col="purple",pch=17)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.1),col="purple",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,0.6),col="orange",pch=18)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.6),col="orange",lty=1)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()



sharktemp=calc.sharktemp(t,t_w,wmin,wmax,r) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)

png("Fig1b.png",width=1200,height = 1200)
# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
#plot(t,V,type="n",xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),xlab="time (t)",ylab="burst speed")
plot(t,V,type="n",xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),xlab="",ylab="",xaxt="n")
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}
segments(-100,0,100,0)
points(t,V,col="red",pch=16)
lines(t,V,col="red",lty=1)
#	text(15,2.3,bquote('v'['t']))
points(t,U,col="blue",pch=15)
lines(t,U,col="blue",lty=1)
#	text(18,1.9,bquote('u'['t']))
points(t,V-U,col="purple",pch=17)
lines(t,V-U,col="purple",lty=1)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()

png("Fig1c.png",width=1200,height = 1200)
# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
plot(t,calc.light_effect(t,1e-10, lk, ld),type="n",xlim=c(0,tnum),xlab="",ylab="",xaxt="n")
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}
lines(t,calc.light_effect(t,0.1, lk, ld),col="red",lty=1)
points(t,calc.light_effect(t,0.1, lk, ld),col="red",pch=15)
lines(t,calc.light_effect(t,1.0, lk, ld),col="blue",lty=1)
points(t,calc.light_effect(t,1.0, lk, ld),col="blue",pch=16)
lines(t,calc.light_effect(t,0.5, lk, ld),col="black",lty=1)
points(t,calc.light_effect(t,0.5, lk, ld),col="black",pch=17)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()
