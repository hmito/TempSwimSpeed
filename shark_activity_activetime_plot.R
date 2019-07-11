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
mass = 5 #effective body size in the context of heat balance
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




v0 = 1.5  #average swim speed (prey is always 1.0)
my = 0.0  #predation by sharks
FigName = "Fig4c"


#start draw figure -------------------------
sharktemp=calc.bodytemp(watertemp,mass) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)
L = calc.light_effect(t,l0,l_shape,0.3)

Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, rep(predcost,length=length(V)), L, my, phi, omega, beta, h, mb,mx)
pred_eff = L*(V-U)^beta
pred_thr = predcost/(1-predcost*h)
pred_sthr = predcost/phi/(1-predcost*h)
prey_eff= alpha*(1+omega*U)/((1-phi)*mx+my*pred_eff/(1+h*pred_eff))
prey_peff= alpha/(1-phi) * alpha*(1+omega*U)/(mx+my*(1-h*predcost)*pred_eff^2/(1+h*pred_eff)/((1-h*predcost)*pred_eff-predcost))
prey_reff= prey_eff*(Ans$ThresholdPreyFreq>0.99)+prey_peff*(Ans$ThresholdPreyFreq<=0.99)
prey_thr = Ans$PreyW


pred_pfo = pred_eff-pred_thr
prey_epfo = prey_eff - prey_thr
prey_pfo = prey_reff-prey_thr

dt=(-1):26 - 0.5
pred_pfo = pred_pfo[c(23,24,1:24,1,2)]
prey_epfo= prey_epfo[c(23,24,1:24,1,2)]
prey_pfo= prey_pfo[c(23,24,1:24,1,2)]

png(paste(FigName,"1.png",sep=""),height=1200,width=1600)
par(cex=4.0,mex=1.0,bg=rgb(0,0,0,0))
plot(rep(dt,times=3),c(pred_pfo,prey_epfo,prey_pfo),type="n",col="red",xaxt="n",xlim=c(0,24),ylim=c(-0.5,1.2),lwd=3,xlab="",ylab="")
lines(c(-100,100),c(0,0))
lines(c(-100,100),c(pred_sthr,pred_sthr-pred_thr),col="red",lwd=3)
lines(dt,prey_epfo,type="l",col="skyblue",lwd=8,lty="dotted")
lines(dt,prey_pfo,type="l",col="blue",lwd=8)
lines(dt,pred_pfo,type="l",col="red",lwd=8)
pred_pfo[pred_pfo<0]=0
prey_pfo[prey_pfo<0]=0
polygon(c(dt,rev(dt)),c(pred_pfo,rep(0,length=length(dt))),col=rgb(1,0,0,0.3),border=rgb(0,0,0,0))
polygon(c(dt,rev(dt)),c(prey_pfo,rep(0,length=length(dt))),col=rgb(0,0,1,0.3),border=rgb(0,0,0,0))
dev.off()


#plot simulation results
png(paste(FigName,"2.png",sep=""),height=1200,width=1600)
par(cex=4.0,mex=1.0,bg=rgb(0,0,0,0))

Prey= Ans$Prey[c(23,24,1:24,1,2)]
Predator = Ans$Predator[c(23,24,1:24,1,2)]

plot(0,0,type="n",
	  xlab="",ylab="",
	  xlim=c(0,tnum),ylim=c(-0.02,1.02),xaxt="n")
lines(dt,Prey,col="blue",lwd=8,lty="dashed")
lines(dt,Predator*0.99,col="red",lwd=8)
axis(1,at=c(0,6,12,18,24))
dev.off()


