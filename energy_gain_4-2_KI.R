#load library
library("Rcpp")
library("BH")
sourceCpp("TempSwimSpeedOptim.cpp")


#mathematically same with calc.shark.temp function if we set
#   mass = 1/(2*0.6*/(sharkradius*log(sharkradius/(sharkradius-skinthickness)))*0.031593/60 )
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

#Summarized figure
plot_sim_result = function(Ans,main){
	Prey= Ans$Prey
	Predator = Ans$Predator
	plot(0,0,type="n",
		  xlab="time (t)",ylab="foraging predator (red), prey (blue)",
		  xlim=c(0,tnum),ylim=c(-0.02,1.02),       
		  main=main
	)
	lines(t,Prey,col="blue",lwd=3)
	lines(t,Predator*0.99,col="red",lwd=3,lty="dashed")
}

pi = acos(-1)
tnum = 24 # time of day
t = 1:tnum-0.5

# influence of light on the predation rate
lmin=1.0 # lowest
lmax=1.0 # highest
lwave=exp(-1.0*cos(2*pi*(tnum/2-t)/tnum))/exp(1.0)
L = lmin+(lmax-lmin)*lwave

# temperature of the water
tempmaxt=15
tmin = 25
tmax = 30
watertemp = tmin+(tmax-tmin)*(cos(2*pi*(t-tempmaxt)/tnum)+1)/2

# work out the shark's body temperature
mass = 5
sharktemp=calc.bodytemp(watertemp,mass) #just forsimplify the parameters
#bodytemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)

# prey immediately track temperature
uave = 1.0
uamp = 0.5
U = uave + uamp*(watertemp-(tmax+tmin)/2)/(tmax-tmin)

# sharks track their own temperature
vave = 1.2
vamp = 1.0
V = vave + vamp*(sharktemp-(tmax+tmin)/2)/(tmax-tmin)

# assumption figure
par(mfrow=c(1,2))
plot(t,watertemp,col="blue",pch=15,,xlim=c(0,tnum),ylim=c(20,max(c(watertemp,sharktemp))),
	  xlab="time (t)",ylab="temperature")
lines(t,watertemp,col="blue",lty=1)
text(7,25.0,bquote('w'['t']))
points(t,sharktemp,col="red",pch=16)
lines(t,sharktemp,col="red",lty=1)
#points(t,calc.bodytemp(watertemp,4),col="purple",pch=16)
#lines(t,calc.bodytemp(watertemp,4),col="purple",lty=1)
#points(t,calc.bodytemp(watertemp,25),col="orange",pch=16)
#lines(t,calc.bodytemp(watertemp,25),col="orange",lty=1)
text(4,28,bquote('s'['t']))
par(new =T)
plot(t,L,col="orange",pch=17, type = "p", axes = FALSE, ylab = "",ylim=c(0,1.2))
lines(t,L,col="orange",lty=1)
text(5,0.25,bquote(lambda['t']))
axis(4)

plot(t,V,col="red",pch=16,xlim=c(0,tnum),ylim=c(0,max(c(V,U))),
	  xlab="time (t)",ylab="burst speed")
lines(t,V,col="red",lty=1)
text(15,2.3,bquote('v'['t']))
points(t,U,col="blue",pch=15)
lines(t,U,col="blue",lty=1)
text(18,1.9,bquote('u'['t']))
points(t,V-U,col="purple",pch=2)
lines(t,V-U,col="purple",lty=1)
text(17,0.9,bquote(list('v'['t'],'- u'['t'])))


#=== simulation parameters ===
#amount of food availability for prey
alpha = rep(1.0, length=tnum)
#metabolic cost for predators when they go out for predation
C = rep(0.05, length=tnum)	
#influence of swim speed on foraging efficiency for prey
omega = 0.1	#obtained reward is alpha*(1+omega*u)
#relative risk of predation for resting prey
phi = 0.1

#the predation rate: L*F*(v-u)^beta / {1 + h*L*F*(v-u)^beta}
#	L: influence of light on the predation rate
#	F: effective prey density, i.e., f + phi*(1-f)
beta = 1.0	#non-linear influence of speed difference on the predation rate 
h = 0.2	   #handling time for predation a prey

#following twp parameters determine the cost of prey
mb = 0.01	# baseline cost for prey (should pay both for resting and foraging)
mx = 0.05	# foraging cost for prey (should pay only for foraging)
my = 0.025 	# relative density of predator/prey

x11()
par(mfrow=c(7,7),mex=0.5)
for(mx in c(0.01,0.02,0.04,0.08,0.16,0.32,0.64)){
for(my in c(0.01,0.02,0.04,0.08,0.16,0.32,0.64)){
	Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C, L, my, phi, omega, beta, h, mb,mx)
	plot_sim_result(Ans,bquote(list(omega==.(omega),beta==.(beta),'m'['y']==.(my),'m'['x']==.(mx))))
}
}
