#load library
source("shark_activity_functions.R")
source("shark_activity_plotfunctions.R")

# constant parameters ============
#set time
pi = acos(-1)
tnum = 24 # time of day
t = 1:tnum-0.5
t_w = 15

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
mu = 1.0			#predation rate at midnight
rho = 0.0		#effect of light on predation rate at noon
kappa = 0.5		#determines the sensitivity for small light
sigma = 0.3		#duration of twilight

#mortality rate of prey by predation
alpha = 1.0 #amount of food for prey
omega = 1.0	#foraging efficiency increment by speed 
phi = 0.1  	#probability of failing to hide in safe place
mb = 0.2	#baseline mortality rate
mx = 1.0	#predation by other predators
my = 1.0 #predation by sharks

#=== default parameter values ===
r = 1.0  	#predator's body radius (meter)
cost=0.15	#predation cost
beta = 1.0 	#predation efficiency
h = 1.0    	#handling time

plot_legend = FALSE#TRUE
light_mode = FALSE

# fig 1 and fig 5a ======================================================================
png("fig1a.png",width=1200,height = 1200)
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
points(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",pch=16)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",pch=17)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,4),col="orange",pch=18)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,4),col="orange",lty=1)
axis(3,at=c(0,4,8,12,16,20,24))
par(new=TRUE)
plot(0,0,type="n",xlab="",ylab="",xlim=c(0,1),xaxt="n",yaxt="n")
axis(1)
dev.off()



sharktemp=calc.sharktemp(t,t_w,wmin,wmax,r) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)

png("fig1b.png",width=1200,height = 1200)
# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
#plot(t,V,type="n",xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),xlab="time (t)",ylab="burst speed")
plot(t,V,type="n",xlim=c(0,tnum),ylim=c(0,2),xlab="",ylab="",xaxt="n",yaxs ="i")
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
axis(3,at=c(0,4,8,12,16,20,24))
par(new=TRUE)
plot(0,0,type="n",xlab="",ylab="",xlim=c(0,1),xaxt="n",yaxt="n")
axis(1)
dev.off()

png("fig5a.png",width=1200,height = 1200)
# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
plot(t,calc.light_effect(t,mu,-0.9, kappa,sigma),type="n",xlim=c(0,tnum),xlab="",ylab="",xaxt="n",xaxs="i",ylim=c(0,1))
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}
plt = hmRLib::cud.colorset()
lines(t,calc.light_effect(t,mu,-0.9, kappa,sigma),col=plt$red,lty=1,lwd=2)
points(t,calc.light_effect(t,mu,-0.9, kappa,sigma),col=plt$red,pch=15)
lines(t,calc.light_effect(t,mu,-0.5, kappa,sigma),col=plt$orange,lty=1,lwd=2)
points(t,calc.light_effect(t,mu,-0.5, kappa,sigma),col=plt$orange,pch=16)
lines(t,calc.light_effect(t,0.1,0.9, kappa,sigma),col=plt$blue,lty=1,lwd=2)
points(t,calc.light_effect(t,0.1,0.9, kappa,sigma),col=plt$blue,pch=17)
lines(t,calc.light_effect(t,0.5,0.5, kappa,sigma),col=plt$green,lty=1,lwd=2)
points(t,calc.light_effect(t,0.5,0.5, kappa,sigma),col=plt$green,pch=18)
lines(t,calc.light_effect(t,mu,0, kappa,sigma),col="black",lty=1,lwd=2)
points(t,calc.light_effect(t,mu,0, kappa,sigma),col="black",pch=4)
axis(3,at=c(0,4,8,12,16,20,24))
par(new=TRUE)
plot(0,0,type="n",xlab="",ylab="",xlim=c(0,1),xaxt="n",yaxt="n")
axis(1)
dev.off()

# fig 3 ================================================================================================
plot.vb.my.figures("fig3a", t, tw, wmin, wmax, ub, 0, vb, 0, 
						 mu,rho,kappa,sigma, alpha, 0.0, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("fig3d", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, 0.0, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("fig3h", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)


plot_and_save.sim_result_with_wave("fig3b", t, tw, wmin, wmax, ub, 0, 1.2, 0, 
											  mu,rho, kappa, sigma, alpha, 0, phi, mb, mx, 1.0, r, cost, beta, h)

plot_and_save.sim_result_with_wave("fig3c", t, tw, wmin, wmax, ub, 0, 1.6, 0, 
											  mu,rho, kappa, sigma, alpha, 0, phi, mb, mx, 1.0, r, cost, beta, h)

plot_and_save.sim_result_with_wave("fig3e", t, tw, wmin, wmax, ub, uk, 1.2, vk, 
											  mu,rho, kappa, sigma, alpha, 0, phi, mb, mx, 0.0, r, cost, beta, h)

plot_and_save.sim_result_with_wave("fig3f", t, tw, wmin, wmax, ub, uk, 1.2, vk, 
											  mu,rho, kappa, sigma, alpha, 0, phi, mb, mx, 2.0, r, cost, beta, h)

plot_and_save.sim_result_with_wave("fig3g", t, tw, wmin, wmax, ub, uk, 1.65, vk, 
											  mu,rho, kappa, sigma, alpha, 0, phi, mb, mx, 2.0, r, cost, beta, h)

plot_and_save.sim_result_with_wave("fig3i", t, tw, wmin, wmax, ub, uk, 1.2, vk, 
											  mu,rho, kappa, sigma, alpha, omega, phi, mb, mx, 0.5, r, cost, beta, h)

plot_and_save.sim_result_with_wave("fig3j", t, tw, wmin, wmax, ub, uk, 1.6, vk, 
											  mu,rho, kappa, sigma, alpha, omega, phi, mb, mx, 1.2, r, cost, beta, h)

plot_and_save.sim_result_with_wave("fig3k", t, tw, wmin, wmax, ub, uk, 1.5, vk, 
											  mu,rho, kappa, sigma, alpha, omega, phi, mb, mx, 0.4, r, cost, beta, h)

# fig 4 ================================================================================================
plot.vb.my.figures("fig4a", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r=0.3, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4b", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r=4.0, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4c", t, tw, wmin, wmax, ub, uk=0.1, vb, vk=0.3, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4d", t, tw, wmin, wmax, ub, uk=0.3, vb, vk=0.1, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4e", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi=0.0, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4f", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi=0.2, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4g", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h=0.25,plot_legend = plot_legend)

plot.vb.my.figures("fig4h", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h=4,plot_legend = plot_legend)

# fig 5 ================================================================================================
plot.vb.my.figures("fig5b", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.5,0.5,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("fig5c", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.1,0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("fig5d", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,-0.5,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("fig5e", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,-0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)


# fig S1 ================================================================================================
coef_k = 6.36 * 1e-6*3600

delay = function(r,vk,uk){
	m = 2*acos(-1)*r^2/coef_k/24
	rv = vk/uk
	acos((rv - 1 - m^2)/sqrt(((rv-1)^2+m^2)*(1+m^2)))/acos(-1)*12
}
ampl = function(r,vk,uk){
	m = 2*acos(-1)*r^2/coef_k/24
	rv = vk/uk
	sqrt(((rv-1)^2+m^2)/(1+m^2))
}
vk= 0.2
uk = 0.2

r = seq(0,1,length=101)
png("figS1a.png", height = 1200,width=1600)
par(cex=4.0)
plot(r,delay(r,vk,uk),type="l", lwd=6,ylim=c(0,12),xlab="",ylab="")
lines(r,delay(r,2.0*vk,uk),lwd=6,col="blue")
lines(r,delay(r,1.25*vk,uk),lwd=6,col="forestgreen")
lines(r,delay(r,0.8*vk,uk),lwd=6,col="orange")
lines(r,delay(r,0.5*vk,uk),lwd=6,col="red")
dev.off()

png("figS1b.png", height = 1200,width=1600)
par(cex=4.0)
r = seq(0,1,length=101)
plot(r,ampl(r,1.0*vk,uk),type="l", lwd=6,ylim=c(0,2),xlab="",ylab="")
lines(r,ampl(r,3.0*vk,uk),lwd=6,col="purple")
lines(r,ampl(r,2.0*vk,uk),lwd=6,col="blue")
lines(r,ampl(r,1.25*vk,uk),lwd=6,col="forestgreen")
lines(r,ampl(r,0.8*vk,uk),lwd=6,col="orange")
lines(r,ampl(r,0.5*vk,uk),lwd=6,col="red")
dev.off()
# constant parameters (again) ============
#set time
pi = acos(-1)
tnum = 24 # time of day
t = 1:tnum-0.5
t_w = 15

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
mu = 1.0			#predation rate at midnight
rho = 0.0		#effect of light on predation rate at noon
kappa = 0.5		#determines the sensitivity for small light
sigma = 0.3		#duration of twilight

#mortality rate of prey by predation
alpha = 1.0 #amount of food for prey
omega = 1.0	#foraging efficiency increment by speed 
phi = 0.1  	#probability of failing to hide in safe place
mb = 0.2	#baseline mortality rate
mx = 1.0	#predation by other predators
my = 1.0 #predation by sharks

#=== default parameter values ===
r = 1.0  	#predator's body radius (meter)
cost=0.15	#predation cost
beta = 1.0 	#predation efficiency
h = 1.0    	#handling time

plot_legend = FALSE#TRUE
light_mode = FALSE
# fig S3  ========================
cost.seq = c(0.1,0.20)
for(ix in 1:length(cost.seq)){
	plot.vb.my.figures(sprintf("figS3[%d]",ix),
							 t, tw, wmin, wmax, ub, uk, vb, vk, 
							 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost.seq[ix], beta, h,plot_legend = plot_legend)
}

# fig S4  ========================
mx.seq = c(0.5,2.0)
for(ix in 1:length(mx.seq)){
	plot.vb.my.figures(sprintf("figS4[%d]",ix),
							 t, tw, wmin, wmax, ub, uk, vb, vk, 
							 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx.seq[ix], my, r, cost, beta, h,plot_legend = plot_legend)
}

# fig S5  ========================
vk.seq = c(0.1,0.2,0.3)
uk.seq = c(0.1,0.2,0.3)
for(ix in 1:length(vk.seq)){
	for(iy in 1:length(uk.seq)){
		plot.vb.my.figures(sprintf("figS5[%d,%d]",ix,iy),
								 t, tw, wmin, wmax, ub, uk.seq[ix], vb, vk.seq[iy], 
								 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
	}
}

# fig S6  ========================
plot_legend = FALSE#TRUE
rho.seq = c(-0.9,-0.5,0.0,+0.5,+0.9)
mu.seq = c(1.0,1.0,1.0,0.5,0.1)
#alpha.seq = c(2.0,1.0,2.0,1.0)
omega.seq = c(0.0,0.0,1.0)
wmin.seq = c(27.5,25,25)
wmax.seq = c(27.5,30,30)
for(iy in 1:length(rho.seq)){
	for(ix in 1:length(omega.seq)){
		plot.vb.my.figures(sprintf("figS6[%d,%d]",ix,iy), t, tw, wmin.seq[ix], wmax.seq[ix], 
								 ub, uk, vb, vk, 
								 mu.seq[iy],rho.seq[iy], kappa,sigma, alpha, omega.seq[ix], 0, 
								 mb, mx, my, r, cost, beta, h,plot_legend=plot_legend )
	}
}

# fig S7  ========================
plot_legend = FALSE#TRUE
rho.seq = c(-0.9,-0.5,0.0,+0.5,+0.9)
mu.seq = c(1.0,1.0,1.0,0.5,0.1)
kappa.seq = c(0.1,1.0,2.0)
for(iy in 1:length(rho.seq)){
	for(ix in 1:length(omega.seq)){
		plot.vb.my.figures(sprintf("figS7[%d,%d]",ix,iy), t, tw, wmin, wmax, 
								 ub, uk, vb, vk, 
								 mu.seq[iy],rho.seq[iy], kappa.seq[ix],sigma, alpha, omega, 0, 
								 mb, mx, my, r, cost, beta, h,plot_legend=plot_legend )
	}
}

# fig S8 ========================
r.seq = c(0.3,1,4)
phi.seq = c(0.0,0.1,0.3)
for(ix in 1:length(r.seq)){
	for(iy in 1:length(phi.seq)){
		plot.vb.my.figures(sprintf("figS8[%d,%d]",ix,iy), 
								 t, tw, 27, 28, ub, uk, vb, vk, 
								 0.5,0.5,kappa,sigma,alpha, omega, phi.seq[iy], mb, mx, my, r.seq[ix], cost, beta, h,plot_legend = plot_legend)
	}
}





