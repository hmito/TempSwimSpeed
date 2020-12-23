#load library
source("shark_activity_functions.R")
source("shark_activity_plotfunctions.R")

# constant parameters ============
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
mu = 1.0			#predation rate at midnight
rho = -0.5		#effect of light on predation rate at noon
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


# Figure 1 default parameter image =======================
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
points(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",pch=16)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",pch=17)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,4),col="orange",pch=18)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,4),col="orange",lty=1)
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
plot(t,calc.light_effect(t,1e-10, lk, ld),type="n",xlim=c(0,tnum),xlab="",ylab="",xaxt="n",ylim=c(0,1))
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
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()

# Fig3 ========================
plot_and_save.sim_result_with_wave("Fig3c", t, tw, wmin, wmax, ub, uk, 1.5, vk, mu, rho, kappa, sigma, alpha, omega, phi, mb, mx, 0.0, r, cost, beta, h)
plot_and_save.sim_result_with_wave("Fig3d", t, tw, wmin, wmax, ub, uk, 2.0, vk, mu, rho, kappa, sigma, alpha, omega, phi, mb, mx, 1.5, r, cost, beta, h)
plot_and_save.sim_result_with_wave("Fig3b", t, tw, wmin, wmax, ub, uk, 1.5, vk, mu, rho, kappa, sigma, alpha, omega, phi, mb, mx, 0.5, r, cost, beta, h)
plot_and_save.sim_result_with_wave("Fig3a", t, tw, wmin, wmax, ub, uk, 1.5, vk, mu, rho, kappa, sigma, alpha, omega, phi, mb, mx, 0.2, r, cost, beta, h)


plot.vb.my.figures("fig3xa", t, tw, 27.5, 27.5, ub, uk, vb, vk, 
						 1.0,0.0,kappa,sigma, alpha, 0.0, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("fig3xb", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,0.0,kappa,sigma, alpha, 0.0, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("fig3xc", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,0.0,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)


plot_and_save.sim_result_with_wave("Fig3ya", t, tw, 27.5, 27.5, ub, uk, 1.4, vk, 
											  1.0, 0.0, kappa, sigma, alpha, 0, phi, mb, mx, 1.0, r, predcost, beta, h)

plot_and_save.sim_result_with_wave("Fig3yb", t, tw, 27.5, 27.5, ub, uk, 1.0, vk, 
											  1.0, 0.0, kappa, sigma, alpha, 0, phi, mb, mx, 1.0, r, predcost, beta, h)

plot_and_save.sim_result_with_wave("Fig3yc", t, tw, wmin, wmax, ub, uk, 1.2, vk, 
											  1.0, 0.0, kappa, sigma, alpha, 0, phi, mb, mx, 0.0, r, cost, beta, h)

plot_and_save.sim_result_with_wave("Fig3yd", t, tw, wmin, wmax, ub, uk, 1.2, vk, 
											  1.0, 0.0, kappa, sigma, alpha, 0, phi, mb, mx, 2.0, r, cost, beta, h)

plot_and_save.sim_result_with_wave("Fig3ye", t, tw, wmin, wmax, ub, uk, 1.65, vk, 
											  1.0, 0.0, kappa, sigma, alpha, 0, phi, mb, mx, 2.0, r, cost, beta, h)

plot_and_save.sim_result_with_wave("Fig3yf", t, tw, wmin, wmax, ub, uk, 1.2, vk, 
											  1.0, 0.0, kappa, sigma, alpha, omega, phi, mb, mx, 0.5, r, cost, beta, h)

plot_and_save.sim_result_with_wave("Fig3yg", t, tw, wmin, wmax, ub, uk, 1.6, vk, 
											  1.0, 0.0, kappa, sigma, alpha, omega, phi, mb, mx, 1.0, r, cost, beta, h)

plot_and_save.sim_result_with_wave("Fig3yh", t, tw, wmin, wmax, ub, uk, 1.4, vk, 
											  1.0, 0.0, kappa, sigma, alpha, omega, phi, mb, mx, 0.3, r, cost, beta, h)

plot_and_save.sim_result_with_wave("Fig3yi", t, tw, wmin, wmax, ub, uk, 1.45, vk, 
											  1.0, 0.0, kappa, sigma, alpha, omega, phi, mb, mx, 0.50, r, predcost, beta, h)

# FIGURE 4 basic zoneplot  ========================
plot_legend = FALSE#TRUE
plot.vb.my.figures("fig4", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)


# FIGURE 5 phi-r plot  ========================
plot.r.phi.figures("fig5a",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.r.phi.figures("fig5b",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma,alpha, omega, phi, mb, 0.0, my, r, cost, beta, h,plot_legend = plot_legend)
plot.r.phi.figures("fig5c",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, 0.5, r, cost, beta, h,plot_legend = plot_legend)
plot.r.phi.figures("fig5d",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, 2.0, r, cost, beta, h,plot_legend = plot_legend)
plot.r.phi.figures("fig5e",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu,-0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.r.phi.figures("fig5f",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 mu, 0.0,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.r.phi.figures("fig5g",t, tw, wmin, wmax, ub, 0.4, vb, vk, 
						 mu, rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.r.phi.figures("fig5h",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.1, 0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.r.phi.figures("fig5h",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.5, 0.5,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)


# FIGURE A2 simple models  ========================
rho.seq = c(-0.9,-0.5,0.0,+0.5,+0.9)
mu.seq = c(1.0,1.0,1.0,0.5,0.1)
#alpha.seq = c(2.0,1.0,2.0,1.0)
omega.seq = c(0.0,0.0,1.0)
wmin.seq = c(27.5,25,25)
wmax.seq = c(27.5,30,30)
for(iy in 1:length(rho.seq)){
	for(ix in 1:length(omega.seq)){
		plot.vb.my.figures(sprintf("figA2[%d,%d]",ix,iy), t, tw, wmin.seq[ix], wmax.seq[ix], 
								 ub, uk, vb, vk, 
								 mu.seq[iy],rho.seq[iy], kappa,sigma, alpha, omega.seq[ix], 0, 
								 mb, mx, my, r, cost, beta, h,plot_legend=FALSE )
	}
}


# Figure A3  ========================
r.seq = c(0.3,1.0,4.0)
phi.seq = c(0.0,0.1,0.3)
for(ix in 1:length(r.seq)){
	for(iy in 1:length(phi.seq)){
		plot.vb.my.figures(sprintf("figA3[%d,%d]",ix,iy), t, tw, wmin, wmax, ub, uk, vb, vk, 
								 mu,rho,kappa,sigma,alpha, omega, phi.seq[iy], mb, mx, my, r.seq[ix], cost, beta, h,plot_legend = plot_legend)
	}
}


# Figure A4  ========================
mx.seq = c(0.5,1.0,2.0)
my.seq = c(2.0,1.0,0.5)
for(ix in 1:length(mx.seq)){
	for(iy in 1:length(my.seq)){
		plot.r.phi.figures(sprintf("figA4[%d,%d]",ix,iy),
								 t, tw, wmin, wmax, ub, uk, vb, vk, 
								 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx.seq[ix], my.seq[iy], r, cost, beta, h,plot_legend = plot_legend)
	}
}


# Figure A5  ========================
rho.seq = c(-0.9,-0.5,0.0)
my.seq = c(2.0,1.0,0.5)
for(ix in 1:length(rho.seq)){
	for(iy in 1:length(my.seq)){
		plot.r.phi.figures(sprintf("figA5[%d,%d]",ix,iy),
								 t, tw, wmin, wmax, ub, uk, vb, vk, 
								 mu,rho.seq[ix],kappa,sigma, alpha, omega, phi, mb, mx, my.seq[iy], r, cost, beta, h,plot_legend = plot_legend)
	}
}

# Figure A6  ========================
rho.seq = c(-0.9,-0.5,0.0)
my.seq = c(2.0,1.0,0.5)
for(ix in 1:length(rho.seq)){
	for(iy in 1:length(my.seq)){
		plot.r.phi.figures(sprintf("figA6[%d,%d]",ix,iy),
								 t, tw, wmin, wmax, ub, uk, vb, vk, 
								 mu,rho.seq[ix],kappa,sigma, alpha, 0.0, phi, mb, mx, my.seq[iy], r, cost, beta, h,plot_legend = plot_legend)
	}
}

# Figure A7  ========================
cost.seq = c(0.10,0.15,0.20)
h.seq = c(2.0,1.0,0.5)
for(ix in 1:length(cost.seq)){
	for(iy in 1:length(h.seq)){
		plot.r.phi.figures(sprintf("figA7[%d,%d]",ix,iy),
								 t, tw, wmin, wmax, ub, uk, vb, vk, 
								 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost.seq[ix], beta, h.seq[iy],plot_legend = plot_legend)
	}
}

# FIGURE A8 vk-uk effect in r-phi plot  ========================
plot_legend = FALSE#TRUE
vk.seq = c(0.1,0.2,0.4)
uk.seq = c(0.1,0.2,0.4)
for(ix in 1:length(vk.seq)){
	for(iy in 1:length(uk.seq)){
		plot.r.phi.figures(sprintf("figA8[%d,%d]",ix,iy),
								 t, tw, wmin, wmax, ub, uk.seq[ix], vb, vk.seq[iy], 
								 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
	}
}
plot_legend = FALSE#TRUE
vk.seq = c(0.1,0.2,0.4)
uk.seq = c(0.1,0.2,0.4)
for(ix in 1:length(vk.seq)){
	for(iy in 1:length(uk.seq)){
		plot.vb.my.figures(sprintf("figA8x[%d,%d]",ix,iy),
								 t, tw, wmin, wmax, ub, uk.seq[ix], vb, vk.seq[iy], 
								 mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
	}
}

# Figure A9  ========================


png("FigA9a.png",width=1200,height = 1200)
# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
plot(t,calc.light_effect(t,1e-10, lk, ld),type="n",xlim=c(0,tnum),xlab="",ylab="",xaxt="n",ylim=c(0,1))
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}
lines(t,calc.light_effect(t,0.1,0.9, kappa,sigma),col="red",lty=1)
points(t,calc.light_effect(t,0.1,0.9, kappa,sigma),col="red",pch=15)
lines(t,calc.light_effect(t,0.5,0.5, kappa,sigma),col="blue",lty=1)
points(t,calc.light_effect(t,0.5,0.5, kappa,sigma),col="blue",pch=16)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()


plot.vb.my.figures("figA9b", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.5,0.5,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("figA9c", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.1,0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.r.phi.figures("figA9d",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.5, 0.5,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot_and_save.sim_result_with_wave("FigXc", t, tw, wmin, wmax, ub, uk, 1.5, vk, mu, rho, kappa, sigma, alpha, omega, phi, mb, mx, 0.0, r, cost, beta, h)

plot.r.phi.figures("figA9e",t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.1, 0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)


# Figure A10 ========================
r.seq = c(0.3,1,4)
phi.seq = c(0.0,0.1,0.3)
for(ix in 1:length(r.seq)){
	for(iy in 1:length(phi.seq)){
		plot.vb.my.figures(sprintf("figA10[%d,%d]",ix,iy), t, tw, wmin, wmax, ub, uk, vb, vk, 
								 0.5,0.5,kappa,sigma,alpha, omega, phi.seq[iy], mb, mx, my, r.seq[ix], cost, beta, h,plot_legend = plot_legend)
	}
}


# Figure A11 ========================
mx.seq = c(0.5,1.0,2.0)
my.seq = c(2.0,1.0,0.5)
for(ix in 1:length(mx.seq)){
	for(iy in 1:length(my.seq)){
		plot.r.phi.figures(sprintf("figA11[%d,%d]",ix,iy),
								 t, tw, wmin, wmax, ub, uk, vb, vk, 
								 0.5,0.5,kappa,sigma, alpha, omega, phi, mb, mx.seq[ix], my.seq[iy], r, cost, beta, h,plot_legend = plot_legend)
	}
}


# Figure A12 ocean ========================
r.seq = c(0.3,1,4)
phi.seq = c(0.0,0.1,0.3)
for(ix in 1:length(r.seq)){
	for(iy in 1:length(phi.seq)){
		plot.vb.my.figures(sprintf("figA12[%d,%d]",ix,iy), t, tw, 27, 28, ub, uk, vb, vk, 
								 0.5,0.5,kappa,sigma,alpha, omega, phi.seq[iy], mb, mx, my, r.seq[ix], cost, beta, h,plot_legend = plot_legend)
	}
}



# light effect with line ========================
plot_and_save.sim_result_with_wave("FigXa", t, tw, wmin, wmax, ub, 0, 1.3, 0, 1, -0.5, kappa, sigma, alpha, omega, phi, mb, 0, 0.5, r, cost, beta, h)
plot_and_save.sim_result_with_wave("FigXc", t, tw, wmin, wmax, ub, 0, 1.3, 0, 0.75, 0, kappa, sigma, alpha, omega, phi, mb, 0, 0.5, r, cost, beta, h)
plot_and_save.sim_result_with_wave("FigXe", t, tw, wmin, wmax, ub, 0, 1.3, 0, 0.5, 0.5, kappa, sigma, alpha, omega, phi, mb, 0, 0.5, r, cost, beta, h)
plot_and_save.sim_result_with_wave("FigXb", t, tw, wmin, wmax, ub, 0, 2, 0, 1, -0.5, kappa, sigma, alpha, omega, phi, mb, 0, 2, r, cost, beta, h)
plot_and_save.sim_result_with_wave("FigXd", t, tw, wmin, wmax, ub, 0, 2, 0, 0.75, 0, kappa, sigma, alpha, omega, phi, mb, 0, 2, r, cost, beta, h)
plot_and_save.sim_result_with_wave("FigXf", t, tw, wmin, wmax, ub, 0, 2, 0, 0.5, 0.5, kappa, sigma, alpha, omega, phi, mb, 0, 2, r, cost, beta, h)


# FIGURE 4 basic zoneplot  ========================
plot_legend = FALSE#TRUE
plot.vb.my.figures("FigYa", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1,-0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("FigYb", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1,-0.5,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("FigYc", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1,0,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("FigYd", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.5,0.5,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("FigYe", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 0.1,0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("FigZa", t, tw, wmin+2, wmax-2, ub, uk, vb, vk, 
						 1,-0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("FigZb", t, tw, wmin+2, wmax-2, ub, uk, vb, vk, 
						 1,-0.5,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("FigZc", t, tw, wmin+2, wmax-2, ub, uk, vb, vk, 
						 1,0,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("FigZd", t, tw, wmin+2, wmax-2, ub, uk, vb, vk, 
						 0.5,0.5,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
plot.vb.my.figures("FigZe", t, tw, wmin+2, wmax-2, ub, uk, vb, vk, 
						 0.1,0.9,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)
