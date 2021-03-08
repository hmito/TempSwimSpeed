plot.vb.my.figures("fig4a", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,0.0,kappa,sigma, alpha, omega, phi, mb, mx, my, r=0.3, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4b", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,0.0,kappa,sigma, alpha, omega, phi, mb, mx, my, r=4.0, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4c", t, tw, wmin, wmax, ub, uk=0.1, vb, vk=0.3, 
						 1.0,0.0,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4d", t, tw, wmin, wmax, ub, uk=0.3, vb, vk=0.1, 
						 1.0,0.0,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4e", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,0.0,kappa,sigma, alpha, omega, phi=0.0, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4f", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,0.0,kappa,sigma, alpha, omega, phi=0.2, mb, mx, my, r, cost, beta, h,plot_legend = plot_legend)

plot.vb.my.figures("fig4g", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,0.0,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h=0.25,plot_legend = plot_legend)

plot.vb.my.figures("fig4h", t, tw, wmin, wmax, ub, uk, vb, vk, 
						 1.0,0.0,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h=4,plot_legend = plot_legend)


lwd = 9
light_mode = TRUE

# Figure 1 default parameter image =======================
png("fig1aS.png",width=1200,height = 1200)
# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
#plot(t,watertemp,col="blue",pch=15, xlab="time (t)",ylab="temperature")
plot(t,watertemp,col="blue",pch=15,type ="n", xlab="",ylab="",xaxt = "n",xaxs="i")
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}
points(t,watertemp,col="blue",pch=15, xlab="",ylab="",xaxt = "n")
lines(t,watertemp,col="blue",lty=1,lwd=lwd)
#	text(7,25.0,bquote('w'['t']))
#points(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",pch=16)
#lines(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",lty=1)
#points(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",pch=17)
#lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="red",pch=16)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="red",lty=1,lwd=lwd)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()

png("fig1aM.png",width=1200,height = 1200)
# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
#plot(t,watertemp,col="blue",pch=15, xlab="time (t)",ylab="temperature")
plot(t,watertemp,col="blue",pch=15,type ="n", xlab="",ylab="",xaxt = "n",xaxs="i")
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}
points(t,watertemp,col="blue",pch=15, xlab="",ylab="",xaxt = "n")
lines(t,watertemp,col="blue",lty=1,lwd=lwd)
#	text(7,25.0,bquote('w'['t']))
#points(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",pch=16)
#lines(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",lty=1)
#points(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",pch=17)
#lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",pch=16)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",lty=1,lwd=lwd)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()

png("fig1aL.png",width=1200,height = 1200)
# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
#plot(t,watertemp,col="blue",pch=15, xlab="time (t)",ylab="temperature")
plot(t,watertemp,col="blue",pch=15,type ="n", xlab="",ylab="",xaxt = "n",xaxs="i")
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}
points(t,watertemp,col="blue",pch=15, xlab="",ylab="",xaxt = "n")
lines(t,watertemp,col="blue",lty=1,lwd=lwd)
#	text(7,25.0,bquote('w'['t']))
#points(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",pch=16)
#lines(t,calc.sharktemp(t,t_w,wmin,wmax,1),col="red",lty=1)
#points(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",pch=17)
#lines(t,calc.sharktemp(t,t_w,wmin,wmax,0.3),col="purple",lty=1)
points(t,calc.sharktemp(t,t_w,wmin,wmax,4),col="red",pch=16)
lines(t,calc.sharktemp(t,t_w,wmin,wmax,4),col="red",lty=1,lwd=lwd)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()



png("fig1bS.png",width=1200,height = 1200)
sharktemp=calc.sharktemp(t,t_w,wmin,wmax,0.3) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)

# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
#plot(t,V,type="n",xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),xlab="time (t)",ylab="burst speed")
plot(t,V,type="n",xlim=c(0,tnum),ylim=c(0,2),xlab="",ylab="",xaxt="n",xaxs="i")
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}

segments(-100,0,100,0)
points(t,V,col="red",pch=16)
lines(t,V,col="red",lty=1,lwd=lwd)
#	text(15,2.3,bquote('v'['t']))
points(t,U,col="blue",pch=15)
lines(t,U,col="blue",lty=1,lwd=lwd)
#	text(18,1.9,bquote('u'['t']))
points(t,V-U,col="purple",pch=17)
lines(t,V-U,col="purple",lty=1,lwd=lwd)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()



png("fig1bM.png",width=1200,height = 1200)
sharktemp=calc.sharktemp(t,t_w,wmin,wmax,1.0) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)

# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
#plot(t,V,type="n",xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),xlab="time (t)",ylab="burst speed")
plot(t,V,type="n",xlim=c(0,tnum),ylim=c(0,2),xlab="",ylab="",xaxt="n",xaxs="i")
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}

segments(-100,0,100,0)
points(t,V,col="red",pch=16)
lines(t,V,col="red",lty=1,lwd=lwd)
#	text(15,2.3,bquote('v'['t']))
points(t,U,col="blue",pch=15)
lines(t,U,col="blue",lty=1,lwd=lwd)
#	text(18,1.9,bquote('u'['t']))
points(t,V-U,col="purple",pch=17)
lines(t,V-U,col="purple",lty=1,lwd=lwd)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()



png("fig1bL.png",width=1200,height = 1200)
sharktemp=calc.sharktemp(t,t_w,wmin,wmax,4.0) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)

# assumption figure
par(cex=4.0,bg=rgb(0,0,0,0))
#plot(t,V,type="n",xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),xlab="time (t)",ylab="burst speed")
plot(t,V,type="n",xlim=c(0,tnum),ylim=c(0,2),xlab="",ylab="",xaxt="n",xaxs="i")
if(light_mode){
	polygon(c(-10,-10,100,100),c(-100,100,100,-100),col="white",border=rgb(0,0,0,0))
	polygon(c(-10,-10,4,4),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
	polygon(c(4,4,8,8),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(16,16,20,20),c(-100,100,100,-100),col="grey90",border=rgb(0,0,0,0))
	polygon(c(20,20,100,100),c(-100,100,100,-100),col="grey75",border=rgb(0,0,0,0))
}

segments(-100,0,100,0)
points(t,V,col="red",pch=16)
lines(t,V,col="red",lty=1,lwd=lwd)
#	text(15,2.3,bquote('v'['t']))
points(t,U,col="blue",pch=15)
lines(t,U,col="blue",lty=1,lwd=lwd)
#	text(18,1.9,bquote('u'['t']))
points(t,V-U,col="purple",pch=17)
lines(t,V-U,col="purple",lty=1,lwd=lwd)
axis(1,at=c(0,4,8,12,16,20,24))
dev.off()