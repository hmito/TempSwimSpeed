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

# Fig C1 --------------------
vk= 0.2
uk = 0.2

r = seq(0,1,length=101)
png("FigC1a.png", height = 1200,width=1600)
par(cex=4.0)
plot(r,delay(r,vk,uk),type="l", lwd=6,ylim=c(0,12),xlab="",ylab="")
lines(r,delay(r,2.0*vk,uk),lwd=6,col="blue")
lines(r,delay(r,1.25*vk,uk),lwd=6,col="forestgreen")
lines(r,delay(r,0.8*vk,uk),lwd=6,col="orange")
lines(r,delay(r,0.5*vk,uk),lwd=6,col="red")
dev.off()

png("FigC1b.png", height = 1200,width=1600)
par(cex=4.0)
r = seq(0,1,length=101)
plot(r,ampl(r,1.0*vk,uk),type="l", lwd=6,ylim=c(0,2),xlab="",ylab="")
lines(r,ampl(r,3.0*vk,uk),lwd=6,col="purple")
lines(r,ampl(r,2.0*vk,uk),lwd=6,col="blue")
lines(r,ampl(r,1.25*vk,uk),lwd=6,col="forestgreen")
lines(r,ampl(r,0.8*vk,uk),lwd=6,col="orange")
lines(r,ampl(r,0.5*vk,uk),lwd=6,col="red")
dev.off()

t = seq(-24,48,by=0.1)
tw = 15
wmax=30
wmin=25
v0 = 1.5
vk = 0.2
h = 1
c = 0.15

u0 = 1.0
uk =0.2
a = 1.0
omega = 1.0
mb = 0.2
mx = 1.0
my = 0.5
b = 1
l0 = 0.5
lk = 5.0
ll = 0.3
phi = 0.1

swimdif = function(t){
  v0-u0+(wmax-wmin)/2*uk*ampl(r,vk,uk)*cos((t-tw-delay(r,vk,uk))/12*acos(-1))
}
light_effect = function(t){
  #Around twilight_coef = 0.3 seems to be reasonable because
  #	- the definition of astronominal twilight is that the sun is less than -18 degree below the horizon. 
  #	- When the Culmination altitute = 90 degree, 18 degree passes for 1.2 hours, so twilight start from t = 4.8 and finish at t=19.2.
  #	- This time becomes longer when the Culmination altitude is less than 90 degree.
  #	- At twilight_coef = 0.3,
  #		- the sea is perfectly dark when t is from 20.5 to 3.5.
  #		- the sea is slightly bright (3% of noon) at t=4.5, 19.5
  #		- the light level rach to 20% and 40% of noon at t=5.5, 6.5.
  twilight_coef = 0.3
  l_min = l0
  light_influence = lk
  lwave=(1-twilight_coef)*cos(2*pi*(t-12)/24) + twilight_coef
  lwave[lwave<0] = 0
  alpha = log(1.0/l_min)
  return(exp(-alpha*lwave^(1/light_influence)))
}



plot_daycycle = function(ylim){
  plot(0:24,0:24,xaxs = "i", ylim=ylim,xaxt="n",yaxt="n",type="n",xlab="",ylab="")
  axis(1,at=c(0,3,6,9,12,15,18,21,24))
  polygon(c(-100,-100,6,6),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
  polygon(c(18,18,100,100),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
  lines(c(0,0),c(-100,100))
  lines(c(12,12),c(-100,100))
  lines(c(24,24),c(-100,100))
}
plot_range_box = function(x,y,range,col="black"){
  if(all(range)){
    start = x[1]
    stop = x[length(x)]
  }else{
    start = x[c(range[-1],FALSE)&c(!range[-length(range)],FALSE)]
    stop = x[c(FALSE,!range[-1])&c(FALSE,range[-length(range)])]
  }
  if(length(start)>=1){
    if(start[1]>stop[1]){
      start = start[-length(start)]
      stop = stop[-1]
    }
    for(i in 1:length(start)){
      polygon(c(start[i],start[i],stop[i],stop[i]),y+c(0,0.6,0.6,0),col=col,border=col)
    }
  }
}

# Fig D1 -------------------------------------
vk = uk*2
l0=1.0
lty="solid"
r = 0.6
png("FigD1a.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot(t,swimdif(t),type="n",xlim=c(0,24),xaxs="i",xlab="",ylab="",xaxt="n")
axis(1,at=c(0,3,6,9,12,15,18,21,24))
polygon(c(-100,-100,6,6),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
polygon(c(18,18,100,100),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
lines(c(0,0),c(-100,100))
lines(c(12,12),c(-100,100))
lines(c(24,24),c(-100,100))
lines(c(-100,100),c(0.8,0.8),lty="dotted",lwd=4,col="black")
lines(c(-100,100),c(0.4,0.4),lty="dotted",lwd=4,col="black")
r = 0.6
lines(t,swimdif(t)*light_effect(t),lwd=6,col="blue",lty=lty)
r = 0.3
lines(t,swimdif(t)*light_effect(t),lwd=6,col="forestgreen",lty=lty)
r = 0.1
lines(t,swimdif(t)*light_effect(t),lwd=6,col="red",lty=lty)
dev.off()

vk = uk*0.5
l0=1.0
lty="solid"
r = 0.6
png("FigD1b.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot(t,swimdif(t),type="n",xlim=c(0,24),xaxs="i",xlab="",ylab="",xaxt="n")
axis(1,at=c(0,3,6,9,12,15,18,21,24))
polygon(c(-100,-100,6,6),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
polygon(c(18,18,100,100),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
lines(c(0,0),c(-100,100))
lines(c(12,12),c(-100,100))
lines(c(24,24),c(-100,100))
lines(c(-100,100),c(0.8,0.8),lty="dotted",lwd=4,col="black")
lines(c(-100,100),c(0.4,0.4),lty="dotted",lwd=4,col="black")

r = 0.6
lines(t,swimdif(t)*light_effect(t),lwd=6,col="blue",lty=lty)
r = 0.3
lines(t,swimdif(t)*light_effect(t),lwd=6,col="forestgreen",lty=lty)
r = 0.1
lines(t,swimdif(t)*light_effect(t),lwd=6,col="red",lty=lty)
dev.off()

png("FigD1c.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot_daycycle(c(-0.4,6))
x = seq(-24,48,by=0.1)
y=5
vk=uk*2
thr = 0.8
r = 0.1
l0=1.0
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="red")
y=4
r =0.3
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="forestgreen")
y =3
r =0.6
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="blue")
y=2
thr = 0.4
r = 0.1
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="red")
y=1
r =0.3
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="forestgreen")
y =0
r =0.6
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="blue")
dev.off()

png("FigD1d.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot_daycycle(c(-0.4,6))
x = seq(-24,48,by=0.1)
y=5
thr = 0.8
vk=uk*0.5
r=0.1
l0=1.0
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="red")
y=4
r=0.3
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="forestgreen")
y=3
r = 0.6
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="blue")
y=2
thr = 0.4
r = 0.1
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="red")
y=1
r=0.3
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="forestgreen")
y=0
r=0.6
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr,col="blue")
dev.off()
# Fig D2 -------------------------------------
kai = 3.0
phi = 0.1
mx = 1.0
omega = 1.9
ub = 1.0
uk = 0.2
ut_thr = (kai*(1-phi)*mx - 1)/omega
ut_range = acos((ut_thr-ub)/uk)/acos(-1)*12

vk = uk*2
v0=1.5
l0=1.0
lty="solid"
r=0.6
png("FigD2a.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot(t,swimdif(t),type="n",xlim=c(0,24),xaxs="i",xlab="",ylab="",xaxt="n")
polygon(c(-100,-100,6,6),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
polygon(c(18,18,100,100),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
polygon(c(15-ut_range,15-ut_range,15+ut_range,15+ut_range),c(-100,100,100,-100),col=rgb(1,0,0,0.2),border=rgb(0,0,0,0))
lines(c(0,0),c(-100,100))
lines(c(12,12),c(-100,100))
lines(c(24,24),c(-100,100))
lines(c(-100,100),c(0.2,0.2),lty="dotted",lwd=4,col="black")
lines(c(-100,100),c(0.5,0.5),lty="dotted",lwd=4,col="black")
r=0.6
lines(t,swimdif(t)*light_effect(t),lwd=6,col="blue",lty=lty)
r=0.3
lines(t,swimdif(t)*light_effect(t),lwd=6,col="forestgreen",lty=lty)
r=0.1
lines(t,swimdif(t)*light_effect(t),lwd=6,col="red",lty=lty)
dev.off()

png("FigD2c.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot_daycycle(c(-0.4,6))
x = seq(-24,48,by=0.1)
y=5
vk=uk*2
thr = 0.5
IsX = (15-ut_range<=x & x<=15+ut_range)
r=0.1
l0=1.0
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="red")
y=4
r=0.3
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="forestgreen")
y=3
r=0.6
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="blue")
y=2
thr = 0.2
r=0.1
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="red")
y=1
r=0.3
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="forestgreen")
y=0
r=0.6
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="blue")
dev.off()

kai = 3.0
phi = 0.0
mx = 1.0
omega = 1.9
ub = 1.0
uk = 0.2
ut_thr = (kai*(1-phi)*mx - 1)/omega
ut_range = acos((ut_thr-ub)/uk)/acos(-1)*12

vk = uk*2
v0=1.5
l0=1.0
lty="solid"
r=0.6
png("FigD2b.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot(t,swimdif(t),type="n",xlim=c(0,24),xaxs="i",xlab="",ylab="",xaxt="n")
polygon(c(-100,-100,6,6),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
polygon(c(18,18,100,100),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
polygon(c(15-ut_range,15-ut_range,15+ut_range,15+ut_range),c(-100,100,100,-100),col=rgb(1,0,0,0.2),border=rgb(0,0,0,0))
lines(c(0,0),c(-100,100))
lines(c(12,12),c(-100,100))
lines(c(24,24),c(-100,100))
lines(c(-100,100),c(0.2,0.2),lty="dotted",lwd=4,col="black")
lines(c(-100,100),c(0.5,0.5),lty="dotted",lwd=4,col="black")
r=0.6
lines(t,swimdif(t)*light_effect(t),lwd=6,col="blue",lty=lty)
r=0.3
lines(t,swimdif(t)*light_effect(t),lwd=6,col="forestgreen",lty=lty)
r=0.1
lines(t,swimdif(t)*light_effect(t),lwd=6,col="red",lty=lty)
dev.off()

png("FigD2d.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot_daycycle(c(-0.4,6))
x = seq(-24,48,by=0.1)
y=5
vk=uk*2
thr = 0.5
IsX = (15-ut_range<=x & x<=15+ut_range)
r=0.1
l0=1.0
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="red")
y=4
r=0.3
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="forestgreen")
y=3
r=0.6
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="blue")
y=2
thr = 0.2
r=0.1
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="red")
y=1
r=0.3
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="forestgreen")
y=0
r=0.6
plot_range_box(x,y,swimdif(t)*light_effect(t)>thr&IsX,col="blue")
dev.off()





# Fig D3 -------------------------------------
kai = 0.9
phi = 0.1
mx = 1.0
my = 0.5
h = 1.0
pred_thr = (1-kai*(1-phi)*mx)/ (h*kai*(1-phi)*mx+kai*my-h)
vk = uk*2
l0=1.0
lty="solid"
v0=1.5
r=0.6
png("FigD3a.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot(t,swimdif(t),type="n",xlim=c(0,24),xaxs="i",xlab="",ylab="",xaxt="n")
axis(1,at=c(0,3,6,9,12,15,18,21,24))
polygon(c(-100,-100,6,6),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
polygon(c(18,18,100,100),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
lines(c(0,0),c(-100,100))
lines(c(12,12),c(-100,100))
lines(c(24,24),c(-100,100))
lines(c(-100,100),c(pred_thr,pred_thr),lty="dashed",lwd=4,col="purple")
lines(c(-100,100),c(0.4,0.4),lty="dotted",lwd=4,col="black")
lines(c(-100,100),c(0.1,0.1),lty="dotted",lwd=4,col="black")
r=0.6
lines(t,swimdif(t)*light_effect(t),lwd=6,col="blue",lty=lty)
r=0.3
lines(t,swimdif(t)*light_effect(t),lwd=6,col="forestgreen",lty=lty)
r=0.1
lines(t,swimdif(t)*light_effect(t),lwd=6,col="red",lty=lty)
dev.off()


png("FigD3c.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot_daycycle(c(-0.4,6))
x = seq(-24,48,by=0.1)
y=5
vk=uk*2
v0 = 1.5
thr = pred_thr
r=0.1
l0=1.0
plot_range_box(x,y,0.4<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="red")
y=4
r = 0.3
plot_range_box(x,y,0.4<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="forestgreen")
y=3
r = 0.6
plot_range_box(x,y,0.4<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="blue")
y=2
r=0.1
plot_range_box(x,y,0.1<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="red")
y=1
r = 0.3
plot_range_box(x,y,0.1<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="forestgreen")
y=0
r=0.6
plot_range_box(x,y,0.1<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="blue")
dev.off()


kai = 0.9
phi = 0.1
mx = 1.0
my = 1.0
h = 1.0
pred_thr = (1-kai*(1-phi)*mx)/ (h*kai*(1-phi)*mx+kai*my-h)
vk = uk*2
l0=1.0
lty="solid"
v0=1.5
r=0.6
png("FigD3b.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot(t,swimdif(t),type="n",xlim=c(0,24),xaxs="i",xlab="",ylab="",xaxt="n")
axis(1,at=c(0,3,6,9,12,15,18,21,24))
polygon(c(-100,-100,6,6),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
polygon(c(18,18,100,100),c(-100,100,100,-100),col="lightgrey",border=rgb(0,0,0,0))
lines(c(0,0),c(-100,100))
lines(c(12,12),c(-100,100))
lines(c(24,24),c(-100,100))
lines(c(-100,100),c(pred_thr,pred_thr),lty="dashed",lwd=4,col="purple")
lines(c(-100,100),c(0.4,0.4),lty="dotted",lwd=4,col="black")
lines(c(-100,100),c(0.1,0.1),lty="dotted",lwd=4,col="black")
r=0.6
lines(t,swimdif(t)*light_effect(t),lwd=6,col="blue",lty=lty)
r=0.3
lines(t,swimdif(t)*light_effect(t),lwd=6,col="forestgreen",lty=lty)
r=0.1
lines(t,swimdif(t)*light_effect(t),lwd=6,col="red",lty=lty)
dev.off()

png("FigD3d.png",width=1600,height=1200)
par(cex=4.0,bg=rgb(0,0,0,0))
plot_daycycle(c(-0.4,6))
x = seq(-24,48,by=0.1)
y=5
vk=uk*2
v0 = 1.5
thr = pred_thr
r=0.1
l0=1.0
plot_range_box(x,y,0.4<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="red")
y=4
r = 0.3
plot_range_box(x,y,0.4<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="forestgreen")
y=3
r = 0.6
plot_range_box(x,y,0.4<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="blue")
y=2
r=0.1
plot_range_box(x,y,0.1<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="red")
y=1
r = 0.3
plot_range_box(x,y,0.1<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="forestgreen")
y=0
r=0.6
plot_range_box(x,y,0.1<swimdif(t)*light_effect(t)&swimdif(t)*light_effect(t)<thr,col="blue")
dev.off()
