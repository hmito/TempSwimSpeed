m = 0.5
rv = 10

delay = function(m,rv){
  acos((rv - 1 - m^2)/sqrt(((rv-1)^2+m^2)*(1+m^2)))/acos(-1)*12
}


M = seq(0,10,length=101)
m = 2*acos(-1)*M/24
png("Y_delay.png", height = 1200,width=1600)
par(cex=4.0)
plot(M,delay(m,1.0),type="l", lwd=6,ylim=c(0,12),xlab="",ylab="")
lines(M,delay(m,2.0),lwd=6,col="blue")
lines(M,delay(m,1.25),lwd=6,col="forestgreen")
lines(M,delay(m,0.8),lwd=6,col="orange")
lines(M,delay(m,0.5),lwd=6,col="red")
dev.off()

m = 0.5
rv = 10


ampl = function(m,rv){
  sqrt(((rv-1)^2+m^2)/(1+m^2))
}

png("Y_ampl.png", height = 1200,width=1600)
par(cex=4.0)
M = seq(0,10,length=101)
m = 2*acos(-1)*M/24
plot(M,ampl(m,1.0),type="l", lwd=6,ylim=c(0,2),xlab="",ylab="")
lines(M,ampl(m,3.0),lwd=6,col="purple")
lines(M,ampl(m,2.0),lwd=6,col="blue")
lines(M,ampl(m,1.25),lwd=6,col="forestgreen")
lines(M,ampl(m,0.8),lwd=6,col="orange")
lines(M,ampl(m,0.5),lwd=6,col="red")
dev.off()



M = seq(0,10,length=101)
m = 2*acos(-1)*M/24
rv = seq(0,3,length=101)
plot(M,ampl(0.2,rv),type="l", lwd=3,ylim=c(0,2),xlab="M",ylab="delay (hours)")
lines(M,ampl(0.2,rv),lwd=3,col="purple")
lines(M,ampl(0.5,rv),lwd=3,col="blue")
lines(M,ampl(1.0,rv),lwd=3,col="forestgreen")
lines(M,ampl(2.0,rv),lwd=3,col="orange")
#lines(M,ampl(m,0.5),lwd=3,col="red")


plot_daycycle = function(ylim){
  plot(0:24,0:24,ylim=ylim,xaxt="n",yaxt="n",type="n",xlab="",ylab="")
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

png("case_1.png",width=1600,height=900)
par(cex=3.0)
plot_daycycle(c(-0.4,6))
x = seq(-24,48,by=0.1)
y=5
thr = 0.4
m =0.1
rv=2
plot_range_box(x,y,(ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))>thr)
y=4
m =1.0
rv=2
plot_range_box(x,y,(ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))>thr)
y=3
m =10.0
rv=2
plot_range_box(x,y,(ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))>thr)
y=2
m =0.1
rv=0.2
plot_range_box(x,y,(ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))>thr)
y=1
m =1.0
rv=0.2
plot_range_box(x,y,(ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))>thr)
y=0
m =10.0
rv=0.2
plot_range_box(x,y,(ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))>thr)
dev.off()


png("case_2.png",width=1600,height=1200)
par(cex=3.0)
plot_daycycle(c(-0.4,9))
x = seq(-24,48,by=0.1)

thru = 10
thrl = -0.4
y=8
m =0.1
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrl<dv, col="skyblue")
plot_range_box(x,y,dv<thru, col="pink")
plot_range_box(x,y,thrl<dv&dv<thru)
y=5
m =1.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrl<dv, col="skyblue")
plot_range_box(x,y,dv<thru, col="pink")
plot_range_box(x,y,thrl<dv&dv<thru)
y=2
m =10.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrl<dv, col="skyblue")
plot_range_box(x,y,dv<thru, col="pink")
plot_range_box(x,y,thrl<dv&dv<thru)

thru = 0.4
thrl = -10
y=7
m =0.1
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrl<dv, col="skyblue")
plot_range_box(x,y,dv<thru, col="pink")
plot_range_box(x,y,thrl<dv&dv<thru)
y=4
m =1.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrl<dv, col="skyblue")
plot_range_box(x,y,dv<thru, col="pink")
plot_range_box(x,y,thrl<dv&dv<thru)
y=1
m =10.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrl<dv, col="skyblue")
plot_range_box(x,y,dv<thru, col="pink")
plot_range_box(x,y,thrl<dv&dv<thru)

thru = 0.4
thrl = -0.4
y=6
m =0.1
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrl<dv, col="skyblue")
plot_range_box(x,y,dv<thru, col="pink")
plot_range_box(x,y,thrl<dv&dv<thru)
y=3
m =1.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrl<dv, col="skyblue")
plot_range_box(x,y,dv<thru, col="pink")
plot_range_box(x,y,thrl<dv&dv<thru)
y=0
m =10.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrl<dv, col="skyblue")
plot_range_box(x,y,dv<thru, col="pink")
plot_range_box(x,y,thrl<dv&dv<thru)

dev.off()



png("case_3.png",width=1600,height=1200)
par(cex=3.0)

plot_daycycle(c(-0.4,9))
x = seq(-24,48,by=0.1)

thru = -0.5
thrdv = 0.4
y=8
m =0.1
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
u = cos((x-15)/24*2*acos(-1))
plot_range_box(x,y,thrdv<dv,col="skyblue")
plot_range_box(x,y,thru<u,col="pink")
plot_range_box(x,y,thrdv<dv&thru<u)
y=7
m =1.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrdv<dv,col="skyblue")
plot_range_box(x,y,thru<u,col="pink")
plot_range_box(x,y,thrdv<dv&thru<u)
y=6
m =10.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrdv<dv,col="skyblue")
plot_range_box(x,y,thru<u,col="pink")
plot_range_box(x,y,thrdv<dv&thru<u)

thru = 0.3
thrdv =-0.5
y=5
m =0.1
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
u = cos((x-15)/24*2*acos(-1))
plot_range_box(x,y,thrdv<dv,col="skyblue")
plot_range_box(x,y,thru<u,col="pink")
plot_range_box(x,y,thrdv<dv&thru<u)
y=4
m =1.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrdv<dv,col="skyblue")
plot_range_box(x,y,thru<u,col="pink")
plot_range_box(x,y,thrdv<dv&thru<u)
y=3
m =10.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrdv<dv,col="skyblue")
plot_range_box(x,y,thru<u,col="pink")
plot_range_box(x,y,thrdv<dv&thru<u)

thru = -0.5
thrdv =-0.5
y=2
m =0.1
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
u = cos((x-15)/24*2*acos(-1))
plot_range_box(x,y,thrdv<dv,col="skyblue")
plot_range_box(x,y,thru<u,col="pink")
plot_range_box(x,y,thrdv<dv&thru<u)
y=1
m =1.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrdv<dv,col="skyblue")
plot_range_box(x,y,thru<u,col="pink")
plot_range_box(x,y,thrdv<dv&thru<u)
y=0
m =10.0
rv=2
dv = (ampl(m,rv)*cos((x-15-delay(m,rv))/24*2*acos(-1)))
plot_range_box(x,y,thrdv<dv,col="skyblue")
plot_range_box(x,y,thru<u,col="pink")
plot_range_box(x,y,thrdv<dv&thru<u)

dev.off()
