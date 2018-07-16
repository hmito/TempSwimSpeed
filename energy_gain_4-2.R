hist.find_peaks = function(hist,min,max,n=101){
	threshold.mx = matrix(rep(seq(min,max,length=n),each=length(hist)+2),length(hist)+2,n)
	hist.mx = matrix(rep(c(0,hist,0),times=n),length(hist)+2,n)
	exist.mx = hist.mx>threshold.mx
	
	result = apply(exist.mx[-1,]!=exist.mx[-nrow(exist.mx),],2,sum)/2
	result.table = table(result)
	peak.num = as.integer(names(result.table)[order(result.table,decreasing=TRUE)[1]])
	
	threshold.no = min((1:n)[result==peak.num])
	exist.seq = exist.mx[,threshold.no]
	
	peaks = data.frame("lower"=NA,"upper"=NA,"top"=NA,"freq"=NA)

	if(peak.num!=0){
		boundary = 1
		
		lower = (1:(length(exist.seq)-1))[exist.seq[-1]&(!exist.seq[-length(exist.seq)])]
		upper = (0:(length(exist.seq)))[(!exist.seq[-1])&exist.seq[-length(exist.seq)]]
		
		for(peak.pos in 1:peak.num){
			this.lower = min((boundary:lower[peak.pos])[hist[boundary:lower[peak.pos]]>0])
			if(peak.pos<peak.num){
				boundary = (upper[peak.pos]:lower[peak.pos+1])[order(hist[upper[peak.pos]:lower[peak.pos+1]])[1]]
			}else{
				boundary = length(hist)
			}
			this.upper = max((upper[peak.pos]:(boundary-1))[hist[upper[peak.pos]:(boundary-1)]>0])
			this.top = (this.lower:this.upper)[order(hist[this.lower:this.upper],decreasing = TRUE)[1]]
			this.freq = sum(hist[this.lower:this.upper])
			
			peaks = rbind(peaks,data.frame("lower"=this.lower,"upper"=this.upper,"top"=this.top,"freq"=this.freq))
		}
	}

	return(peaks[-1,])
}

calc.shark.temp = function(tnum,watertemp,sharkradius,skinthickness){
  # calculate shark temp by iteration
  sharktemp=seq(25,25,length.out=tnum)
  for (ii in c(1:20)) {
    for (jj in c(1:24)) {
      ii
      jj
      # from conduction physics
      Q=-2*0.6*(sharktemp[jj]-watertemp[jj])/(sharkradius*log(sharkradius/(sharkradius-skinthickness)))
      # result is in Watts per minute
      # change to celsius per hour
      Q=Q*0.031593/60 
      if (jj==24){sharktemp[1]=sharktemp[jj]+Q}
      if (jj<24){sharktemp[jj+1]=sharktemp[jj]+Q}
      #    print(jj)
      #    print(Q)
      #    print(sharktemp[jj]+Q)
      #    Sys.sleep(0.01)
    }
  }
  return(sharktemp)
}

# 2x2x2 comparison
library("Rcpp")
library("BH")
sourceCpp("TempSwimSpeedOptim.cpp")

pi = acos(-1)
tnum = 24 # time of day
t = 0:23+0.5
# influence of light on the predation rate
lmin=0.0 # lowest
lmax=1.0 # highest
lwave=exp(-1.0*cos(2*pi*(tnum/2-t)/tnum))/exp(1.0)
# temperature of the water
tempmaxt=15
tmin = 25
tmax = 30
watertemp=tmin+(tmax-tmin)*(cos(2*pi*(t-tempmaxt)/tnum)+1)/2
#watertemp=seq(0,0,length.out=length(watertemp))+tmin+(tmax-tmin)/2
# prey immediately track temperature
umin = 1.0
utemp = 0.25
U = umin + utemp*(watertemp-min(watertemp))
# work out the shark's body temperature
sharkradius=0.2
skinthickness=sharkradius/20
sharktemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)
# sharks track their own temperature
vmin = 2.0
vtemp = utemp
V = vmin + vtemp*(sharktemp-min(sharktemp))
#amount of food availability for prey
alpha = rep(0.02, length=tnum)
#[NEW 18/07/16] influence of swim speed on foraging efficiency
omega = 0	#obtained reward is alpha*((1-omega)+omega*U)
#metabolic cost for predators when they go out for predation
C = rep(0.02, length=tnum)		
#following three parameters determine the prey traits
phi = 0.15		#relative risk of predation for resting prey
#following three parameters determine the predation rate: a*(v-u)^beta / {1 + h*a*(v-u)^beta} 
beta = 1.0		#non-linear influence of speed difference
h = 2.0	  #handling time for predation a prey
#following twp parameters determine the cost of prey
mb = 0.01	# basleline cost for prey (should pay both for resting and foraging)
mx	= 0.01	# foraging cost for prey (should pay only for foraging)
my = 0.1	#relative density of predator/prey

# assumption figure
par(mfrow=c(1,2))
plot(t,watertemp,col="blue",pch=15,,xlim=c(0,tnum),ylim=c(20,max(c(watertemp,sharktemp))),
     xlab="time (t)",ylab="temperature")
lines(t,watertemp,col="blue",lty=1)
text(5,24.5,bquote('w'['t']))
points(t,sharktemp,col="red",pch=16)
lines(t,sharktemp,col="red",lty=1)
text(4,28,bquote('s'['t']))
points(t,lmin+(lmax-lmin)*lwave*10+20,col="orange",pch=17)
lines(t,lmin+(lmax-lmin)*lwave*10+20,col="orange",lty=1)
text(4.5,21.5,bquote(lambda['t']))

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
#points(t,(V-U)^3,col="green",pch=17)
#lines(t,(V-U)^3,col="green",lty=1)
#text(20,0.4,bquote(list('(v'['t'],'- u'['t'],')'^3)))

x11()
par(mfrow=c(4,2))
# loop over the options 
for (i in 1:8) {
  #if we set beta=0, the swim speed (i.e., temperature) has no influence on the predation rate.
	beta = 1.0*(i>4)
#	lmin = 0.5-0.5*((i+1) %% 2)
#	lmax = 0.5+0.5*((i+1) %% 2)
  # vary the magnitude of other predation
#	mx = 0.01*((i+1) %% 2)  
	alpha = rep(0.01+0.09*(((i+1) %% 2)), length=tnum)
	# vary the magnitude of shark predation
	my = 0.1*(((i+1) %% 4)<2)
	# find solution
	Ans = tss_probforage_energygain_optimize(V, U, alpha, C, lmin+(lmax-lmin)*lwave, my, phi, omega, beta, h,mb,mx)
	#Summarized figure
	Prey= Ans$Prey
	Predator = Ans$Predator
	plot(0,0,type="n",
		  xlab="time (t)",ylab="foraging predator (red), prey (blue)",
		  xlim=c(0,tnum),ylim=c(-0.02,1.02),       
	    main=bquote(list('l'[1] == .(lmax-lmin), beta==.(beta),phi==.(phi),'m'['y']==.(my),'m'['b']==.(mb),'m'['x']==.(mx))))
	lines(t,Predator*0.99,col="red",lwd=3,lty="dashed")
	lines(t,Prey,col="blue",lwd=3)
} # end of loop 

#=== analysis of peaks3 ===
#	example of the analysis based on the pattern of the predator active time.
#	Two parameter (lmax and V-U) is changed in x-axis and y-axis.
DotNum=31

x.seq = seq(0,0.1,length=DotNum)
y.seq = seq(0,0.1,length=DotNum)

PredatorMode.mx = matrix(0,DotNum,DotNum)
PredatorFitness.mx = matrix(0,DotNum,DotNum)
# plot some of the results as whole strategy
png("test.png",height=2500,width=2500)
par(oma=c(0.2,0.2,0.2,0.2),mfrow=c(4,4))

for (y.i in 1:length(y.seq)) {
	for (x.i in 1:length(x.seq)) {
	  # reset baseline
	  lmin = 0.0
	  lmax = 1.0
	  umin = 1.0
	  utemp = 0.25
	  vmin = 2.0
	  vtemp = utemp
	  phi = 0.15
	  my = 0.1
	  mb = 0.01	
	  mx	= 0.0	
	  sharkradius=0.2
	  
	  # FIGURE 4 range of light and difference in speed (x=0:1, y=0:1)
#	  sharkradius=x.seq[x.i]##
#	  vmin = y.seq[y.i]
#	  thisxlab="shark radius"
#	  thisylab="shark speed minimum"
	  
	  # FIGURE 5 influence of speed sensitivity to temperature
#	  	  lmin = 0.5 - x.seq[x.i]/2
#	     lmax = 0.5 + x.seq[x.i]/2
#	  C = rep(x.seq[x.i], length=tnum)	#
#      utemp=y.seq[y.i]
#	    vtemp = utemp
#	    thisxlab="cost of hunting"
#	      thisylab="sensitivity of speed to temperature"
	  
	  # FIGURE 6 influence of mortality rates (x=0:0.09, y=0:0.09)
   # phi = x.seq[x.i]##
	  beta = y.seq[y.i]##
	  my=x.seq[x.i]##
	  thisxlab="predation risk from sharks"
    thisylab="predation depends on speed difference "
	  
	  # correct the things
    skinthickness=sharkradius/20
    sharktemp=calc.shark.temp(tnum,watertemp,sharkradius,skinthickness)
    U = umin + utemp*(watertemp-min(watertemp))
    V = vmin + vtemp*(sharktemp-min(sharktemp))
	  
	  # find solution
    Ans = tss_probforage_energygain_optimize(V, U, alpha, C, lmin+(lmax-lmin)*lwave, my, phi, omega, beta, h,mb,mx)
	  
	  if (round((y.i-1)/((length(y.seq)-1)/3))==((y.i-1)/((length(y.seq)-1)/3)))
	  { if (round((x.i-1)/((length(x.seq)-1)/3))==((x.i-1)/((length(x.seq)-1)/3)))
	  {
	    x.i
	    y.i
	    plot(0,0,type="n",
	       xlab="time (t)",ylab="foraging predator (red), prey (blue)",
	       xlim=c(0,tnum),ylim=c(-0.02,1.02),       
	       main=bquote(paste(x,'=',.(x.seq[x.i]),', y=', .(y.seq[y.i]))))  #  main=figlabs[i]
	    lines(t,Ans$Predator*0.98,col="red",lwd=3)
	    lines(t,Ans$Prey,col="blue",lwd=3)
	  }
	  }
	  
	  #find_peaks function return the table of the peaks
		#	Each row of the table shows infomation of each peak
		#		Therefore, nrow (number of row) is the number of peaks 
		#	"lower" col and "upper" cik show the lower and upper value (in this case, time) of the focal peak
		#		Therefore, upper-lower is the width of each peak
		#	"top" and "freq" cols will be useless in this model
		#		top is the maximum point (time) of the peak
		#		freq is the total value inside of the peak
		PredatorPeaks = hist.find_peaks(Ans$Predator,0.00,0.05)
		Mode=0
		if(nrow(PredatorPeaks)>0) {
	      Mode=nrow(PredatorPeaks)
	      if(Mode>3){Mode=3}
		    if(Ans$Predator[1]>0.5){
		      if(Ans$Predator[tnum]>0.5) {
            Mode=Mode+2
		      } else { 
		        Mode=Mode+3
		      }
		    }
		}
		if (sum(Ans$Predator)>tnum-1) {Mode=7}
		PredatorMode.mx[x.i,y.i] = Mode
		PredatorFitness.mx[x.i,y.i] = sum(Ans$Predator*(Ans$PredatorReward1-Ans$PredatorCost))
	}
}
dev.off()
#Mode
#	 0 "black": no foraging
#	 1 "blue":	one peak during day 
#	 2 "cyan": two peaks during day
#	 3 "green":	three peaks during day
#	 4 "red": one peak over night
#	 5 "orange": two peaks, one overnight
#	 6 "yellow": three peaks, one overnight 
#  7 "white": constant foraging
x11()
par(mfrow=c(1,2))
image(x.seq,y.seq,PredatorMode.mx,zlim=c(0,7),xlab=thisxlab,ylab=thisylab,col=c("black","blue","cyan","green","red","orange","yellow","white"))
image(x.seq,y.seq,PredatorFitness.mx,xlab=thisxlab,ylab=thisylab, axes = FALSE)
contour(x.seq,y.seq,PredatorFitness.mx, add = TRUE, drawlabels = TRUE,levels = seq(min(PredatorFitness.mx), max(PredatorFitness.mx), by = 0.1))
