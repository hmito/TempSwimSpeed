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

#TempSwimSpeed_v1_02
#   Predator efficiency 
library("Rcpp")
library("BH")
sourceCpp("TempSwimSpeedOptim.cpp")

pi = acos(-1)

tnum = 24
t = 0:23+0.5

vmin = 1.5
vmax = 2.5
vmaxt = 18

umin = 0.5
umax = 2.0
umaxt = 15

#Define the speed of predator (V) and prey (U) at each time step
V = vmin + (vmax - vmin)*(cos(2*pi*(t-vmaxt)/tnum)+1)/2
U = umin + (umax - umin)*(cos(2*pi*(t-umaxt)/tnum)+1)/2


K = rep(0.02, length=tnum)		#amount of food availability for prey
C = rep(0.015, length=tnum)		#metabolic cost for predators when they go out for predation

#following three parameters determine the prey traits
e = 0.3		#relative risk of predation for resting prey
d = 0.001	#relative density of predator/prey
#following three parameters determine the predation rate: a*(v-u)^b / {1 + h*a*(v-u)^b} 
b = 3.0		#non-linear influence of speed difference
h = 2.0	  #handling time for predation a prey
lmean=0.5 # minimum risk (in the dark)
lamp=1.0 # amplitude of light

#influence of light on the predation rate
L = 0.2+0.8*exp(-2.0*cos(2*pi*t/tnum))/exp(2.0)	
L = lmean+lamp*(exp(-2.0*cos(2*pi*t/tnum))/exp(2.0)-0.5)

#following twp parameters determine the cost of prey
cb = 0.0001	#metabolic cost for prey (should pay both for resting and foraging)
cf	= 0.0001	#foraging cost for prey (should pay only for foraging)

# assumption figure
plot(t,V,col="red",xlim=c(0,tnum),ylim=c(0,max(c(vmax,umax))))
lines(t,V,col="red")
points(t,U,col="blue")
lines(t,U,col="blue")
points(t,V-U,col="purple")
lines(t,V-U,col="purple")
points(t,(V-U)^3,col="green")
lines(t,(V-U)^3,col="green")
points(t,L,col="black")
lines(t,L,col="black")
x11()

par(mfrow=c(4,2))
figlabs=c('safe rest, no effect of temperature','dangerous rest, no effect of temperature','safe rest, effect of temperature','dangerous rest & effect of temperature')

# loop over the options 
for (i in 1:8) {
	
	#Variation of model
	#if we set d=0, the prey decide their behaviour only based on the food availability for them regardless of the predation.
	# d = 0.0
	#
	#if we set b=0, the swim speed (i.e., temperature) has no influence on the predation rate.
	# b = 0.0
	#
	#if fix the V and U as constant values, the influence of temperature is ignored in this model.
	# V=seq(1.5, length=tnum)
	# U=seq(1.2, length=tnum)
	
	lamp=0.6*((i+1) %% 2)
	b=3.0*(i>4)
	
	e=0.125+0.125*(((i+1) %% 4)<2)
	#d = 0.001+0.01*i
	
	
	#influence of light on the predation rate
	L = lmean+lamp*(exp(-2.0*cos(2*pi*t/tnum))/exp(2.0)-0.5)
	
	# find solution
	Ans = tss_probforage_energygain_optimize(V, U, K, C, L, d, e, b, h,cb,cf)
	
	#Summarized figure
	Prey= Ans$Prey
	Predator = Ans$Predator

	plot(0,0,type="n",
		  xlab="time (t)",ylab="foraging predator (red), prey (blue)",
		  xlim=c(0,tnum),ylim=c(-0.02,1.02),       
		  main=bquote(paste('l'['1']*'=',.(lamp),', b=', .(b),', e=', .(e))))  #  main=figlabs[i]
	lines(t,Predator*0.98,col="red",lwd=3)
	lines(t,Prey,col="blue",lwd=3)
} # end of loop 


dev.off()

#=== analysis of peaks1 ===
#	example of the analysis counting the number of peaks of predator.
#	Two parameter (lamp and V-U) is changed in x-axis and y-axis.
DotNum=41

x.seq = seq(0,0.5,length=DotNum)
y.seq = seq(0,1,length=DotNum)

PreyPeaks.mx = matrix(0,DotNum,DotNum)
PredatorPeaks.mx = matrix(0,DotNum,DotNum)

for (y.i in 1:length(y.seq)) {
	for (x.i in 1:length(x.seq)) {
		#set the parameter valus
		# x: effect of lamp
		# y: difference of vmax - umax, vmin-umin
		lamp = x.seq[x.i]
		vudif = y.seq[y.i]
		
		#influence of light on the predation rate
		L = lmean+lamp*(exp(-2.0*cos(2*pi*t/tnum))/exp(2.0)-0.5)
		
		#influence of speed
		umin = 0.5
		umax = 2.0
		umaxt = 15
		
		vmin = umin + vudif
		vmax = umax + vudif
		vmaxt = 18
		
		V = vmin + (vmax - vmin)*(cos(2*pi*(t-vmaxt)/tnum)+1)/2
		U = umin + (umax - umin)*(cos(2*pi*(t-umaxt)/tnum)+1)/2
		
		# find solution
		Ans = tss_probforage_energygain_optimize(V, U, K, C, L, d, e, b, h,cb,cf)

		#find_peaks function return the table of the peaks
		#	Each row of the table shows infomation of each peak
		#		Therefore, nrow (number of row) is the number of peaks 
		#	"lower" col and "upper" cik show the lower and upper value (in this case, time) of the focal peak
		#		Therefore, upper-lower is the width of each peak
		#	"top" and "freq" cols will be useless in this model
		#		top is the maximum point (time) of the peak
		#		freq is the total value inside of the peak
		PreyPeaks = hist.find_peaks(Ans$Prey,0.00,0.05)
		PredatorPeaks = hist.find_peaks(Ans$Predator,0.00,0.05)
		PreyPeaks.mx[x.i,y.i] = nrow(PreyPeaks) + ifelse(PreyPeaks$lower[1]==1 & PreyPeaks$upper[nrow(PreyPeaks)]==tnum, -1, 0)
		PredatorPeaks.mx[x.i,y.i] = nrow(PredatorPeaks) + ifelse(PredatorPeaks$lower[1]==1 & PredatorPeaks$upper[nrow(PredatorPeaks)]==tnum, -1, 0)
		#ifelse is used for the reduction of the number of peaks because time:0 == time:24, i.e, the two peaks are actually only one peak.
	}
}

#colour shows the number of peaks
filled.contour(x.seq,y.seq,PredatorPeaks.mx,col=grey.colors(20))


#=== analysis of peaks2 ===
#	Map of the predator-prey strategy plot.
#	output is png file with defined name
DotNum=10

x.seq = seq(0,0.5,length=DotNum)
y.seq = seq(0,1,length=DotNum)

PredatorMode.mx = matrix(0,DotNum,DotNum)
png("test.png",height=2500,width=2500)
par(oma=c(0.2,0.2,0.2,0.2),mfrow=c(DotNum,DotNum))
for (y.i in length(y.seq):1) {
	for (x.i in 1:length(x.seq)) {
		#set the parameter valus
		# x: effect of lamp
		# y: difference of vmax - umax, vmin-umin
		lamp = x.seq[x.i]
		vudif = y.seq[y.i]
		
		#influence of light on the predation rate
		L = lmean+lamp*(exp(-2.0*cos(2*pi*t/tnum))/exp(2.0)-0.5)
		
		#influence of speed
		umin = 0.5
		umax = 2.0
		umaxt = 15
		
		vmin = umin + vudif
		vmax = umax + vudif
		vmaxt = 18
		
		V = vmin + (vmax - vmin)*(cos(2*pi*(t-vmaxt)/tnum)+1)/2
		U = umin + (umax - umin)*(cos(2*pi*(t-umaxt)/tnum)+1)/2
		
		# find solution
		Ans = tss_probforage_energygain_optimize(V, U, K, C, L, d, e, b, h,cb,cf)
		
		#Summarized figure
		Prey= Ans$Prey
		Predator = Ans$Predator
		
		plot(0,0,type="n",
			  xlab="time (t)",ylab="foraging predator (red), prey (blue)",
			  xlim=c(0,tnum),ylim=c(-0.02,1.02),       
			  main=bquote(paste('l'['1']*'=',.(lamp),', b=', .(b),', e=', .(e))))  #  main=figlabs[i]
		lines(t,Predator*0.98,col="red",lwd=3)
		lines(t,Prey,col="blue",lwd=3)
	}
}
dev.off()



#=== analysis of peaks3 ===
#	example of the analysis based on the pattern of the predator active time.
#	Two parameter (lamp and V-U) is changed in x-axis and y-axis.
DotNum=41

x.seq = seq(0,0.5,length=DotNum)
y.seq = seq(0,1,length=DotNum)

PredatorMode.mx = matrix(0,DotNum,DotNum)

for (y.i in 1:length(y.seq)) {
	for (x.i in 1:length(x.seq)) {
		#set the parameter valus
		# x: effect of lamp
		# y: difference of vmax - umax, vmin-umin
		lamp = x.seq[x.i]
		vudif = y.seq[y.i]
		
		#influence of light on the predation rate
		L = lmean+lamp*(exp(-2.0*cos(2*pi*t/tnum))/exp(2.0)-0.5)
		
		#influence of speed
		umin = 0.5
		umax = 2.0
		umaxt = 15
		
		vmin = umin + vudif
		vmax = umax + vudif
		vmaxt = 18
		
		V = vmin + (vmax - vmin)*(cos(2*pi*(t-vmaxt)/tnum)+1)/2
		U = umin + (umax - umin)*(cos(2*pi*(t-umaxt)/tnum)+1)/2
		
		# find solution
		Ans = tss_probforage_energygain_optimize(V, U, K, C, L, d, e, b, h,cb,cf)
		
		#find_peaks function return the table of the peaks
		#	Each row of the table shows infomation of each peak
		#		Therefore, nrow (number of row) is the number of peaks 
		#	"lower" col and "upper" cik show the lower and upper value (in this case, time) of the focal peak
		#		Therefore, upper-lower is the width of each peak
		#	"top" and "freq" cols will be useless in this model
		#		top is the maximum point (time) of the peak
		#		freq is the total value inside of the peak
		PredatorPeaks = hist.find_peaks(Ans$Predator,0.00,0.05)

		#Mode
		#	-1 "black": unknown 
		#	 0 "grey": No predation
		#	 1 "blue":	afternoon 
		#	 2 "cyan": afternoon-midnight(over 12)
		#	 3 "green":	beforenoon-midnight (rest only morning)
		#	 4 "orange": afternoon with short rest on everning
		#	 5 "yellow": afternoon-midnight(over 12) with short rest on everning 
		#	 6 "red": afternoon-midnight(over 12) with short rest on morning 
		Mode=-1
		if(nrow(PredatorPeaks)==0){
			Mode = 0
		}else if(nrow(PredatorPeaks)==1){
			if(12 < PredatorPeaks$lower[1] && PredatorPeaks$upper[1] <= 24){
				#predation is start after 0 pm and finish before 0am
				Mode = 1
			}
		}else if(nrow(PredatorPeaks)==2){
			if(PredatorPeaks$lower[1]==1 && PredatorPeaks$upper[1] < 9 &&
				12 < PredatorPeaks$lower[2] && PredatorPeaks$upper[2] == 24){
				#predation is start after 0pm and finish between [0pm : 9am].
				Mode = 2
			}else if(PredatorPeaks$lower[1]==1 && PredatorPeaks$upper[1] < 9 &&
				PredatorPeaks$upper[1] < PredatorPeaks$lower[2] && PredatorPeaks$upper[2] == 24){
				#predation is start beforenoon[9am:0pm], and finish before 9 am
				Mode = 3
			}else if(12 < PredatorPeaks$lower[1] && PredatorPeaks$upper[1] < 24 &&
				PredatorPeaks$upper[1] < PredatorPeaks$lower[2] && PredatorPeaks$upper[2] <= 24){
				#predation is start after 0pm, and finish before 0am with short rest.
				Mode = 4
			}
		}else if(nrow(PredatorPeaks)==3){
			if(PredatorPeaks$lower[1]==1 && PredatorPeaks$upper[1] < 9 &&
				12 < PredatorPeaks$lower[2] && PredatorPeaks$upper[2] < 24&&
				PredatorPeaks$upper[2] < PredatorPeaks$lower[3] && PredatorPeaks$upper[3] == 24){
				#predation is start after 0pm, and finish between [0pm : 9am]. short rest at between [0pm : 0am].
				Mode = 5
			}else if(PredatorPeaks$lower[1]==1 && PredatorPeaks$upper[1] < 9 &&
				PredatorPeaks$upper[1] < PredatorPeaks$lower[2] && PredatorPeaks$upper[2] < 9&&
				12 < PredatorPeaks$lower[3] && PredatorPeaks$upper[3] == 24){
				#predation is start after 0pm, and finish between [0pm : 9am]. short rest at between [0pm : 9am].
				Mode = 6
			}
		}
		PredatorMode.mx[x.i,y.i] = Mode
	}
}

image(x.seq,y.seq,PredatorMode.mx,zlim=c(-1,6),col=c("black","grey","blue","cyan","green","orange","yellow","red"))
