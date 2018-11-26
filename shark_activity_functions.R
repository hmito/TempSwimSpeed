#load library
library("Rcpp")
library("BH")
sourceCpp("TempSwimSpeedOptim.cpp")

#find peaks from histgram data 
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

#calculate bodytemp based on the body mass (effective body size on the watertemp dynamics)
#param
#	watertemp: sequence of water temperature
#	mass: effective body size in the context of heat balance
#		e.g., basically, it is similar with the body size; i.e., small/large body fish is quickly/slowley heated
#		e.g., thin body shape reduce effective body size because their surface area is relatively larger
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

#calculate mass from sharkradius & skinthickness
#	mathematically same with calc.shark.temp function if we use
#  	calc.bodytemp(watertemp, calc.mass_for_bodytemp(sharkradius,skinthickness))
calc.mass_for_bodytemp = function(sharkradius,skinthickness){
	1/(2*0.6*0.031593/60/(sharkradius*log(sharkradius/(sharkradius-skinthickness))))
}

#calculate light effect with old definition
#	the predation rate is the exponential of cosin curve
calc.light_effect_old = function(t,dark,light){
	lwave=exp(-1.0*cos(2*pi*(length(t)/2-t)/length(t)))/exp(1.0)
	return(dark+(light-dark)*lwave)
}

#calculate light effect with new definition
#	the light level is the cosin curve + constant value for twilight
#	the light level less than zero is ignored (just considred as zero)
#	the predation rate is exponential of light level
#param
#	t: time sequence
#	l_min: minimum predation rate
#	l_max: maximum predation rate
#	light_influence: strength of light influence on predation (how strong light is required for the reduction of predation rate)
#	twilight_coef: influence of twilight. 0.3 seems good (see comments inside of the function)
#return
#	sequence of predatation rate at each time
calc.light_effect = function(t, l_min, l_max, light_influence, twilight_coef=0.3){
	#Around twilight_coef = 0.3 seems to be reasonable because
	#	- the definition of astronominal twilight is that the sun is less than -18 degree below the horizon. 
	#	- When the Culmination altitute = 90 degree, 18 degree passes for 1.2 hours, so twilight start from t = 4.8 and finish at t=19.2.
	#	- This time becomes longer when the Culmination altitude is less than 90 degree.
	#	- At twilight_coef = 0.3,
	#		- the sea is perfectly dark when t is from 20.5 to 3.5.
	#		- the sea is slightly bright (3% of noon) at t=4.5, 19.5
	#		- the light level rach to 20% and 40% of noon at t=5.5, 6.5.
	
	lwave=(1-twilight_coef)*cos(2*pi*(t-12)/length(t)) + twilight_coef
	lwave[lwave<0] = 0
	alpha = log(l_max/l_min)
	return(l_max*exp(-alpha*lwave^(1/light_influence)))
}

#plot figure of the simulation result
#	x:time
#	y:foraging frequency of prey(blue) and predator(red))
#parameters
#	Ans: return value of tss_probforage_energygain_optimize function
#	title: title of the figure
#return
#	none
plot.sim_result = function(Ans,title){
	Prey= Ans$Prey
	Predator = Ans$Predator
	plot(0,0,type="n",
		  #		  xlab="time (t)",ylab="foraging predator (red), prey (blue)",
		  xlab="",ylab="",
		  xlim=c(0,tnum),ylim=c(-0.02,1.02),       
		  main=title
	)
	lines(t,Prey,col="blue",lwd=3)
	lines(t,Predator*0.99,col="red",lwd=3,lty="dashed")
}

#plot figure summrizing assumption
#	water temp and shark temp in left panel
#	V, U and L in right panel
#parameters
#	t: sequence of time
#	watertemp: water temperature
#	sharktemp: shark body temperature
#	V: predator swim speed
#	U: prey swim speed
#	L: light effect
#return
#	none
plot.assumption=function(t,watertemp,sharktemp,V,U,L){
	# assumption figure
	par(mfrow=c(1,2))
	plot(t,watertemp,col="blue",pch=15, xlab="time (t)",ylab="temperature")
	lines(t,watertemp,col="blue",lty=1)
	#	text(7,25.0,bquote('w'['t']))
	points(t,sharktemp,col="red",pch=16)
	lines(t,sharktemp,col="red",lty=1)
	#points(t,calc.bodytemp(watertemp,4),col="purple",pch=16)
	#lines(t,calc.bodytemp(watertemp,4),col="purple",lty=1)
	#points(t,calc.bodytemp(watertemp,25),col="orange",pch=16)
	#lines(t,calc.bodytemp(watertemp,25),col="orange",lty=1)
	#	text(4,28,bquote('s'['t']))
	par(new =T)
	plot(t,L,col="orange",pch=17, type = "p", axes = FALSE, ylab = "",ylim=c(0,max(L)),main="temperature")
	lines(t,L,col="orange",lty=1)
	#	text(5,0.25,bquote(lambda['t']))
	axis(4)
	
	plot(t,V,type="n",xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),
		  xlab="time (t)",ylab="burst speed",main="swim speed")
	segments(-100,0,100,0)
	points(t,V,col="red",pch=16)
	lines(t,V,col="red",lty=1)
	#	text(15,2.3,bquote('v'['t']))
	points(t,U,col="blue",pch=15)
	lines(t,U,col="blue",lty=1)
	#	text(18,1.9,bquote('u'['t']))
	points(t,V-U,col="purple",pch=2)
	lines(t,V-U,col="purple",lty=1)
	#	text(17,0.9,bquote(list('v'['t'],'- u'['t'])))
}

#major-active-time based categorization with 4 time zone
#parameters
#	Ans: return value of tss_probforage_energygain_optimize function
#	MajorProb: threshold for the definition of major-acitve time
#		major-active time := [foraging frequency during the focal duration]/[total foraging frequency] >= MajorProb
majortime.get_category = function(Ans, MajorProb = 0.20){
	ans = 0
	
	if(sum(Ans$Predator)!=0){
		#active time
		#	t=0-6: group 1
		#	t=6-12: group 2
		#	t=12-18: group 3
		#	t=18-24: group 4
		# 1, 2, 4, 8 means the active time is group 1,2,3,4
		# multiple group use is just shown by the summation
		#	gruop 1+3: 1+4
		#	group 2+3+4: 2+4+8
		#	no use: 0
		#	full use: 15
		{
			p = Ans$Predator
			group = rep(1:4,times=c(6,6,6,6))
			
			for(igroup in unique(group)){
				if(sum(p[group==igroup])/sum(p) > MajorProb){
					ans = ans + 2^(igroup-1)
				}
			}
		}
		
		#peak number
		#	16 * (peaknum - 1)
		#	0-15 : one peak
		#	16-31: two peaks
		#	32-47: three peaks
		#	...
		{
			thr = 0.5			
			peaks = hist.find_peaks(p,0,thr)
			
			#connect 24-1
			if(nrow(peaks)>1 && peaks$lower[1]==1 && peaks$upper[nrow(peaks)]==24){
				peaks$lower[1] = peaks$lower[nrow(peaks)]
				peaks$freq[1] = peaks$freq[1] + peaks$freq[nrow(peaks)]
				peaks = peaks[-nrow(peaks),]
			}
			ans = ans + 16*(nrow(peaks)-1)
		}
	}
	
	return(ans)
}

#transform majortime category to plotmode, which define the colour of plot
#parameters
#	catergory: return value of majortime.get_category
#return
#	plotmode (used in image.plotmode)
majortime.get_plotmode = function(category){
	#colour for categories which is not allocated by any colour 
	errclr = "grey"
	
	clr = category
	clr[category==category] = errclr #error
	clr[category== 0]  = "black"	#nothing
	clr[category== 15] = "white"	#all
	clr[category== 6]  = "red"	#day
	clr[category== 8]  = "purple"	#early-night
	clr[category== 12] = "yellow"	#pm
	clr[category== 11] = "skyblue"	#except pm-day
	clr[category== 14] = "orange"	#except am-night
	clr[category== 9]  = "blue"	#night
	clr[category== 13] = "mediumpurple1"	#except am-day
	clr[category== 7] = "saddlebrown"	#except pm-night
	clr[category== 1] = "blue4"	#am-night
	clr[category== 4] = "red3"	#pm-day
	clr[category>=16 & category<32]  = "forestgreen"	#two peaks
	clr[category>=32]  = "green"	#three peaks
	
	dotclr = category
	dotclr[category==category] = "none"
	
	#list of error category (not allocated by any colour )
	err = sort(unique(as.vector(category[clr == errclr])))
	
	return(list(clr=clr,dotclr=dotclr,err_category=err))
}

#major-active-time based categorization with 5 time zone
#parameters
#	Ans: return value of tss_probforage_energygain_optimize function
#	MajorProb: threshold for the definition of major-acitve time
#		major-active time := [foraging frequency during the focal duration]/[total foraging frequency] >= MajorProb
majortime5.get_category = function(Ans, MajorProb = 0.20){
	ans = 0
	
	if(sum(Ans$Predator)!=0){
		#active time
		#	t=0-4: group 1
		#	t=5-8: group 2
		#	t=9-15: group 3
		#	t=16-19: group 4
		#	t=20-23: group 5
		# 1, 2, 4, 8, 16 means the active time is group 1,2,3,4, 5
		# multiple group use is just shown by the summation
		#	gruop 1+3: 1+4
		#	group 2+3+4: 2+4+8
		#	no use: 0
		#	full use: 31
		{
			p = Ans$Predator
			group = c(rep(1,5),rep(2,4),rep(3,7),rep(4,4),rep(5,4))
				
			for(igroup in unique(group)){
				if(sum(p[group==igroup])/sum(p) > MajorProb){
					ans = ans + 2^(igroup-1)
				}
			}
		}
		
		#peak number
		#	100 * (peaknum - 1)
		#	0-31 : one peak
		#	100-131: two peaks
		#	200-231: three peaks
		#	...
		{
			thr = 0.5			
			peaks = hist.find_peaks(p,0,thr)
			
			#connect 24-1
			if(nrow(peaks)>1 && peaks$lower[1]==1 && peaks$upper[nrow(peaks)]==24){
				peaks$lower[1] = peaks$lower[nrow(peaks)]
				peaks$freq[1] = peaks$freq[1] + peaks$freq[nrow(peaks)]
				peaks = peaks[-nrow(peaks),]
			}
			ans = ans + 100*(nrow(peaks)-1)
		}
	}
	
	return(ans)
}

#transform majortime5 category to plotmode, which define the colour of plot
#parameters
#	catergory: return value of majortime.get_category
#return
#	plotmode (used in image.plotmode)
majortime5.get_plotmode = function(category){
	#colour for categories which is not allocated by any colour 
	errclr = "grey"
	
	#remove information of peak numbers
	majortime=(category%%100)

	#definition of image colours
	clr = category
	clr[category==category] = errclr #error
	#add allocation of colors in the following lines
	clr[majortime== 0]  = "black"	#nothing
	clr[majortime== 31] = "white"	#all
	clr[majortime== 1]  = "blue"	#nocturnal  (0-4)
	clr[majortime== 17]  = "blue"	#nocturnal (20-4)
	clr[majortime== 16]  = "purple"	#early night (20-23)
	#etc...

	#definition of dot colours
	dotclr = category
	dotclr[category==category] = "none"
	#add allocation of dot colours in the following lines
	dotclr[100<=category & category<200] = "black"
	dotclr[200<=category & category<300] = "white"
	
	#list of error category (not allocated by any colour )
	err = sort(unique(as.vector(category[clr == errclr])))
	
	return(list(clr=clr,dotclr=dotclr, err_category=err))
}

#plot category with different colours following the definition in plotmode
#	x.seq: sequence of xaxis
#	y.seq: sequence of yaxis
#	plotmode: matrix of plotmode (return value of majortime.get_plotmode or majortime5.get_plotmode)
#	dotcex: change the scale of dots 
image.plotmode=function(x.seq,y.seq,plotmode,dotcex = 0.33,...){
	fz = factor(plotmode$clr)
	z = matrix(as.integer(fz),length(x.seq),length(y.seq))
	clr = levels(fz)
	clr[clr=="none"] = rgb(0,0,0,0)
	image(x.seq,y.seq,z,col=levels(fz),...)
	
	fd = factor(plotmode$dotclr)
	for(clr in levels(fd)){
		if(clr == "none")next
		
		x.mx = matrix(rep(x.seq,times = length(y.seq)), length(x.seq), length(y.seq))
		y.mx = matrix(rep(y.seq,each  = length(x.seq)), length(x.seq), length(y.seq))
		
		z = (fd==clr)
		xpos = x.mx[z]
		ypos = y.mx[z]
		
		points(xpos,ypos,col = clr, pch=16, cex=dotcex)
	}
}
