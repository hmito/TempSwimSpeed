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

#calculate bodytemp based on the body mass by iterated numerical calculation (effective body size on the watertemp dynamics)
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

#analytically calc bodytemp based on the body mass (effective body size on the watertemp dynamics)
#param
#	t: sequence of time
#	tw: time of peak (0-24)
#	wmin: minimum water temperature
#	wmax: maximum water temperature
#	mass: effective body size in the context of heat balance
#		e.g., basically, it is similar with the body size; i.e., small/large body fish is quickly/slowley heated
#		e.g., thin body shape reduce effective body size because their surface area is relatively larger
get.bodytemp = function(t, tw, wmin, wmax, mass){
	V = 1 / sqrt(1 + (2*acos(-1)*mass/24)^2)
	return((wmax+wmin)/2 + (wmax-wmin)/2*V*cos(2*acos(-1)*(t-tw)/24 - acos(V)))
}

#calculate mass from sharkradius & skinthickness
#	mathematically same with calc.shark.temp function if we use
#  	calc.bodytemp(watertemp, calc.mass_for_bodytemp(sharkradius,skinthickness))
calc.mass_for_bodytemp = function(sharkradius,skinthickness){
	1/(2*0.6*0.031593/60/(sharkradius*log(sharkradius/(sharkradius-skinthickness))))
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
calc.light_effect = function(t, l_min, light_influence, twilight_coef){
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
	alpha = log(1.0/l_min)
	return(exp(-alpha*lwave^(1/light_influence)))
}

#plot figure of the simulation result
#	x:time
#	y:foraging frequency of prey(blue) and predator(red))
#parameters
#	Ans: return value of tss_probforage_energygain_optimize function
#	title: title of the figure
#return
#	none
plot.sim_result = function(Ans,title,L){
	Prey= Ans$Prey
	Predator = Ans$Predator
	plot(0,0,type="n",
		  #		  xlab="time (t)",ylab="foraging predator (red), prey (blue)",
		  xlab="",ylab="",
		  xlim=c(0,tnum),ylim=c(-0.02,1.02),       
		  main=title
	)
	lines(t,L,col="gray",lwd=3)
	lines(t,Prey,col="blue",lwd=3,lty="dashed")
	lines(t,Predator*0.99,col="red",lwd=3)
}

#plot and save figures of the simulation result with performance wave
#	x:time
#	y in upper panel:predation performance and foraging efficiency of predator and prey
#	y in lower panel:foraging frequency of prey(blue) and predator(red))
#parameters
#	FigName: Figure name
#	Ans: return value of tss_probforage_energygain_optimize function
#	title: title of the figure
#return
#	none
plot_and_save.sim_result_with_wave = function(FigName,V,U,L,alpha,beta,mx,my,mb,phi,omega,h){
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
	
	png(paste(FigName,"_upper.png",sep=""),height=1200,width=1600)
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
	png(paste(FigName,"_lower.png",sep=""),height=1200,width=1600)
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
plot.assumption=function(t,watertemp,sharktemp,V,U){
	# assumption figure
	par(mfrow=c(1,2))
	plot(t,watertemp,col="blue",pch=15, xlab="time (t)",ylab="temperature")
	lines(t,watertemp,col="blue",lty=1)
	#	text(7,25.0,bquote('w'['t']))
	points(t,sharktemp,col="red",pch=16)
	lines(t,sharktemp,col="red",lty=1)
	points(t,calc.bodytemp(watertemp,1),col="purple",pch=17)
	lines(t,calc.bodytemp(watertemp,1),col="purple",lty=1)
	points(t,calc.bodytemp(watertemp,10),col="orange",pch=18)
	lines(t,calc.bodytemp(watertemp,10),col="orange",lty=1)
	#	text(4,28,bquote('s'['t']))
	#	par(new =T)
	#	plot(t,L,col="orange",pch=17, type = "p", axes = FALSE, ylab = "",ylim=c(0,max(L)),main="temperature")
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
	points(t,V-U,col="purple",pch=17)
	lines(t,V-U,col="purple",lty=1)
	points(t,(V-U)^3.0,col="green",pch=18)
	lines(t,(V-U)^3.0,col="green",lty=1)
	#	text(17,0.9,bquote(list('v'['t'],'- u'['t'])))
}

#major-active-time based categorization with 5 time zone
#parameters
#	Ans: return value of tss_probforage_energygain_optimize function
#	MajorProb: threshold for the definition of major-acitve time
#		major-active time := [foraging frequency during the focal duration]/[total foraging frequency] >= MajorProb
activetime6.get_category = function(Ans){
	ans = 0
	
	if(sum(Ans$Predator)!=0){
		{
			p = Ans$Predator
			
			group = c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4),rep(6,4))
			
			for(igroup in unique(group)){
				if(sum(p[group==igroup])>0){
					ans = ans + 10^(6-igroup)
					#					ans = ans + 2^(igroup-1)
				}
			}
		}
	}
	
	return(ans)
}

#transform majortime5 category to plotmode, which define the colour of plot
#parameters
#	catergory: return value of majortime.get_category
#return
#	plotmode (used in image.plotmode)
#		$clr: colour for image plot. "none" will be ignored.
#		$dotclr: colours for points. "none" will be ignored.
#		$err_category: category which is not assigned by any colour for image. 
activetime6.get_plotmode = function(category){
	#colour for categories which is not assigned by any colour 
	errclr = "white"
	
	#remove information of peak numbers
	majortime=(category)
	
	#definition of image colours
	clr = category
	clr[category==category] = errclr #error
	#add assignment of colours in the following lines
	clr[majortime== 000000]  = "black"	#nothing
	clr[majortime== 111111]  = "grey"	#asynchronous
	clr[majortime== 111110]  = "grey"	#asynchronous
	clr[majortime== 111101]  = "grey"	#asynchronous
	clr[majortime== 111011]  = "grey"	#asynchronous
	clr[majortime== 110111]  = "grey"	#asynchronous
	clr[majortime== 101111]  = "grey"	#asynchronous
	clr[majortime== 011111]  = "grey"	#asynchronous

	# only last  periods	
	clr[majortime== 000011]  = "blue"		#early nocturnal
	clr[majortime== 000001]  = "blue"		#early nocturnal
	# only first  periods	
	clr[majortime== 100000]  = "darkviolet"#late nocturnal
	clr[majortime== 110000]  = "darkviolet"#late nocturnal
	# both firs and last periods
	clr[majortime== 100001]  = "navyblue"	#nocturnal
	clr[majortime== 100011]  = "navyblue"	#nocturnal
	clr[majortime== 110001]  = "navyblue"	#nocturnal
	clr[majortime== 110011]  = "navyblue"	#nocturnal

	# only sunrise & morning
	clr[majortime== 001000]  = "darkpink"	#early diurnal
	clr[majortime== 010000]  = "darkpink"	#early diurnal
	clr[majortime== 011000]  = "darkpink"	#early diurnal
	# only afternoon & sunset
	clr[majortime== 000100]  = "orange"	#late diurnal
	clr[majortime== 000010]  = "orange"	#late diurnal
	clr[majortime== 000110]  = "orange"	#late diurnal
	# diurnal
	clr[majortime== 001100]  = "yellow"	#diurnal
	clr[majortime== 001110]  = "yellow"	#diurnal
	clr[majortime== 011100]  = "yellow"	#diurnal
	clr[majortime== 011110]  = "yellow"	#diurnal

	# all am
	clr[majortime== 111100]  = "skyblue"	#am
	clr[majortime== 111000]  = "skyblue"	#am
	clr[majortime== 111001]  = "skyblue"	#am
	# all pm
	clr[majortime== 001111]  = "plum"	#pm
	clr[majortime== 000111]  = "plum"	#pm
	clr[majortime== 100111]  = "plum"	#pm

	# 2nd + 5th, not after dark
	clr[majortime== 010010]  = "gold"	#crepuscular
	clr[majortime== 011010]  = "gold"	#crepuscular
	clr[majortime== 010110]  = "gold"	#crepuscular
	clr[majortime== 001010]  = "gold"	#crepuscular
	
	clr[majortime== 100010]  = "red"	#pre-sunrise crepuscular
	clr[majortime== 100100]  = "red"	#pre-sunrise crepuscular
	clr[majortime== 101000]  = "red"	#pre-sunrise crepuscular
	clr[majortime== 100110]  = "red"	#pre-sunrise crepuscular
	clr[majortime== 101100]  = "red"	#pre-sunrise crepuscular
	clr[majortime== 101110]  = "red"	#pre-sunrise crepuscular
	clr[majortime== 110010]  = "red"	#pre-sunrise crepuscular
	clr[majortime== 110100]  = "red"	#pre-sunrise crepuscular
	clr[majortime== 110110]  = "red"	#pre-sunrise crepuscular
	clr[majortime== 111010]  = "red"	#pre-sunrise crepuscular
	
	# crepuscular + after dark
	clr[majortime== 010001]  = "green"	#after sunset crepuscular
	clr[majortime== 001001]  = "green"	#after sunset crepuscular
	clr[majortime== 000101]  = "green"	#after sunset crepuscular
	clr[majortime== 011001]  = "green"	#after sunset crepuscular
	clr[majortime== 001101]  = "green"	#after sunset crepuscular
	clr[majortime== 011101]  = "green"	#after sunset crepuscular
	clr[majortime== 010011]  = "green"	#after sunset crepuscular
	clr[majortime== 001011]  = "green"	#after sunset crepuscular
	clr[majortime== 011011]  = "green"	#after sunset crepuscular
	clr[majortime== 010111]  = "green"	#after sunset crepuscular
	clr[majortime== 101001]  = "green"	#after sunset crepuscular
	clr[majortime== 100101]  = "green"	#after sunset crepuscular
	clr[majortime== 101101]  = "green"	#after sunset crepuscular
	clr[majortime== 110101]  = "green"	#after sunset crepuscular
	clr[majortime== 101011]  = "green"	#after sunset crepuscular
	
	#definition of dot colours
	dotclr = category
	dotclr[category==category] = "none"
	
	#list of error category (not assigned by any colour )
	err = sort(unique(as.vector(category[clr == errclr])))
	
	return(list(clr=clr,dotclr=dotclr, err_category=err))
}

#plot category with different colours following the definition in plotmode
#	x.seq: sequence of xaxis
#	y.seq: sequence of yaxis
#	plotmode: matrix of plotmode (return value of majortime.get_plotmode or majortime5.get_plotmode)
#	dotcex: change the scale of dots 
image.plotmode=function(x.seq,y.seq,plotmode,dotcex = 0.33,plot_legend=FALSE,...){
	
	#colstrat=c("nothing","asynchronous","early_nocturnal","late_nocturnal","long_nocturnal","sunrise","sunset","diurnal","afternoon","pm","crepuscular","post-sunrise crepuscular","post-sunset crepuscular")
	#collabs=c("black","grey","blue","darkviolet","navyblue","deeppink","pink","yellow","orange","plum","gold","red","seagreen2")
	
	fz = factor(plotmode$clr)
	z = matrix(as.integer(fz),length(x.seq),length(y.seq))
	clr = levels(fz)
	clr[clr=="none"] = rgb(0,0,0,0)
	image(x.seq,y.seq,z,col=levels(fz),...)
	if(plot_legend)legend("topright",legend=colstrat, pch=16, col=collabs, cex=0.8, pt.cex = 2.2)
	
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
