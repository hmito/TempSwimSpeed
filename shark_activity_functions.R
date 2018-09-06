#load library
library("Rcpp")
library("BH")
sourceCpp("TempSwimSpeedOptim.cpp")

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

#light effect
calc.light_effect = function(t,lmin,lmax){
	lwave=exp(-1.0*cos(2*pi*(length(t)/2-t)/length(t)))/exp(1.0)
	return(lmin+(lmax-lmin)*lwave)
}

#Summarized figure
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

#plot pair of temp and V,U,L
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
	
	plot(t,V,col="red",pch=16,xlim=c(0,tnum),ylim=c(min(V-U),max(c(V,U))),
		  xlab="time (t)",ylab="burst speed",main="swim speed")
	segments(-100,0,100,0)
	lines(t,V,col="red",lty=1)
	#	text(15,2.3,bquote('v'['t']))
	points(t,U,col="blue",pch=15)
	lines(t,U,col="blue",lty=1)
	#	text(18,1.9,bquote('u'['t']))
	points(t,V-U,col="purple",pch=2)
	lines(t,V-U,col="purple",lty=1)
	#	text(17,0.9,bquote(list('v'['t'],'- u'['t'])))
}

#peak-based pattern categorization
allpeak_no.sim_result=function(Ans,ThrPrb=0.20){
	ans = 0
	
	p = Ans$Predator
	p = (p + c(p[-1],p[1]) + c(p[length(p)],p[-length(p)]))/3
	
	if(all(p<=0.25)){
		return(0)
	}else if(all(p>0.25)){
		return(44)
	}
	
	peaks = hist.find_peaks(p,0,0.05)
	#connect 24-1
	if(nrow(peaks)>1 && peaks$lower[1]==1 && peaks$upper[nrow(peaks)]==24){
		peaks$lower[1] = peaks$lower[nrow(peaks)]
		peaks$freq[1] = peaks$freq[1] + peaks$freq[nrow(peaks)]
		peaks = peaks[-nrow(peaks),]
	}
	
	peaks = peaks[peaks$freq/sum(peaks$freq)>ThrPrb,]
	peaks = peaks[order(peaks$freq),]
	
	#--time division--
	#1:night (20-3)
	#2:before-noon(4-12)
	#3:after-noon(12-19)
	#evening:4:16-21
	#timediv = c(1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,1,1,1)
	timediv = c(rep(1,3),rep(2,8),rep(3,8),rep(1,5))
	ans = 0
	for(i in 1:nrow(peaks)){
		ans = ans*100 + 
			timediv[peaks$lower[i]]*10 + 
			ifelse(peaks$freq[i]>18 && timediv[peaks$lower[i]]==timediv[peaks$upper[i]],4,timediv[peaks$upper[i]])
	}
	
	return(ans)
}
allpeak_no.mode = function(original_no){
	no = original_no%%100
	mode = no
	mode[no==mode]=0	#reset all
	mode[no==0] =1 #black	: no use
	mode[no==44]=2 #white	: all use
	mode[no==32]=3 #cyan		: from afternoon  to beforenoon (rest around noon)
	mode[no==31]=4 #blue		: from afternoon  to mid-night
	mode[no==21]=5 #red		: from beforenoon to mid-night (rest early morning)
	mode[no==22]=6 #darkgreen	: beforenoon
	mode[no==33]=7 #darkred	: afternoon
	mode[no==11]=8 #darkblue: midnight
	mode[no==24]=9 #purple	: late beforenoon to early beforenoon (rest befornoon)
	mode[no==12]=10 #green: midnight to beforenoon (rest befornoon)
	mode[no==34]=11 #yellowgreen	: late afternoon to early afternoon (rest afternoon)
	mode[no==23]=12 #yellow	: beforenoon to afternoon (daytime)
	mode[no==14]=13 #orange : late midnight to early midnight (rest midnight)
	
	sub = (original_no-no)/100
	submode = sub
	submode[sub==submode]=0	#reset all
	submode[sub==22] =1 #green: beforenoon
	submode[sub==12] =2 #green4: midnight-beforenoon
	submode[sub>100] =3 #black: more than one sub-peaks
	
	return(list(main=mode,sub=submode))
}
allpeak_no.color=function(){
	return(
		c("grey",
		  "black",
		  "white",
		  "cyan",
		  "blue",
		  "red",
		  "green4",
		  "red4",
		  "darkblue",
		  "purple",
		  "green",
		  "yellowgreen",
		  "yellow",
		  "orange"
		)
	)
}
allpeak_no.point=function(x.seq,y.seq,mode,...){
	submode = mode$sub
	x = matrix(rep(x.seq,times=length(y.seq)),length(x.seq),length(y.seq))
	y = matrix(rep(y.seq,each =length(x.seq)),length(x.seq),length(y.seq))
	
	x.pos = as.vector(x[submode==1])
	y.pos = as.vector(y[submode==1])
	points(x.pos,y.pos,col="green",pch=16,...)
	
	x.pos = as.vector(x[submode==2])
	y.pos = as.vector(y[submode==2])
	points(x.pos,y.pos,col="green4",pch=16,...)
	
	x.pos = as.vector(x[submode==3])
	y.pos = as.vector(y[submode==3])
	points(x.pos,y.pos,col="black",pch=17,...)
}
allpeak_no.image=function(x.seq,y.seq,mode,...){
	main = mode$main
	col = allpeak_no.color()
	image(x.seq,y.seq,main,zlim=c(0,length(col)-1),col=col,...)
}

