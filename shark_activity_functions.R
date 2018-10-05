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
calc.light_effect_old = function(t,dark,light){
	lwave=exp(-1.0*cos(2*pi*(length(t)/2-t)/length(t)))/exp(1.0)
	return(dark+(light-dark)*lwave)
}
calc.light_effect = function(t, l_min, l_max, light_influence, twilight_coef=0.3){
	#Around twilight_coef = 0.3 seems to be reasonable because
	#	- the definition of astronominal twilight is that the sun is less than -18 degree below the horizon. 
	#	- On the Equator, 18 degree passes for 1.2 hours, so twilight start from t = 4.8 or end at t=19.2.
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

#peak based pattern categorization
allpeak.get_category = function(Ans, ThrPrb=0.20){
	ans = 0
	thr = 0.5
	
	p = Ans$Predator
#	p = (p + c(p[-1],p[1]) + c(p[length(p)],p[-length(p)]))/3
	p[c(p[-1],p[1])>thr&c(p[length(p)],p[-length(p)]) > thr]= 1
		
	if(all(p<=thr)){
		return(0)
	}else if(all(p>thr)){
		return(44)
	}
	
	p[p>0.05]=1
	peaks = hist.find_peaks(p,0,thr)
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
	#2:before-noon(4-11)
	#3:after-noon(12-19)
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
allpeak.get_plotmode = function(category){
	no = category%%100
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
	
	sub = (category-no)/100
	submode = sub
	submode[sub==submode]=0	#reset all
	submode[sub==22] =1 #green: beforenoon
	submode[sub==12] =2 #green4: midnight-beforenoon
	submode[sub>100] =3 #black: more than one sub-peaks
	
	return(list(main=mode,sub=submode))
}
detail.allpeak.no.color=function(){
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
allpeak.point=function(x.seq,y.seq,plotmode,...){
	submode = plotmode$sub
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
allpeak.image=function(x.seq,y.seq,plotmode,...){
	main = plotmode$main
	col = detail.allpeak.no.color()
	image(x.seq,y.seq,main,zlim=c(0,length(col)-1),col=col,...)
}

#major-active-time based categorization
majortime.get_category2 = function(Ans, MajorProb = 0.25){
	thr = 0.5
	
	p = Ans$Predator
	#	p = (p + c(p[-1],p[1]) + c(p[length(p)],p[-length(p)]))/3
	p[c(p[-1],p[1])>thr&c(p[length(p)],p[-length(p)]) > thr]= 1
	
	if(all(p<=thr)){
		return(0)
	}
	
	p[p>0.05]=1
	peaks = hist.find_peaks(p,0,thr)
	#connect 24-1
	if(nrow(peaks)>1 && peaks$lower[1]==1 && peaks$upper[nrow(peaks)]==24){
		peaks$lower[1] = peaks$lower[nrow(peaks)]
		peaks$freq[1] = peaks$freq[1] + peaks$freq[nrow(peaks)]
		peaks = peaks[-nrow(peaks),]
	}
	
	peaks = peaks[peaks$freq/sum(peaks$freq)>0.00,]
	peaks = peaks[order(peaks$freq),]
	
	#--time division--
	timediv = c(rep(1,length=6),rep(2,length=6),rep(3,length=6),rep(4,length=6))
	#timediv = c(rep(1,3),rep(2,8),rep(3,8),rep(1,5))
	ans = 0
	for(i in 1:nrow(peaks)){
		subp = numeric(length(p))
		if(peaks$lower[i] <= peaks$upper[i]){
			subp[peaks$lower[i]:peaks$upper[i]] = p[peaks$lower[i]:peaks$upper[i]]
		}else{
			subp[ 1:peaks$upper[i]] = p[ 1:peaks$upper[i]]
			subp[peaks$lower[i]:24] = p[peaks$lower[i]:24]
		}
		for(j in unique(timediv)){
			if(sum(subp[timediv==j])/sum(subp) >= MajorProb){
				ans = ans + 2^(j-1)*100^(i-1)
			}
		}
	}
	return(ans)
}
majortime.get_plotmode2 = function(category){
	no = category
	no[category==category] = 0
	no[category== 0] = 1	#nothing
	no[category==15] = 2	#all
	no[category== 1] = 3	#night-am
	no[category== 8] = 4	#night-pm
	no[category== 9] = 5	#night
	no[category== 2] = 6	#day-am
	no[category== 4] = 7	#day-pm
	no[category== 6] = 8	#day
	no[category== 3] = 9	#pm
	no[category==12] =10	#am
	no[category==11] =11	#night + day-am
	no[category==14] =12	#day + night-pm
	no[category>=100]=13	#two peaks
	
	return(no)
}


#major-active-time based categorization
majortime.get_category3 = function(Ans, MajorProb = 0.25){
	thr = 0.5
	
	p = Ans$Predator
	p[c(p[-1],p[1])>thr&c(p[length(p)],p[-length(p)]) > thr]= 1
	
	if(sum(p>=thr)==0){
		return(0)
	}else if(sum(p>=thr)>=18){
		return(1)
	}
	
	p[p>0.05]=1
	peaks = hist.find_peaks(p,0,thr)
	#connect 24-1
	if(nrow(peaks)>1 && peaks$lower[1]==1 && peaks$upper[nrow(peaks)]==24){
		peaks$lower[1] = peaks$lower[nrow(peaks)]
		peaks$freq[1] = peaks$freq[1] + peaks$freq[nrow(peaks)]
		peaks = peaks[-nrow(peaks),]
	}
	
	peaks = peaks[peaks$freq/sum(peaks$freq)>0.00,]
	peaks = peaks[order(peaks$freq,decreasing=TRUE),]
	
	#--time division--
	timediv = c(rep(1,length=3),rep(2,length=6),rep(3,length=6),rep(4,length=6),rep(1,length=3))
	#timediv = c(rep(1,3),rep(2,8),rep(3,8),rep(1,5))
	ans = 0
	for(i in 1:nrow(peaks)){
		if(peaks$lower[i] <= peaks$upper[i]){
			mean = mean(peaks$lower[i]:peaks$upper[i])
		}else{
			mean = mean(c(peaks$lower[i]:24, 1:peaks$upper[i]+24))%%24
		}
		
		ans = ans + (
			2*(1*(3<=mean & mean<9)+
				2*(9<=mean & mean<15)+
				3*(15<=mean& mean<21)+
				4*(21<=mean | mean<3)
			)+1*(peaks$freq[i]>=12)
		)*10^(i-1)
	}
	return(ans)
}
majortime.get_plotmode3 = function(category){
	no = category
	no[category==category] = 0
	no[category== 0] = 1	#nothing
	no[category== 1] = 2	#all
	no[category== 2] = 3	#morning short
	no[category== 3] = 4	#morning long
	no[category== 4] = 5	#day short
	no[category== 5] = 6	#day long
	no[category== 6] = 7	#evening short
	no[category== 7] = 8	#evening long
	no[category== 8] = 9	#night short
	no[category== 9] =10	#night long
	no[category >10] =11	#two peaks
	
	return(no)
}
detail.majortime.no.color3 = function(){
	return(c("yellow","black","white","green", "forestgreen","pink","red","yellow","gold","skyblue","blue","grey"))
}

#major-active-time based categorization
majortime.get_category = function(Ans, MajorProb = 0.20){
	ans = 0
	#p = rep(c(1,0,1),times=c(5,10,9))
	
	if(sum(Ans$Predator)!=0){
		#major active time 
		{
			p = Ans$Predator
			group = rep(1:4,times=c(6,6,6,6))
#			group = rep(c(1:4,1),times=c(3,6,6,6,3))
			
			for(igroup in unique(group)){
				if(sum(p[group==igroup])/sum(p) > MajorProb){
					ans = ans + 2^(igroup-1)
				}
			}
		}
		
		#peak number
		{
			thr = 0.5
		#	p[c(p[-1],p[1])>thr&c(p[length(p)],p[-length(p)]) > thr]= 1
		#	p[c(p[-1],p[1])<thr&c(p[length(p)],p[-length(p)]) < thr]= 0
			
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
majortime.get_plotmode = function(category){
	no = category
	no[category==category] = 0
	no[category== 0]  = 1	#nothing
	no[category== 15] = 2	#all
	no[category== 6]  = 3	#day
	no[category== 8]  = 4	#early-night
	no[category== 12] = 5	#pm
	no[category== 11] = 6	#except pm-day
	no[category== 14] = 7	#except am-night
	no[category== 9]  = 8	#night
	no[category== 13] = 9	#except am-day
	no[category== 7] = 10	#except pm-night
	no[category== 1] = 11	#am-night
	no[category== 4] = 12	#pm-day
	no[category>=16 & category<32]  = 13	#two peaks
	no[category>=32]  = 14					#three peaks
	return(no)
}
majortime.image=function(x.seq,y.seq,plotmode,legend=FALSE,...){
	#black	: nothing
	#white	: all time use
	#red		: day time (6-18)
	#purple	: early-night (18-24)
	#yellow	: pm (12-24)
	#skyblue	: except pm-day (18-12)
	#orange	: except am-night (6-24)
	#blue		: night time (18-6)
	#pale purple: except am-day (12-6)
	#brown	: except pm-night (0-18)
	#darkblue: am-night
	#darkred	: pm-day
	col=c("grey","black","white","red","purple","yellow","skyblue","orange","blue","mediumpurple1","saddlebrown","blue4","red3","forestgreen","green")
	image(x.seq,y.seq,plotmode,zlim=c(0,length(col)-1),col=col,...)
	if(legend)legend("topright", pch=19,legend = c("error","no","all","day","erl-night","pm","ex pm-day","ex am-night","night", "ex am-day","ex pm-night", "am-night","pm-day","two peaks","three peaks"), col = col)
}
majortime.image_black=function(x.seq,y.seq,plotmode,...){
	col=c("white","black","white","white","white","white","white","white","white","darkgrey","grey")
	image(x.seq,y.seq,plotmode,zlim=c(0,length(col)-1),col=col,...)
	for(no in 3:8){
		hmRLib::image_polygon(x.seq,y.seq,plotmode==no,...)
	}
}

majortime.get_category4 = function(Ans, MajorProb = 0.20){
	ans = 0
	#p = rep(c(1,0,1),times=c(5,10,9))
	
	if(sum(Ans$Predator)!=0){
		thr = 0.5
		#	p[c(p[-1],p[1])>thr&c(p[length(p)],p[-length(p)]) > thr]= 1
		#	p[c(p[-1],p[1])<thr&c(p[length(p)],p[-length(p)]) < thr]= 0
		
		peaks = hist.find_peaks(Ans$Predator,0,thr)
		
		#connect 24-1
		if(nrow(peaks)>1 && peaks$lower[1]==1 && peaks$upper[nrow(peaks)]==24){
			peaks$lower[1] = peaks$lower[nrow(peaks)]
			peaks$freq[1] = peaks$freq[1] + peaks$freq[nrow(peaks)]
			peaks = peaks[-nrow(peaks),]
		}
		peaks = peaks[order(peaks$freq,decreasing = TRUE),]
		
		for(ipeak in 1:nrow(peaks)){
			time = logical(24)
			if(peaks$lower[ipeak] <= peaks$upper[ipeak]){
				time[peaks$lower[ipeak]:peaks$upper[ipeak]]=TRUE
			}else{
				time[peaks$lower[ipeak]:24]=TRUE
				time[ 1:peaks$upper[ipeak]]=TRUE
			}
			p = Ans$Predator
			p[!time] = 0

			group = rep(1:4,times=c(6,6,6,6))
			#group = rep(c(1:4,1),times=c(3,6,6,6,3))
			
			for(igroup in unique(group)){
				if(sum(p[group==igroup])/sum(p) > MajorProb){
					ans = ans + 2^(igroup-1)*100^(ipeak-1)
				}
			}
		}
	}
	
	return(ans)
}


detail.majortime.no.color3 = function(){
	return(c("yellow","black","white","green", "forestgreen","pink","red","yellow","gold","skyblue","blue","grey"))
}


