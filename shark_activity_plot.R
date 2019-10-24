#load library
source("shark_activity_functions.R")

plot.r.phi.figures=function(nameIn,vk,beta,l0,mx,my,v0,r,phi,watertemp,wmax,wmin,u0,uk,alpha,predcost,predCostMass,L,omega,h,mb,l_min,light_influence,twilight_coef,t){
  x.seq = seq(0.05,1.00,length=4)
  y.seq = seq(0.02,0.18,length=5)
  name = paste(nameIn,"_r-phi","_B",beta,"_uk",10*uk,"_vK",10*vk,"_L",10*l0,"_mX",10*mx,"_mY",10*my,"_v0",v0*10,"_r",r,"_c",predCostMass*100,"_phi",phi*100,sep = "")
  png(paste("examples_",name,".png",sep=""),height=2000,width=2000)
  par(mfrow=c(length(y.seq),length(x.seq)),cex=2.0,mex=0.3)
  for(phi in rev(y.seq)){
    for(r in x.seq){
      # cost over the day
	    C_day = predcost * (1 + predCostMass*(r^3)*(38-watertemp))
      sharktemp=calc.sharktemp(t,t_w,wmin,wmax,r) 
      U = u0 + uk*(watertemp-(wmax+wmin)/2)
      V = v0 + vk*(sharktemp-(wmax+wmin)/2)
      L = calc.light_effect(t,l_min, light_influence, twilight_coef)
      #run simulation
      Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C_day, L, my, phi, omega, beta, h, mb,mx)
      #plot simulation results
      plot.sim_result(Ans,bquote(list("r"==.(r),'phi'==.(phi))),L)
    }
  }
  print(Ans$Predator)
  dev.off()
  grid = 101
  x.ax = seq(0.01,1.00,length=grid)
  y.ax = seq(0.0,0.2,length=grid)
  #plot categories of simulation results with changing v0 and phi
  no = matrix(0,grid,grid)
  for(y in 1:length(y.ax)){
    phi = y.ax[y]
    for(x in 1:length(x.ax)){
      r = x.ax[x]
      # cost over the day
      C_day = predcost * (1 + predCostMass*(r^3)*(38-watertemp))
      sharktemp=calc.sharktemp(t,t_w,wmin,wmax,r) 
      U = u0 + uk*(watertemp-(wmax+wmin)/2)
      V = v0 + vk*(sharktemp-(wmax+wmin)/2)
      L = calc.light_effect(t,l_min, light_influence, twilight_coef)
      #run simulation
      Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C_day, L, my, phi, omega, beta, h, mb,mx)
      #calculate category
      no[x,y] = activetime6.get_category(Ans)
    }
  }
  #plot category image
  plotmode = activetime6.get_plotmode(no)
  png(paste("zone_",name,".png",sep=""),height=1600,width=1600)
  par(mfrow=c(1,1),cex=4.0,bg=rgb(0,0,0,0))
  image.plotmode(x.ax,y.ax,plotmode,xlab="r",ylab="phi")
  dev.off()
  #list of categorization error (grey colors) 
  plotmode$err_category
  #list of categorization error (grey colors) by ignoring the number of peaks
  sort(unique((plotmode$err_category)%%100))
}

plot.my.v0.figures=function(nameIn,vk,beta,l0,mx,my,v0,r,phi,watertemp,wmax,wmin,u0,uk,alpha,predcost,predCostMass,L,omega,h,mb,l_min,light_influence,twilight_coef,t){
  x.seq = seq(1.0,2.0,length=5)
  y.seq = seq(0.0,2.0,length=5)
  name = paste(nameIn,"_v0-my","_B",beta,"_uk",10*uk,"_vK",10*vk,"_L",10*l0,"_mX",10*mx,"_mY",10*my,"_v0",v0*10,"_r",r,"_c",predCostMass*100,"_phi",phi*100,sep = "")
  #plot multiple results of simulations with changing v0 and r
  png(paste("examples_",name,".png",sep=""),height=2000,width=2000)
  par(mfrow=c(length(y.seq),length(x.seq)),cex=2.0,mex=0.3)
  for(my in rev(y.seq)){
    for(v0 in x.seq){
      # cost over the day
      C_day = predcost * (1 + predCostMass*(r^3)*(38-watertemp))
      sharktemp=calc.sharktemp(t,t_w,wmin,wmax,r) 
      U = u0 + uk*(watertemp-(wmax+wmin)/2)
      V = v0 + vk*(sharktemp-(wmax+wmin)/2)
      L = calc.light_effect(t,l_min, light_influence, twilight_coef)
      #run simulation
      Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C_day, L, my, phi, omega, beta, h, mb,mx)
      #plot simulation results
      plot.sim_result(Ans,bquote(list('v'['0']==.(v0),'m'['y']==.(my))),L)
    }
  }
  dev.off()
  grid = 101
  x.ax = seq(1.0,2.0,length=grid)
  y.ax = seq(0.0,2.0,length=grid)
  #y.ax  = seq(1.0,2.0,length=grid)
  #plot categories of simulation results with changing v0 and r
  no = matrix(0,grid,grid)
  for(y in 1:length(y.ax)){
    my = y.ax[y]
    for(x in 1:length(x.ax)){
      v0 = x.ax[x]
      # cost over the day
      C_day = predcost * (1 + predCostMass*(r^3)*(38-watertemp))
      sharktemp=calc.sharktemp(t,t_w,wmin,wmax,r) 
      U = u0 + uk*(watertemp-(wmax+wmin)/2)
      V = v0 + vk*(sharktemp-(wmax+wmin)/2)
      L = calc.light_effect(t,l_min, light_influence, twilight_coef)
      #run simulation
      Ans = tss_probforage_energygain_optimize_linear(V, U, alpha, C_day, L, my, phi, omega, beta, h, mb,mx)
      #calculate category
      no[x,y] = activetime6.get_category(Ans)
    }
  }
  #plot category image
  plotmode = activetime6.get_plotmode(no)
  png(paste("zone_",name,".png",sep=""),height=1600,width=1600)
  par(mfrow=c(1,1),cex=4.0,bg=rgb(0,0,0,0))
  image.plotmode(x.ax,y.ax,plotmode,xlab="v0",ylab="mY")
  dev.off()
  #list of categorization error (grey colors) 
  plotmode$err_category
  #list of categorization error (grey colors) by ignoring the number of peaks
  sort(unique((plotmode$err_category)%%100))
}

#=== constant parameters ===
#set time
pi = acos(-1)
tnum = 24 # time of day
t = 1:tnum-0.5

# temperature of the water
t_w=15
wmin = 25
wmax = 30
watertemp = calc.watertemp(t,t_w,wmin,wmax)

# amount of food availability for prey
alpha = rep(1.0, length=tnum)
# baseline mortality for prey (should pay both for resting and foraging)
mb = 0.2	
# average prey swim speed
u0 = 1.0
# metabolic cost for predators when they go out for predation
predcost=0.15
predCostMass=0.0 # due to endothermy

#=== default parameter values ===
omega = 1.0	#foraging efficiency increment by speed 
beta = 1.0 	#predation efficiency
phi = 0.1  	#probability of failing to hide in safe place
h = 1.0    	#handling time

#prey and predator speed
r = 0.3  #body radius (meter)
uk = 0.2 #influence of bodytemp
v0 = 1.5  #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp

#light effects
l0 = 1.0				#REQUIRE non-zero value! predation efficiency at noon
l1 = 1.0					#predation efficiency in perfect dark
light_influence = 5.0			#determines the sensitivity for small light
l_min=1;
twilight_coef=0.3;
#if light_influence is small, predation efficiency under twilight is mostly same with that in perfect dark.
#if light_influence is large, predation efficiency under twilight is mostly same with that at noon

#mortality rate of prey by predation
mx = 1.0	#predation by other predators
my = 0.5 #predation by sharks

sharktemp=calc.sharktemp(t,t_w,wmin,wmax,r) 
U = u0 + uk*(watertemp-(wmax+wmin)/2)
V = v0 + vk*(sharktemp-(wmax+wmin)/2)
L = calc.light_effect(t,l_min, light_influence, twilight_coef)

# FIGURE 5 & A2  basuc result
plot.my.v0.figures("fig5",vk,beta,l0,mx,my,v0,r,phi,watertemp,wmax,wmin,u0,uk,alpha,predcost,predCostMass,L,omega,h,mb,l_min,light_influence,twilight_coef,t)

# FIGURE A3 simple models
for(l0 in c(1.0,0.2)){# effect of light
  for(mx in c(0.0,1.0)){ # other predators
    for(beta in c(0.0,1.0)){ # effect of speed
    for(uk in c(0.0,0.2)){# effect of temp on speed
      vk=uk
      plot.my.v0.figures("figA3",vk,beta,l0,mx,my,v0,r,phi,watertemp,wmax,wmin,u0,uk,alpha,predcost,predCostMass,L,omega,h,mb,l_min,light_influence,twilight_coef,t)
      }
    }
  }
}

beta=1.0
# FIGURE 6 A4 A5 
for(l0 in c(0.1,0.5,1.0)){
  for(mx in c(0.0,0.5,1.0,2.0)){
    for(my in c(0.0,0.5,1.0)){
      #=== draw r-phi figures === 
      plot.r.phi.figures("fig6",vk,beta,l0,mx,my,v0,r,phi,watertemp,wmax,wmin,u0,uk,alpha,predcost,predCostMass,L,omega,h,mb,l_min,light_influence,twilight_coef,t)
    }
  }
}



  