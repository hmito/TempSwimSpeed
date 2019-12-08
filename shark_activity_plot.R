#load library
source("shark_activity_functions.R")

plot.r.phi.figures=function(nameIn, t, tw, wmin, wmax, ub, uk, vb, vk, lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h, plot_legend ){
  x.seq = seq(0.05,1.00,length=4)
  y.seq = seq(0.02,0.24,length=5)
  name = paste(nameIn,"_r-phi","_B",beta,"_uk",10*uk,"_vk",10*vk,"_lm",10*lm,"_mX",10*mx,"_mY",10*my,"_vb",vb*10,"_omega",omega*10,"_c",cost*100,"_phi",100*phi,sep = "")
  png(paste("examples_",name,".png",sep=""),height=2000,width=2000)
  par(mfrow=c(length(y.seq),length(x.seq)),cex=2.0,mex=0.3)
  for(phi in rev(y.seq)){
    for(r in x.seq){
      # calculate time-depending parameters 
      watertemp=calc.watertemp(t,tw,wmin,wmax)
      sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
      U = ub + uk*(watertemp-(wmax+wmin)/2)
      V = vb + vk*(sharktemp-(wmax+wmin)/2)
      K = rep(alpha,length=length(t))
	    C = rep(cost,length=length(t))
      L = calc.light_effect(t, lm, lk, ld)
      #run simulation
      Ans = tss_probforage_energygain_optimize_linear(V, U, K, C, L, my, phi, omega, beta, h, mb,mx)
      #plot simulation results
      plot.sim_result(Ans,bquote(list("r"==.(r),'phi'==.(phi))),L)
    }
  }
  dev.off()

  
  #plot categories of simulation results with changing v0 and phi
  grid = 101
  x.ax = seq(0.01,1.00,length=grid)
  y.ax = seq(0.0,0.25,length=grid)
  no = matrix(0,grid,grid)
  for(y in 1:length(y.ax)){
    phi = y.ax[y]
    for(x in 1:length(x.ax)){
      r = x.ax[x]
      # calculate time-depending parameters 
      watertemp=calc.watertemp(t,tw,wmin,wmax)
      sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
      U = ub + uk*(watertemp-(wmax+wmin)/2)
      V = vb + vk*(sharktemp-(wmax+wmin)/2)
      K = rep(alpha,length=length(t))
	    C = rep(cost,length=length(t))
      L = calc.light_effect(t, lm, lk, ld)
      #run simulation
      Ans = tss_probforage_energygain_optimize_linear(V, U, K, C, L, my, phi, omega, beta, h, mb,mx)
      #calculate category
      no[x,y] = activetime6.get_category(Ans)
    }
  }
  #plot category image
  plotmode = activetime6.get_plotmode(no)
  png(paste("zone_",name,".png",sep=""),height=1600,width=1600)
  par(mfrow=c(1,1),cex=5.0,bg=rgb(0,0,0,0))
  image.plotmode(x.ax,y.ax,plotmode,xlab="r",ylab="phi",plot_legend = plot_legend)
  dev.off()
  
  #list of categorization error (grey colors) 
  plotmode$err_category
  #list of categorization error (grey colors) by ignoring the number of peaks
  sort(unique((plotmode$err_category)%%100))
}

plot.vb.my.figures=function(nameIn, t, tw, wmin, wmax, ub, uk, vb, vk, lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h, plot_legend ){
  x.seq = seq(1.0,2.0,length=5)
  y.seq = seq(0.0,2.0,length=5)
  name = paste(nameIn,"_vb-my","_B",beta,"_uk",10*uk,"_vK",10*vk,"_lm",10*lm,"_mX",10*mx,"_mY","[y]","_vb",vb*10,"_omega",omega*10,"_c",cost*100,"_phi",phi*100,sep = "")

  #plot multiple results of simulations with changing v0 and r
  png(paste("examples_",name,".png",sep=""),height=2000,width=2000)
  par(mfrow=c(length(y.seq),length(x.seq)),cex=2.0,mex=0.3)
  for(my in rev(y.seq)){
    for(vb in x.seq){
    	# calculate time-depending parameters
      watertemp=calc.watertemp(t,tw,wmin,wmax)
      sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
      U = ub + uk*(watertemp-(wmax+wmin)/2)
      V = vb + vk*(sharktemp-(wmax+wmin)/2)
      K = rep(alpha,length=length(t))
	   C = rep(cost,length=length(t))
      L = calc.light_effect(t, lm, lk, ld)
      #run simulation
      Ans = tss_probforage_energygain_optimize_linear(V, U, K, C, L, my, phi, omega, beta, h, mb,mx)
      #plot simulation results
      plot.sim_result(Ans,bquote(list('v'['b']==.(vb),'m'['y']==.(my))),L)
    }
  }
  dev.off()

  
  #plot categories of simulation results with changing v0 and r
  grid = 101
  x.ax = seq(1.0,2.0,length=grid)
  y.ax = seq(0.0,2.0,length=grid)
  no = matrix(0,grid,grid)
  for(y in 1:length(y.ax)){
    my = y.ax[y]
    for(x in 1:length(x.ax)){
      vb = x.ax[x]

      # calculate time-depending parameters       
      watertemp=calc.watertemp(t,tw,wmin,wmax)
      sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
      U = ub + uk*(watertemp-(wmax+wmin)/2)
      V = vb + vk*(sharktemp-(wmax+wmin)/2)
      K = rep(alpha,length=length(t))
	   C = rep(cost,length=length(t))
      L = calc.light_effect(t, lm, lk, ld)
      #run simulation
      Ans = tss_probforage_energygain_optimize_linear(V, U, K, C, L, my, phi, omega, beta, h, mb,mx)
      #calculate category
      no[x,y] = activetime6.get_category(Ans)
    }
  }
  #plot category image
  plotmode = activetime6.get_plotmode(no)
  png(paste("zone_",name,".png",sep=""),height=1600,width=1600)
  par(mfrow=c(1,1),cex=5.0,bg=rgb(0,0,0,0))
  image.plotmode(x.ax,y.ax,plotmode,xlab="vb",ylab="mY",plot_legend = plot_legend)
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
tw=15
wmin = 25
wmax = 30

#prey and predator speed
ub = 1.0 #[can be fixed] average swim speed
uk = 0.2 #influence of bodytemp
vb = 1.5 #average swim speed (prey is always 1.0)
vk = 0.2 #influence of bodytemp

#light effects
lm = 0.5		#predation efficiency at noon
lk = 0.2		#determines the sensitivity for small light
ld = 0.3		#duration of twilight

#mortality rate of prey by predation
alpha = 1.0 #amount of food for prey
omega = 1.0	#foraging efficiency increment by speed 
phi = 0.1  	#probability of failing to hide in safe place
mb = 0.2	#baseline mortality rate
mx = 1.0	#predation by other predators
my = 0.5 #predation by sharks

#=== default parameter values ===
r = 0.3  	#predator's body radius (meter)
cost=0.15	#predation cost
beta = 1.0 	#predation efficiency
h = 1.0    	#handling time


# FIGURE 4 basic zoneplot
plot.vb.my.figures("fig4", t, tw, wmin, wmax, ub, uk, vb, vk, 
                   lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend = TRUE)
	
# FIGURE A2 simple models
for(lm in c(0.0,4.0)){# effect of light
  for(mx in c(0.0,1.0)){ # other predators
    for(beta in c(0.0,1.0)){ # effect of speed
    for(uk in c(0.0,0.2)){# effect of temp on speed
      vk=uk
      plot_legend = FALSE
      if (lm==0.0 && mx==0.0 && beta==0.0 && uk==0.0) { 
        plot_legend = TRUE 
      }
      plot.vb.my.figures("figA2", t, tw, wmin, wmax, ub, uk, vb, vk, 
                         lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend )
      }
    }
  }
}

beta=1.0
mx = 1.0
# FIGURE 5 A2 A3
for(lm in c(4.0,1.0,0.0)){
  for(omega in c(0.0,1.0,2.0)){
    for(mx in c(0.0,1.0,2.0)){
      for(my in c(0.0,1.0,2.0)){
        plot_legend = FALSE
        if (mx==0.0 & my==0.0) { 
          plot_legend = TRUE 
        }
      	plot.r.phi.figures("figA3",t, tw, wmin, wmax, ub, uk, vb, vk, 
      	                   lm, lk, ld, alpha, omega, phi, mb, mx, my, r, cost, beta, h,plot_legend )
      }
    }
  }
}



  
