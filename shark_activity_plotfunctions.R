#load library
source("shark_activity_functions.R")

plot.r.phi.figures=function(nameIn, t, tw, wmin, wmax, ub, uk, vb, vk, mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h, plot_legend, plot_example=FALSE){
	x.seq = seq(0.2,5.00,length=4)
	y.seq = seq(0.02,0.24,length=5)
	name = paste(nameIn,"_r-phi","_B",beta,"_uk",10*uk,"_vk",10*vk,"_lm",10*lm,"_mX",10*mx,"_mY",10*my,"_vb",vb*10,"_omega",omega*10,"_r","[x]","_c",cost*100,"_phi","[y]",sep = "")
	if(plot_example){
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
				L = calc.light_effect(t, mu,rho,kappa,sigma)
				#run simulation
				Ans = tss_probforage_energygain_optimize_linear(V, U, K, C, L, my, phi, omega, beta, h, mb,mx)
				#plot simulation results
				plot.sim_result(Ans,bquote(list("r"==.(r),'phi'==.(phi))),L)
			}
		}
		dev.off()
	}
	
	#plot categories of simulation results with changing v0 and phi
	grid = 101
	x.ax = seq(0.05,6.00,length=grid)
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
			L = calc.light_effect(t, mu,rho,kappa,sigma)
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
	image.plotmode(x.ax,y.ax,plotmode,xlab="",ylab="",plot_legend = plot_legend)
	dev.off()
	
	#list of categorization error (grey colors) 
	plotmode$err_category
	#list of categorization error (grey colors) by ignoring the number of peaks
	sort(unique((plotmode$err_category)%%100))
}

plot.vb.my.figures=function(nameIn, t, tw, wmin, wmax, ub, uk, vb, vk, mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h, plot_legend, plot_example=FALSE){
	x.seq = seq(1.0,2.0,length=5)
	y.seq = seq(0.0,2.0,length=5)
	name = paste(nameIn,"_vb-my","_B",beta,"_uk",10*uk,"_vK",10*vk,"_lm",10*lm,"_mX",10*mx,"_mY","[y]","_vb","[x]","_omega",omega*10,"_r",r,"_c",cost*100,"_phi",phi*100,sep = "")
	
	#plot multiple results of simulations with changing v0 and r
	if(plot_example){
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
				L = calc.light_effect(t, mu,rho,kappa,sigma)
				#run simulation
				Ans = tss_probforage_energygain_optimize_linear(V, U, K, C, L, my, phi, omega, beta, h, mb,mx)
				#plot simulation results
				plot.sim_result(Ans,bquote(list('v'['b']==.(vb),'m'['y']==.(my))),L)
			}
		}
		dev.off()
	}
	
	
	#plot categories of simulation results with changing v0 and r
	grid = 101
	x.ax = seq(1.0,2.0,length=grid)
	y.ax = seq(0.0,2.5,length=grid)
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
			L = calc.light_effect(t, mu,rho,kappa,sigma)
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
	image.plotmode(x.ax,y.ax,plotmode,xlab="",ylab="",plot_legend = plot_legend)
	dev.off()
	
	#list of categorization error (grey colors) 
	plotmode$err_category
	#list of categorization error (grey colors) by ignoring the number of peaks
	sort(unique((plotmode$err_category)%%100))
}

plot.vb.r.figures=function(nameIn, t, tw, wmin, wmax, ub, uk, vb, vk, mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h, plot_legend, plot_example=FALSE){
	y.seq = seq(0.2,5.00,length=4)
	x.seq = seq(1.0,2.0,length=5)
	
	name = paste(nameIn,"_vb-r","_B",beta,"_uk",10*uk,"_vk",10*vk,"_lm",10*lm,"_mX",10*mx,"_mY",10*my,"_vb","[x]","_omega",omega*10,"_r","[y]","_c",cost*100,"_phi",phi,sep = "")
	if(plot_example){
		png(paste("examples_",name,".png",sep=""),height=2000,width=2000)
		par(mfrow=c(length(y.seq),length(x.seq)),cex=2.0,mex=0.3)
		for(r in rev(y.seq)){
			for(vb in x.seq){
				# calculate time-depending parameters 
				watertemp=calc.watertemp(t,tw,wmin,wmax)
				sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
				U = ub + uk*(watertemp-(wmax+wmin)/2)
				V = vb + vk*(sharktemp-(wmax+wmin)/2)
				K = rep(alpha,length=length(t))
				C = rep(cost,length=length(t))
				L = calc.light_effect(t, mu,rho,kappa,sigma)
				#run simulation
				Ans = tss_probforage_energygain_optimize_linear(V, U, K, C, L, my, phi, omega, beta, h, mb,mx)
				#plot simulation results
				plot.sim_result(Ans,bquote(list("r"==.(r),'phi'==.(phi))),L)
			}
		}
		dev.off()
	}
	
	
	#plot categories of simulation results with changing v0 and phi
	grid = 101
	y.ax = seq(0.05,6.00,length=grid)
	x.ax = seq(1.0,2.5,length=grid)
	no = matrix(0,grid,grid)
	for(y in 1:length(y.ax)){
		r = y.ax[y]
		for(x in 1:length(x.ax)){
			vb = x.ax[x]
			# calculate time-depending parameters 
			watertemp=calc.watertemp(t,tw,wmin,wmax)
			sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
			U = ub + uk*(watertemp-(wmax+wmin)/2)
			V = vb + vk*(sharktemp-(wmax+wmin)/2)
			K = rep(alpha,length=length(t))
			C = rep(cost,length=length(t))
			L = calc.light_effect(t, mu,rho,kappa,sigma)
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
	image.plotmode(x.ax,y.ax,plotmode,xlab="",ylab="",plot_legend = plot_legend)
	dev.off()
	
	#list of categorization error (grey colors) 
	plotmode$err_category
	#list of categorization error (grey colors) by ignoring the number of peaks
	sort(unique((plotmode$err_category)%%100))
}

plot.my.phi.figures=function(nameIn, t, tw, wmin, wmax, ub, uk, vb, vk, mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h, plot_legend, plot_example=FALSE){
	x.seq = seq(0.0,2.0,length=5)
	y.seq = seq(0.02,0.24,length=5)
	name = paste(nameIn,"_my-phi","_B",beta,"_uk",10*uk,"_vK",10*vk,"_lm",10*lm,"_mX",10*mx,"_mY","[x]","_vb",vb,"_omega",omega*10,"_r",r,"_c",cost*100,"_phi","[y]",sep = "")
	
	#plot multiple results of simulations with changing v0 and r
	if(plot_example){
		png(paste("examples_",name,".png",sep=""),height=2000,width=2000)
		par(mfrow=c(length(y.seq),length(x.seq)),cex=2.0,mex=0.3)
		for(phi in rev(y.seq)){
			for(my in x.seq){
				# calculate time-depending parameters
				watertemp=calc.watertemp(t,tw,wmin,wmax)
				sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
				U = ub + uk*(watertemp-(wmax+wmin)/2)
				V = vb + vk*(sharktemp-(wmax+wmin)/2)
				K = rep(alpha,length=length(t))
				C = rep(cost,length=length(t))
				L = calc.light_effect(t, mu,rho,kappa,sigma)
				#run simulation
				Ans = tss_probforage_energygain_optimize_linear(V, U, K, C, L, my, phi, omega, beta, h, mb,mx)
				#plot simulation results
				plot.sim_result(Ans,bquote(list('v'['b']==.(vb),'m'['y']==.(my))),L)
			}
		}
		dev.off()
	}
	
	
	#plot categories of simulation results with changing v0 and r
	grid = 101
	x.ax = seq(0.0,2.0,length=grid)
	y.ax = seq(0.0,0.25,length=grid)
	no = matrix(0,grid,grid)
	for(y in 1:length(y.ax)){
		phi = y.ax[y]
		for(x in 1:length(x.ax)){
			my = x.ax[x]
			
			# calculate time-depending parameters       
			watertemp=calc.watertemp(t,tw,wmin,wmax)
			sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
			U = ub + uk*(watertemp-(wmax+wmin)/2)
			V = vb + vk*(sharktemp-(wmax+wmin)/2)
			K = rep(alpha,length=length(t))
			C = rep(cost,length=length(t))
			L = calc.light_effect(t, mu,rho,kappa,sigma)
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
	image.plotmode(x.ax,y.ax,plotmode,xlab="",ylab="",plot_legend = plot_legend)
	dev.off()
	
	#list of categorization error (grey colors) 
	plotmode$err_category
	#list of categorization error (grey colors) by ignoring the number of peaks
	sort(unique((plotmode$err_category)%%100))
}

plot.uk.my.figures=function(nameIn, t, tw, wmin, wmax, ub, uk, vb, vk, mu,rho,kappa,sigma, alpha, omega, phi, mb, mx, my, r, cost, beta, h, plot_legend, plot_example=FALSE){
	x.seq = seq(0.0,4.0,length=5)
	y.seq = seq(0.0,2.0,length=5)
	name = paste(nameIn,"_uk-my","_B",beta,"_uk","[x]","_vK",10*vk,"_lm",10*lm,"_mX",10*mx,"_mY","[y]","_vb",vb,"_omega",omega*10,"_r",r,"_c",cost*100,"_phi",phi*100,sep = "")
	
	#plot multiple results of simulations with changing v0 and r
	if(plot_example){
		png(paste("examples_",name,".png",sep=""),height=2000,width=2000)
		par(mfrow=c(length(y.seq),length(x.seq)),cex=2.0,mex=0.3)
		for(my in rev(y.seq)){
			for(uk in x.seq){
				# calculate time-depending parameters
				watertemp=calc.watertemp(t,tw,wmin,wmax)
				sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
				U = ub + uk*(watertemp-(wmax+wmin)/2)
				V = vb + vk*(sharktemp-(wmax+wmin)/2)
				K = rep(alpha,length=length(t))
				C = rep(cost,length=length(t))
				L = calc.light_effect(t, mu,rho,kappa,sigma)
				#run simulation
				Ans = tss_probforage_energygain_optimize_linear(V, U, K, C, L, my, phi, omega, beta, h, mb,mx)
				#plot simulation results
				plot.sim_result(Ans,bquote(list('v'['b']==.(vb),'m'['y']==.(my))),L)
			}
		}
		dev.off()
	}
	
	
	#plot categories of simulation results with changing v0 and r
	grid = 101
	x.ax = seq(0.0,5.0,length=grid)
	y.ax = seq(0.0,2.0,length=grid)
	no = matrix(0,grid,grid)
	for(y in 1:length(y.ax)){
		my = y.ax[y]
		for(x in 1:length(x.ax)){
			uk = x.ax[x]
			
			# calculate time-depending parameters       
			watertemp=calc.watertemp(t,tw,wmin,wmax)
			sharktemp=calc.sharktemp(t,tw,wmin,wmax,r) 
			U = ub + uk*(watertemp-(wmax+wmin)/2)
			V = vb + vk*(sharktemp-(wmax+wmin)/2)
			K = rep(alpha,length=length(t))
			C = rep(cost,length=length(t))
			L = calc.light_effect(t, mu,rho,kappa,sigma)
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
	image.plotmode(x.ax,y.ax,plotmode,xlab="",ylab="",plot_legend = plot_legend)
	dev.off()
	
	#list of categorization error (grey colors) 
	plotmode$err_category
	#list of categorization error (grey colors) by ignoring the number of peaks
	sort(unique((plotmode$err_category)%%100))
}
