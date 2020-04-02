server = function(input,output,session){
  hideElement("slider")
  hideTab("tabs","Output")
  hideElement("fracw")
  hideElement("P")
  

  observeEvent(input$fluid,{
    if (input$fluid == "vwe" || input$fluid == "vwd") {
      showElement("fracw")
    } else {
      hideElement("fracw")
    }
  })
  
  
  
  observeEvent(input$start,{
    
    # Step change (m)
    dz = input$L/input$modelsteps
    
    # Inner pipe surface area (m2)
    SA_it = pi*input$tid*input$L
    
    # Outer pipe surface area (m2)
    SA_ot = pi*(2*input$tt + input$tid)*input$L
    
    # Inner shell surface area (m2)
    SA_is = pi*input$sid*input$L
    
    # Outer shell surface area (m2)
    SA_os = pi*(2*input$st + input$sid)*input$L
    
    # Number of tubes per row vector
    Ntubes_row = as.numeric(unlist(strsplit(input$ntubesrow,split = ",")))
    
    # Liquid water properties
    dens_lw = function(Temp) {
      Tau = 1 - Temp/647.096
      y = 17.863 + 58.606*Tau^0.35 - 95.396*Tau^(2/3) + 213.89*Tau - 141.26*Tau^(4/3) # mol/dm^3
      y = y*18.015 # kg/m3
      return(y)
    }
    
    cp_lw = function(Temp) {
      y = 276370 - 2090.1*Temp + 8.125*Temp^2 - 0.014116*Temp^3 + 9.3701*10^(-6)*Temp^4 # J/(kmol K)
      y = y/18.015 # J/(kg K)
      return(y)
    }
    
    mu_lw = function(Temp) {
      y = exp(-52.843 + 3703.6/Temp + 5.866*log(Temp) -5.879*10^(-29)*T^10) # Pa s
      return(y)
    }
    
    k_lw = function(Temp) {
      y = -0.432 + 0.0057255*Temp - 0.000008078*Temp^2 + 1.861*10^(-9)*Temp^3 # W/(m K)
      return(y)
    }
    
    
    # Liquid ethanol properties
    dens_le = function(Temp) {
      
      y = 1.6288/(0.27469^(1+(1-Temp/514)^0.23178)) # mol/dm^3
      y = y*46.068 # kg/m3
      return(y)
    }
    
    cp_le = function(Temp) {
      y = 102640 - 139.63*Temp - 0.030341*Temp^2 + 0.0020386*Temp^3 # J/(kmol K)
      y = y/46.068 # J/(kg K)
      return(y)
    }
    
    mu_le = function(Temp) {
      y = exp(7.875 + 781.98/Temp -3.0418*log(Temp)) # Pa s
      return(y)
    }
    
    k_le = function(Temp) {
      y = 0.2468 - 0.000264*Temp  # W/(m K)
      return(y)
    }
    
    
    #Liquid decane properties
    dens_ld = function(Temp) {
      
      y = 0.41084/(0.25175^(1+(1-Temp/617.7)^0.28571)) # mol/dm^3
      y = y*142.282 # kg/m3
      return(y)
    }
    
    cp_ld = function(Temp) {
      y = 278620 - 197.91*Temp + 1.0737*Temp^2 # J/(kmol K)
      y = y/142.282 # J/(kg K)
      return(y)
    }
    
    mu_ld = function(Temp) {
      y = exp(-9.6489 + 1181.1/Temp -0.24367*log(Temp) + 9.0522*10^34*Temp^-15) # Pa s
      return(y)
    }
    
    k_ld = function(Temp) {
      y = 0.2063 - 0.00025*Temp  # W/(m K)
      return(y)
    }
    
    
    
    # Vapour water properties
    hvap_vw = function(Temp) {
      Tau = 1 - Temp/647.096
      y = 5.2053*10^7*(1-Tau)^(0.3199 - 0.212*Tau + 0.25795*Tau^2) # J/kmol
      y = y/18.015 # J/kg
      return(y)
    }
    
    cp_vw = function(Temp) {
      y = 0.33363*10^5 + 0.2697*10^5*((2.6105*10^3/Temp)/(sinh((2.6105*10^3)/Temp)))^2 + 0.08896*10^5*((1169/Temp)/(cosh((1169)/Temp)))^2 # J/(kmol K)
      y = y/18.015 # J/(kg K)
      return(y)
    }
    
    mu_vw = function(Temp) {
      y = 1.7096*10^-8*Temp^1.1146 # Pa s
      return(y)
    }
    
    k_vw = function(Temp) {
      y = 6.2041*10^-6*Temp^1.3973 # W/m K
      return(y)
    }
    
    
    # Vapour ethanol properties
    hvap_ve = function(Temp) {
      Tau = Temp/514
      y = 5.5789*10^7*(1-Tau)^(0.31245) # J/kmol
      y = y/46.068 # J/kg
      return(y)
    }
    
    cp_ve = function(Temp) {
      y = 0.492*10^5 + 1.4577*10^5*((1.6628*10^3/Temp)/(sinh((1.6628*10^3)/Temp)))^2 + 0.939*10^5*((744.7/Temp)/(cosh((744.7)/Temp)))^2 # J/(kmol K)
      y = y/46.068 # J/(kg K)
      return(y)
    }
    
    mu_ve = function(Temp) {
      y = (1.0613*10^-7*Temp^0.8066)/(1+52.7/Temp) # Pa s
      return(y)
    }
    
    k_ve = function(Temp) {
      y = (-0.010109*Temp^0.6475)/(1-7332/Temp-268000/(Temp^2)) # W /(m K)
      return(y)
    }
    
    
    
    # Vapour decane properties
    hvap_vd = function(Temp) {
      Tau = Temp/617.7
      y = 6.6126*10^7*(1-Tau)^(0.39797) # J/kmol
      y = y/142.28 # J/kg
      return(y)
    }
    
    cp_vd = function(Temp) {
      y = 1.672*10^5 + 5.353*10^5*((1.6141*10^3/Temp)/(sinh((1.6141*10^3)/Temp)))^2 + 3.782*10^5*((742/Temp)/(cosh((742)/Temp)))^2 # J/(kmol K)
      y = y/142.28 # J/(kg K)
      return(y)
    }
    
    mu_vd = function(Temp) {
      y = (2.64*10^-8*Temp^0.9487)/(1+71/Temp) # Pa s
      return(y)
    }
    
    k_vd = function(Temp) {
      y = (-668.4*Temp^0.9323)/(1-4.071*10^9/Temp) # W /(m K)
      return(y)
    }
    
    
    
    # Vapour water + ethanol properties
    hvap_vwe = function(x,y,Temp) {
      h1 = hvap_vw(Temp)
      h2 = hvap_ve(Temp)
      
      a0 = -3.63868*10^5 + 1.83829*Temp - 2.32763*Temp^2
      a1 = 9.25982*10^5 - 4.83586*Temp + 6.37228*Temp^2
      a2 = -14.04894*10^5 + 7.51661*Temp - 10.11280*Temp^2
      a3 = 10.91318*10^5 - 5.89498*Temp + 7.98868*Temp^2
      a5 = -2.79986*10^5 + 1.50557*Temp - 2.03127*Temp^2
      
      # convert to mol fraction
      ymol = y/(18.015)/(y/18.015+(1-y)/46.068)
      xmol = x/(18.015)/(x/18.015+(1-x)/46.068)
      
      hE = xmol*(1-xmol)*(a0+a1*(1-xmol)^0.5+a2*(1-xmol)^1.5+a3*(1-xmol)^2.5+a5*(1-xmol)^4.5) 
      
      hL = (1-x)*cp_le(Temp)*(Temp-273.15)+x*cp_lw(Temp)*(Temp-273.15) + hE
      
      hV = (1-y)*(h2+cp_ve(Temp)*(Temp-273.15))+y*(h1+cp_vw(Temp)*(Temp-273.15))
      
      h = hV - hL # J/kg
      return(h)
    }
    
    cp_vwe = function(y,Temp) {
      c1 = cp_vw(Temp)
      c2 = cp_ve(Temp)
      
      c = y*c1 + (1-y)*c2
      return(c) 
    }
    
    mu_vwe = function(y,Temp) {
      m1 = mu_vw(Temp)
      m2 = mu_ve(Temp)
      
      # convert to mol fraction
      y = y/(18.015)/(y/18.015+(1-y)/46.068)
      
      m = (y*m1*sqrt(18.015)+(1-y)*m2*sqrt(46.068))/(y*sqrt(18.015)+(1-y)*sqrt(46.068))
      return(m)
    }
    
    k_vwe = function(y,Temp) {
      k1 = k_vw(Temp)
      k2 = k_ve(Temp)
      
      phi12 = 1/sqrt(8)*(1+18.015/46.068)^-0.5*(1+sqrt(mu_vw(Temp)/mu_ve(Temp))*(46.068/18.015)^0.25)^2
      phi21 = 1/sqrt(8)*(1+46.068/18.015)^-0.5*(1+sqrt(mu_ve(Temp)/mu_vw(Temp))*(18.015/46.068)^0.25)^2
      
      # convert to mol fraction
      y = y/(18.015)/(y/18.015+(1-y)/46.068)
      
      k = y*k1/((1-y)*phi12 + y)+ (1-y)*k2/(y*phi21 +(1-y))
      return(k)
    }
    
    # Liquid water + ethanol properties
    dens_lwe = function(x,Temp) {
      A = -9.054199 + 0.01572057*Temp 
      B = -4.930763 + 0.01241796*Temp 
      C = 5.286817 - 0.01989784*Temp 
      
      # convert to mol fraction
      x = x/(18.015)/(x/18.015+(1-x)/46.068)
      
      vE = x*(1-x)*(A + B*(1-2*(1-x)) + C*(1-2*(1-x))^2)/1000000 # m3/mol
      
      MW = (x*18.015+(1-x)*46.068)/1000 # kg/mol
      
      d = MW/(vE+(x*0.018015)/dens_lw(Temp)+((1-x)*0.046068)/dens_le(Temp)) # kg/m3
      return(d)
    }
    
    cp_lwe = function(x,Temp) {
      c1 = cp_lw(Temp)
      c2 = cp_le(Temp)
      
      c = x*c1 + (1-x)*c2
      return(c) 
    }
    
    
    mu_lwe = function(x,Temp) {
      m1 = mu_lw(Temp)
      m2 = mu_le(Temp)
      
      # convert to mol fraction
      x = x/(18.015)/(x/18.015+(1-x)/46.068)
      
      xe = 1-x
      
      m = exp(x*log(m1)+xe*log(m2)+724.652*x*xe/Temp+729.357*x*xe*(x-xe)/Temp + 976.050*x*xe*(x-xe)^2/Temp)
      return(m)
    }
    
    k_lwe = function(x,Temp) {
      k1 = k_lw(Temp)
      k2 = k_le(Temp)
      
      k12 = 2/(1/k1 + 1/k2)
      
      
      # convert to mol fraction
      x = x/(18.015)/(x/18.015+(1-x)/46.068)
      
      phi1 = x/dens_lw(Temp)*18.015
      phi2 = (1-x)/dens_le(Temp)*46.068 
      
      phi11 = phi1/(phi1+phi2)
      phi22 = phi2/(phi1+phi2)
      
      k = phi11^2*k1 + 2*phi11*phi22*k12 + phi22^2*k2
      return(k)
    }
    
    
    Tcondy = function(y) {
      
      z = vector(mode = "numeric", length = length(y))
      if (input$fluid == "vwe") {
        for (i in 1:length(y)) {
          if (y[i] < 0.2) {
            z[i] = 88.383*y[i]^3 + 106.69*y[i]^2 - 12.263*y[i] + 78.315 + 273.15
          } else {
            z[i] = 19.169*y[i]^3 - 58.99*y[i]^2 + 71.195*y[i] + 68.631 + 273.15 
          }
        }
      } else if (input$fluid == "vwd") {
        for (i in 1:length(y)) {
          if (y[i] < 0.587583) {
            z[i] = 2654.4*y[i]^6 - 6221.4*y[i]^5 + 5828.1*y[i]^4 - 2967.4*y[i]^3 + 968.55*y[i]^2 - 301.14*y[i] + 174.03 + 273.15
          } else {
            z[i] = 5.8327*y[i] + 94.515 + 273.15 
          }
        }
      }
      return(z)
    }
    
    
    Tcondx = function(x) {
      if (input$fluid == "vwe") {
        z = vector(mode = "numeric", length = length(x))
        for (i in 1:length(x)) {
          z[i] = 330.33*x[i]^6 - 892.47*x[i]^5 + 978.2*x[i]^4 - 535.54*x[i]^3 + 154.3*x[i]^2 - 13.25*x[i] + 78.319 + 273.15 - Tcond[i]
        }
      } else if (input$fluid == "vwd") {
        z = vector(mode = "numeric", length = length(x))
        for (i in 1:length(x)) {
          if (x[i] <= 0.01089) {
            z[i] = 542157*x[i]^2 - 13099*x[i] + 174.86 + 273.15
            if (z[i] < 97.9425 + 273.15) {
              z[i] = 97.9425 + 273.15 - Tcond[i]
            } else {
              z[i] = z[i] - Tcond[i]
            }
          } else if (x[i] <= 0.992156) {
            z[i] = 97.9425 + 273.15 - Tcond[i]
          } else {
            z[i] = 306.67*x[i] - 206.32 - Tcond[i]
          }
        }
      }
      return(z)
    }
    
    
    
    # Vapour water + decane properties
    hvap_vwd = function(x,Temp) {
      h1 = hvap_vw(Temp)
      h2 = hvap_vd(Temp)
      h = (x/dens_lw(Temp)*h1 + (1-x)/dens_ld(Temp)*h2)/(x/dens_lw(Temp)+(1-x)/dens_ld(Temp)) # J/kg
      return(h)
    }
    
    cp_vwd = function(y,Temp) {
      c1 = cp_vw(Temp)
      c2 = cp_vd(Temp)
      
      c = y*c1 + (1-y)*c2
      return(c) 
    }
    
    mu_vwd = function(y,Temp) {
      m1 = mu_vw(Temp)
      m2 = mu_vd(Temp)
      
      # convert to mol fraction
      y = y/(18.015)/(y/18.015+(1-y)/142.282)
      
      m = (y*m1*sqrt(18.015)+(1-y)*m2*sqrt(142.282))/(y*sqrt(18.015)+(1-y)*sqrt(148.282))
      return(m)
    }
    
    k_vwd = function(y,Temp) {
      k1 = k_vw(Temp)
      k2 = k_vd(Temp)
      
      phi12 = 1/sqrt(8)*(1+18.015/142.282)^-0.5*(1+sqrt(mu_vw(Temp)/mu_vd(Temp))*(142.282/18.015)^0.25)^2
      phi21 = 1/sqrt(8)*(1+142.282/18.015)^-0.5*(1+sqrt(mu_vd(Temp)/mu_vw(Temp))*(18.015/142.282)^0.25)^2
      
      # convert to mol fraction
      y = y/(18.015)/(y/18.015+(1-y)/142.282)
      
      k = y*k1/((1-y)*phi12 + y)+ (1-y)*k2/(y*phi21 +(1-y))
      return(k)
    }
    
    
    
    # Liquid water + decane properties
    dens_lwd = function(x,Temp) {
      d1 = dens_lw(Temp)
      d2 = dens_ld(Temp)
      d = x*d1 + (1-x)*d2
      return(d)
    }
    
    cp_lwd = function(x,Temp) {
      c1 = cp_lw(Temp)
      c2 = cp_le(Temp)
      c = x*c1 + (1-x)*c2
      return(c) 
    }
    
    
    mu_lwd = function(x,Temp) {
      m1 = mu_lw(Temp)
      m2 = mu_le(Temp)
      m = x*m1 + (1-x)*m2
      return(m)
    }
    
    k_lwd = function(x,Temp) {
      k1 = k_lw(Temp)
      k2 = k_le(Temp)
      k = x*k1 + (1-x)*k2
      return(k)
    }
    
    
    
    
    # Double-pipe or shell and tube, no baffles --------------------------------------------------------------------------------------------  
    if (input$baffles == 0) {withProgress(value = 0.5, message = "Calculation in progress...",{
      
      # Cross sectional area of the shell (m2)
      CSA_s = pi/4*(input$sid^2-sum(Ntubes_row)*(input$tid+2*input$tt)^2)
      
      # Cross sectional area of the tube (m2)
      CSA_t = pi/4*(input$tid^2)
      
      # Hydraulic diameter shell/pipe
      dhs = ((input$sid)^2-sum(Ntubes_row)*(input$tid+2*input$tt)^2)/(input$sid+sum(Ntubes_row)*(input$tid+2*input$tt))
      dht = input$tid
      
      # Inital temperature estimate vector
      Temp = c(rep(input$Tsi,times = input$modelsteps),
               rep(input$Tti,times = input$modelsteps),
               rep(mean(c(input$Tti,input$Tsi)),times = input$modelsteps))
      SensTempOld = Temp
      CondTempOld = Temp
      
      # Create mass flowrate vector
      Fs = rep(input$Fs,input$modelsteps)
      
      # Create condensed mass vector
      Mcond = rep(0,input$modelsteps)
      
      # Create RH vector
      RH = rep(input$RH,input$modelsteps)
      
      
      err = 100
      count = 1
      while (err > input$error) { 
        # Create matrix containing all calclated temp values
        T_matrix = matrix(data = 0, ncol = input$modelsteps*3,nrow = input$modelsteps*3)
        
        # Properties
        Dens = dens_lw(Temp)
        Cp = cp_lw(Temp)
        Mu = mu_lw(Temp)
        K = k_lw(Temp)
        
        # Reynolds #
        vs = Fs/(Dens[1:input$modelsteps]*CSA_s)
        vt = input$Ft/(Dens[(input$modelsteps+1):(2*input$modelsteps)]*CSA_t*sum(Ntubes_row))
        Re_s = Dens[1:input$modelsteps]*dhs*vs/Mu[1:input$modelsteps]
        Re_t = Dens[(input$modelsteps+1):(2*input$modelsteps)]*dht*vt/Mu[(input$modelsteps+1):(2*input$modelsteps)]
        
        # Prandtl #
        Pr_s = Mu[1:input$modelsteps]*Cp[1:input$modelsteps]/K[1:input$modelsteps]
        Pr_t = Mu[(input$modelsteps+1):(2*input$modelsteps)]*Cp[(input$modelsteps+1):(2*input$modelsteps)]/K[(input$modelsteps+1):(2*input$modelsteps)]
        
        # Graetz #
        Gz_s = dhs/input$L*Re_s*Pr_s
        Gz_t = dht/input$L*Re_t*Pr_t
        
        # Nusselt #
        if (mean(Re_s) > 2100) {
          Nu_s = 0.027*Re_s^0.8*Pr_s^(1/3)*(Mu[1:input$modelsteps]/Mu[(2*input$modelsteps+1):(3*input$modelsteps)])^0.14
        } else {
          Nu_s = 1.86*Gz_s^(1/3)*(Mu[1:input$modelsteps]/Mu[(2*input$modelsteps+1):(3*input$modelsteps)])^0.14
        }
        
        if (mean(Re_t) > 2100) {
          Nu_t = 0.027*Re_t^0.8*Pr_t^(1/3)*(Mu[(input$modelsteps+1):(2*input$modelsteps)]/Mu[(2*input$modelsteps+1):(3*input$modelsteps)])^0.14
        } else {
          Nu_t = 1.86*Gz_t^(1/3)*(Mu[(input$modelsteps+1):(2*input$modelsteps)]/Mu[(2*input$modelsteps+1):(3*input$modelsteps)])^0.14
        }
        
        # Heat transfer coeffcient
        h_s = Nu_s*K[1:input$modelsteps]/dhs
        h_t = Nu_t*K[(input$modelsteps+1):(2*input$modelsteps)]/dht
        
        testhsdp <<- h_s
        
        # Insert coefficents for shell
        T_matrix[1,1] = 1
        for (i in 2:input$modelsteps){
          T_matrix[i,i] = 1
          T_matrix[i,i-1] = -1 + sum(Ntubes_row)*h_s[i-1]*SA_ot*dz/(input$L*input$Fs*Cp[i-1])
          T_matrix[i,i-1+2*input$modelsteps] = -sum(Ntubes_row)*h_s[i-1]*SA_ot*dz/(input$L*input$Fs*Cp[i-1])
        }
        
        # Insert coefficents for tube
        T_matrix[2*input$modelsteps,2*input$modelsteps] = 1
        for (i in (input$modelsteps+1):(2*input$modelsteps-1)) {
          T_matrix[i,i] = 1
          T_matrix[i,i+1] = -1 + sum(Ntubes_row)*h_t[i+1-input$modelsteps]*SA_it*dz/(input$L*input$Ft*Cp[i+1])
          T_matrix[i,i+1+input$modelsteps] = -sum(Ntubes_row)*h_t[i+1-input$modelsteps]*SA_it*dz/(input$L*input$Ft*Cp[i+1])
        }
        
        # Insert coefficents for wall
        for (i in (2*input$modelsteps+1):(3*input$modelsteps)) {
          T_matrix[i,i] = -1
          T_matrix[i,i-input$modelsteps] = h_t[i-2*input$modelsteps]*SA_it/(h_s[i-2*input$modelsteps]*SA_ot + h_t[i-2*input$modelsteps]*SA_it)
          T_matrix[i,i-2*input$modelsteps] = h_s[i-2*input$modelsteps]*SA_ot/(h_s[i-2*input$modelsteps]*SA_ot + h_t[i-2*input$modelsteps]*SA_it)
        }
        
        # Right hand side solution to matrix equation
        rhs_vector = rep(0, input$modelsteps*3)
        rhs_vector[1] = input$Tsi
        rhs_vector[2*input$modelsteps] = input$Tti
        
        # Solve matrix problem using LU decomposition
        x_vector = lusys(T_matrix,rhs_vector)
        
        err = sqrt(mean(((x_vector-Temp)/Temp)^2))*100
        Temp = x_vector
        print(err)
      }
      
      # Initialize heatmap data matrix without baffle columns
      heatmapdata = matrix(ncol = input$modelsteps, nrow = 2*length(Ntubes_row)+1)
      
      for (i in 1:(2*length(Ntubes_row)+1)) {
        if (i %% 2 == 0) {
          heatmapdata[i,1:input$modelsteps] = x_vector[(input$modelsteps+1):(2*input$modelsteps)]
        } else {
          heatmapdata[i,1:input$modelsteps] = x_vector[1:input$modelsteps]        
        }
      }
      
      
      tubelines = list()
      # Create border lines for tubes in heatmap
      for (i in 1:(2*length(Ntubes_row))) {
        tubelines = c(tubelines,
                      list(list(type = 'line', xref = "x", yref = "y",
                                x0 = -0.5, x1 = input$modelsteps-0.5,
                                y0 = i-0.5, y1 = i-0.5,
                                line = list(color = "black", width = 10)
                      )
                      )
        )
      }
      
      heatmaptext = list()
      for (i in 1:nrow(heatmapdata)) {
        heatmaptext[[i]] = paste(round(heatmapdata[i,1:input$modelsteps],digits = 2), "K", sep = " ")
      }
      
      # Create plot and output
      output$plot1 = renderPlotly({
        plot_ly(z = heatmapdata, type = "heatmap", hoverinfo = 'text', 
                text = heatmaptext,
                colorscale = "RdBu", 
                colorbar = list(len = 1, title = list(text = "Temperature (K)",side = "right"))) %>%
          config(displayModeBar = F) %>%
          layout(
            shapes = tubelines,
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            xaxis = list(
              ticks = "",
              showticklabels = F,
              showgrid = F,
              zeroline = F,
              fixedrange = T
            ),
            yaxis = list(
              ticks = "",
              showticklabels = F,
              showgrid = F,
              zeroline = F,
              fixedrange = T
            )
          )
      })
    })
      
      # Double-pipe and shell and tube with baffles version 2 --------------------------------------------------------------------------------------------  
    } else if (input$baffles > 0) {withProgress(value = 0.5, message = "Calculation in progress...",{
      
      # Number of sections seperated by baffles
      sections = input$baffles + 1
      
      # Number of model steps per section
      stepsPerSection = ceiling(input$modelsteps/sections)
      
      # Cross sectional area of the tube (m2)
      CSA_t = pi/4*(input$tid^2)
      
      # Hydraulic diameter tube
      dht = input$tid
      
      # Calculate space between tubes (half a space on either side)
      pitch_y = input$sid/(length(Ntubes_row))
      
      # Chord length and cross-sectional area of shell per row
      chord_length = 2*sqrt(2*pitch_y*((1:length(Ntubes_row))-0.5)*input$sid/2 - (pitch_y*((1:length(Ntubes_row))-0.5))^2)
      CSA_s = c(chord_length*input$L/(sections*stepsPerSection),0)
      
      # Calculate space between tubes (half a space on either side)
      maxrowindex = match(max(Ntubes_row),Ntubes_row)
      pitch_x = min(chord_length[Ntubes_row %in% max(Ntubes_row)])/max(Ntubes_row)
      
      CSA_s_vector = vector()
      for (i in 1:length(CSA_s)) {
        for (j in 1:sections){
          
          if (j %% 2 == 0){
            CSA_s_vector[(1+stepsPerSection*(j-1) + (i-1)*stepsPerSection*sections):(j*stepsPerSection + (i-1)*stepsPerSection*sections)] = rev(CSA_s)[i]
          } else {
            CSA_s_vector[(1+stepsPerSection*(j-1) + (i-1)*stepsPerSection*sections):(j*stepsPerSection + (i-1)*stepsPerSection*sections)] = CSA_s[i]
          }
        }
      }
      
      # Number of tube rows
      nrows = length(Ntubes_row)
      
      # Inital temperature estimate vector
      Temp = c(rep(input$Tsi,times = (nrows+1)*sections*stepsPerSection),
               rep(input$Tti,times = nrows*sections*stepsPerSection),
               rep(mean(c(input$Tti,input$Tsi)),times = nrows*sections*stepsPerSection))
      SensTempOld = Temp
      CondTempOld = Temp
      
      # Create mass flowrate vector
      Fs = rep(input$Fs/stepsPerSection,stepsPerSection*sections*(nrows+1))
      
      # Create condensed mass vector
      Mcond = rep(0,stepsPerSection*sections*(nrows+1))
      
      # Create RH vector
      RH = rep(input$RH,stepsPerSection*sections*(nrows+1))
      
      # Initialize x and y vectors for mixtures
      if (input$fluid == "vwe" || input$fluid == "vwd") {
        yw = input$fracw
        ye = 1 - yw
        y = rep(yw,sections*stepsPerSection*(nrows+1))
        x = rep(0.5,sections*stepsPerSection*(nrows+1))
      }
      
      err = 100
      count = 1
      while (err > input$error && count <= input$iterations) { 
        
        # Create matrix containing all calculated temp values
        T_matrix = matrix(data = 0, ncol = stepsPerSection*sections*(3*nrows+1),nrow = stepsPerSection*sections*(3*nrows+1))
        
        # Initial values + averaging for next section
        for (i in seq(1,sections*stepsPerSection,by = stepsPerSection*2)) {
          
          T_matrix[i:(stepsPerSection+i-1),i:(stepsPerSection+i-1)] = diag(stepsPerSection)
          
        }
        for (i in seq(sections*stepsPerSection*nrows+stepsPerSection+1,sections*stepsPerSection*(nrows+1),by = stepsPerSection*2)) {
          
          T_matrix[i:(stepsPerSection+i-1),i:(stepsPerSection+i-1)] = diag(stepsPerSection)
          
        }
        
        if (sections > 2){
          for (i in seq(1+stepsPerSection*2,sections*stepsPerSection,by = stepsPerSection*2)) {
            
            T_matrix[i:(stepsPerSection+i-1),(i-stepsPerSection):(i-1)] = matrix(data = -1/stepsPerSection,ncol = stepsPerSection, nrow = stepsPerSection)
            
          }
        }
        
        for (i in seq(sections*stepsPerSection*nrows+stepsPerSection+1,sections*stepsPerSection*(nrows+1),by = stepsPerSection*2)) {
          
          T_matrix[i:(stepsPerSection+i-1),(i-stepsPerSection):(i-1)] = matrix(data = -1/stepsPerSection,ncol = stepsPerSection, nrow = stepsPerSection)
          
        }
        # Right hand side solution to matrix equation
        rhs_vector = rep(0, stepsPerSection*sections*(3*nrows+1))
        rhs_vector[1:stepsPerSection] = input$Tsi
        
        # Known inlet tube stuff for matrix
        for (i in seq(sections*stepsPerSection*(nrows+2),sections*stepsPerSection*(2*nrows+1),by = sections*stepsPerSection)) {
          T_matrix[i,i] = 1
          rhs_vector[i] = input$Tti
        }
        
        if (input$fluid == "vw"){
          Tcond = rep(1730.63/(8.07131-log10(input$P/101325*760)) - 233.426 + 273.15,sections*stepsPerSection*(nrows+1))
          for (i in 1:(stepsPerSection*sections*(nrows+1))) {
            if (Temp[i] < Tcond[i]) {
              Temp[i] = Tcond[i]
            }
          }
        } else if (input$fluid == "vwe") {
          
          #Calculate azeotrope
          condFunc = function(z) {
            temp = Tcondy(z) - (330.33*z^6 - 892.47*z^5 + 978.2*z^4 - 535.54*z^3 + 154.3*z^2 - 13.25*z + 78.319 + 273.15)
            return(temp)
          }
          yfinal = nleqslv(0.05, condFunc)$x
          xfinal = yfinal
          Tazeo = Tcondy(yfinal)
          
          Tcond = Tcondy(y)
          testtcond <<- Tcond
          
          MWL = 1/(x/18.015+(1-x)/46.068)
          MWV = 1/(y/18.015+(1-y)/46.068)
          
          for (i in 1:(stepsPerSection*sections*(nrows+1))) {
            if (Temp[i] < Tcond[i]) {
              Temp[i] = Tcond[i]
            }
          }
        } else if (input$fluid == "vwd") {
          
          #Calculate azeotrope
          condFunc = function(z) {
            temp = Tcondy(z) - (97.9425 + 273.15)
            return(temp)
          }
          yfinal = nleqslv(0.5, condFunc)$x
          xfinal = yfinal
          Tazeo = Tcondy(yfinal)
          
          Tcond = Tcondy(y)
          testtcond <<- Tcond
          
          MWL = 1/(x/18.015+(1-x)/142.29)
          MWV = 1/(y/18.015+(1-y)/142.29)
          
          for (i in 1:(stepsPerSection*sections*(nrows+1))) {
            if (Temp[i] < Tcond[i]) {
              Temp[i] = Tcond[i]
            }
          }
        }
        
        
        # Properties
        Dens = dens_lw(Temp)
        Cp = cp_lw(Temp)
        Mu = mu_lw(Temp)
        K = k_lw(Temp)
        
        if (input$fluid == "vw") {
          Dens[1:(stepsPerSection*sections*(nrows+1))] = input$P*18.015/(8.314*Temp[1:(stepsPerSection*sections*(nrows+1))])/1000
          Cp[1:(stepsPerSection*sections*(nrows+1))] = cp_vw(Temp[1:(stepsPerSection*sections*(nrows+1))])
          Mu[1:(stepsPerSection*sections*(nrows+1))] = mu_vw(Temp[1:(stepsPerSection*sections*(nrows+1))])
          K[1:(stepsPerSection*sections*(nrows+1))] = k_vw(Temp[1:(stepsPerSection*sections*(nrows+1))])
          Hvap = hvap_vw(Temp[1:(stepsPerSection*sections*(nrows+1))])
        } else if (input$fluid == "vwe") {
          Dens[1:(stepsPerSection*sections*(nrows+1))] = input$P*MWV/(8.314*Temp[1:(stepsPerSection*sections*(nrows+1))])/1000
          Cp[1:(stepsPerSection*sections*(nrows+1))] = cp_vwe(y,Temp[1:(stepsPerSection*sections*(nrows+1))])
          Mu[1:(stepsPerSection*sections*(nrows+1))] = mu_vwe(y,Temp[1:(stepsPerSection*sections*(nrows+1))])
          K[1:(stepsPerSection*sections*(nrows+1))] = k_vwe(y,Temp[1:(stepsPerSection*sections*(nrows+1))])
          Hvap = hvap_vwe(x,y,Temp[1:(stepsPerSection*sections*(nrows+1))])
        } else if (input$fluid == "vwd") {
          Dens[1:(stepsPerSection*sections*(nrows+1))] = input$P*MWV/(8.314*Temp[1:(stepsPerSection*sections*(nrows+1))])/1000
          Cp[1:(stepsPerSection*sections*(nrows+1))] = cp_vwd(y,Temp[1:(stepsPerSection*sections*(nrows+1))])
          Mu[1:(stepsPerSection*sections*(nrows+1))] = mu_vwd(y,Temp[1:(stepsPerSection*sections*(nrows+1))])
          K[1:(stepsPerSection*sections*(nrows+1))] = k_vwd(y,Temp[1:(stepsPerSection*sections*(nrows+1))]) 
          Hvap = hvap_vwd(x,Temp[1:(stepsPerSection*sections*(nrows+1))])
        }
        
        # Reynolds #
        vs = Fs/(Dens[1:(stepsPerSection*sections*(nrows+1))]*CSA_s_vector)
        vt = input$Ft/(Dens[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*CSA_t*sum(Ntubes_row))
        
        vs_max = vs*pitch_x/(pitch_x-(input$tid + 2*input$tt))
        
        if (input$config == "triangle" && pitch_x > input$tid + 2*input$tt + 2*(sqrt((pitch_x/2)^2 + pitch_y^2) - input$tid - 2*input$tt)) {
          pitch_diag = sqrt((pitch_x/2)^2 + pitch_y^2)
          vs_max = pitch_x/2*vs/(pitch_diag - input$tid - 2*input$tt)
        } 
        Re_s = Dens[1:(stepsPerSection*sections*(nrows+1))]*(dht+2*input$tt)*vs_max/Mu[1:(stepsPerSection*sections*(nrows+1))]
        Re_t = Dens[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*dht*vt/Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]
        Re_d = Dens[1:(stepsPerSection*sections*(nrows+1))]*(dht+2*input$tt)*vs/Mu[1:(stepsPerSection*sections*(nrows+1))]
        
        
        # Prandtl #
        Pr_s = Mu[1:(stepsPerSection*sections*(nrows+1))]*Cp[1:(stepsPerSection*sections*(nrows+1))]/K[1:(stepsPerSection*sections*(nrows+1))]
        Pr_t = Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/K[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]
        
        # Create Temp_w vector that matches the shell vector
        Start = (stepsPerSection*sections*(2*nrows+1)+1)
        End = (stepsPerSection*sections*(3*nrows+1))
        Temp_w_aligned = c(Temp[Start:End],rep(NA,stepsPerSection*sections))
        for (i in seq(2,sections, by = 2)) {
          for (j in (nrows+1):1) {
            if (j == 1){
              Temp_w_aligned[((j-1)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-1)*sections*stepsPerSection + (i)*stepsPerSection)] = NA
            } else {
              Temp_w_aligned[((j-1)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-1)*sections*stepsPerSection + (i)*stepsPerSection)] = 
                Temp[((j-2)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-2)*sections*stepsPerSection + (i)*stepsPerSection)]
            }
          }
        }
        
        
        if (input$fluid == "vw") {
          Pr_w_aligned = mu_vw(Temp_w_aligned)*cp_vw(Temp_w_aligned)/k_vw(Temp_w_aligned)
        } else if (input$fluid == "vwe") {
          Pr_w_aligned = mu_vwe(y,Temp_w_aligned)*cp_vwe(y,Temp_w_aligned)/k_vwe(y,Temp_w_aligned)
        } else if (input$fluid == "vwd") {
          Pr_w_aligned = mu_vwd(y,Temp_w_aligned)*cp_vwd(y,Temp_w_aligned)/k_vwd(y,Temp_w_aligned)
        } else {
          # Create Pr_w vector that matches the shell vector
          Pr_w_aligned = mu_lw(Temp_w_aligned)*cp_lw(Temp_w_aligned)/k_lw(Temp_w_aligned)
        }
        
        # Graetz #
        Gz_t = dht/input$L*Re_t*Pr_t
        
        # Nusselt #
        Nuss_func = function(C,m) {
          Nuss = C*(Re_s^m)*(Pr_s^0.36)*(Pr_s/Pr_w_aligned)^0.25
          return(Nuss)
        } 
        
        # Double pipe vs shell and tube if statement
        if (sum(Ntubes_row) == 1) {
          Nu_s = 0.3 + (0.62*(Re_d^0.5)*Pr_s^(1/3))/((1+(0.4/Pr_s)^(2/3))^0.25)*((1+(Re_d/282000)^0.625)^(0.8))
        } else {
          # Configuration if statement
          if (input$config == "square") {
            
            # Reynolds if statement
            if (mean(Re_s[Re_s != Inf],na.rm = T) < 100) {
              Nu_s = Nuss_func(0.8,0.4)
            } else if (mean(Re_s[Re_s != Inf],na.rm = T)  < 1000) {
              Nu_s = 0.3 + (0.62*(Re_d^0.5)*Pr_s^(1/3))/((1+(0.4/Pr_s)^(2/3))^0.25)*((1+(Re_d/282000)^0.625)^(0.8))
            } else if (mean(Re_s[Re_s != Inf],na.rm = T)  < 200000) {
              Nu_s = Nuss_func(0.27,0.63)
            } else {
              Nu_s = Nuss_func(0.021,0.84)
            }
            
          } else {
            if (mean(Re_s[Re_s != Inf],na.rm = T)  < 100) {
              Nu_s = Nuss_func(0.9,0.4)
            } else if (mean(Re_s[Re_s != Inf],na.rm = T)  < 1000) {
              Nu_s = 0.3 + (0.62*(Re_d^0.5)*Pr_s^(1/3))/((1+(0.4/Pr_s)^(2/3))^0.25)*((1+(Re_d/282000)^0.625)^(0.8))
            } else if (mean(Re_s[Re_s != Inf],na.rm = T)  < 200000) {
              # Pitch if statement
              if (pitch_x/pitch_y < 2) {
                Nu_s = Nuss_func(0.35*(pitch_x/pitch_y)^0.2,0.6)
              } else {
                Nu_s = Nuss_func(0.4,0.6)
              }
            } else {
              Nu_s = Nuss_func(0.022,0.84)
            }
          }
        }
        
        testres <<- Re_s
        testNus <<- Nu_s
        testPrs <<- Pr_s
        testvsm <<- vs_max
        testvs <<- vs
        
        
        if (mean(Re_t) > 2100) {
          Nu_t = 0.027*Re_t^0.8*Pr_t^(1/3)*(Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/Mu[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))])^0.14
        } else {
          Nu_t = 1.86*Gz_t^(1/3)*(Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/Mu[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))])^0.14
        }
        
        # Heat transfer coeffcient
        h_s = Nu_s*K[1:(stepsPerSection*sections*(nrows+1))]/(dht+2*input$tt)
        h_t = Nu_t*K[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/dht

        # Matrix coefficients
        Tube_vector = vector()
        for (i in 1:(nrows+1)) {
          for (j in 1:sections){
            if (j %% 2 == 0){
              Tube_vector[(1+stepsPerSection*(j-1) + (i-1)*stepsPerSection*sections):(j*stepsPerSection + (i-1)*stepsPerSection*sections)] = rev(c(Ntubes_row,NA))[i]
            } else {
              Tube_vector[(1+stepsPerSection*(j-1) + (i-1)*stepsPerSection*sections):(j*stepsPerSection + (i-1)*stepsPerSection*sections)] = c(Ntubes_row,NA)[i]
            }
          }
        }
        
        # Replace hs with film hs if condensation occurs
        for (i in seq(1,sections*stepsPerSection,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            for (k in seq(0,stepsPerSection-1)) { 
              if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && input$fluid == "vw") {
                h_s[(i+j+k)] = 0.943*((k_lw(Tcond[i+j+k])^3*dens_lw(Tcond[i+j+k])^2*9.81*Hvap[i+j+k])/(mu_lw(Tcond[i+j+k])*(input$tid+2*input$tt)*(Tcond[i+j+k]-Temp[(i+j+k+stepsPerSection*sections*(2*nrows+1))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
              } else if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && input$fluid == "vwe") {
                h_s[(i+j+k)] = 0.943*((k_lwe(x[(i+j+k)],Tcond[(i+j+k)])^3*dens_lwe(x[(i+j+k)],Tcond[(i+j+k)])^2*9.81*Hvap[i+j+k])/(mu_lwe(x[(i+j+k)],Tcond[(i+j+k)])*(input$tid+2*input$tt)*(Tcond[(i+j+k)]-Temp[(i+j+k+stepsPerSection*sections*(2*nrows+1))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
              } else if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && input$fluid == "vwd") {
                hw = 0.943*((k_lw(Tcond[(i+j+k)])^3*dens_lw(Tcond[(i+j+k)])^2*9.81*hvap_vw(Temp[i+j+k]))/(mu_lw(Tcond[(i+j+k)])*(input$tid+2*input$tt)*(Tcond[(i+j+k)]-Temp[(i+j+k+stepsPerSection*sections*(2*nrows+1))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
                hd = 0.943*((k_ld(Tcond[(i+j+k)])^3*dens_ld(Tcond[(i+j+k)])^2*9.81*hvap_vd(Temp[i+j+k]))/(mu_ld(Tcond[(i+j+k)])*(input$tid+2*input$tt)*(Tcond[(i+j+k)]-Temp[(i+j+k+stepsPerSection*sections*(2*nrows+1))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
                SAfracw = (x[i+j+k]/dens_lw(Temp[(i+j+k)]))/(x[i+j+k]/dens_lw(Temp[(i+j+k)])+(1-x[i+j+k])/dens_ld(Temp[(i+j+k)]))
                h_s[(i+j+k)] = SAfracw*hw + (1-SAfracw)*hd
              }
            }
          }
        }
        for (i in seq(1+stepsPerSection+sections*stepsPerSection,sections*stepsPerSection*2,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            for (k in seq(0,stepsPerSection-1)) {
              if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && input$fluid == "vw") {
                h_s[(i+j+k)] = 0.943*((k_lw(Tcond[i+j+k])^3*dens_lw(Tcond[i+j+k])^2*9.81*Hvap[i+j+k])/(mu_lw(Tcond[i+j+k])*(input$tid+2*input$tt)*(Tcond[i+j+k]-Temp[(i+j+k+stepsPerSection*sections*(2*nrows))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
              } else if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && input$fluid == "vwe") {
                h_s[(i+j+k)] = 0.943*((k_lwe(x[(i+j+k)],Tcond[(i+j+k)])^3*dens_lwe(x[(i+j+k)],Tcond[(i+j+k)])^2*9.81*Hvap[i+j+k])/(mu_lwe(x[(i+j+k)],Tcond[(i+j+k)])*(input$tid+2*input$tt)*(Tcond[(i+j+k)]-Temp[(i+j+k+stepsPerSection*sections*(2*nrows))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
              } else if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && input$fluid == "vwd") {
                hw = 0.943*((k_lw(Tcond[(i+j+k)])^3*dens_lw(Tcond[(i+j+k)])^2*9.81*hvap_vw(Temp[i+j+k]))/(mu_lw(Tcond[(i+j+k)])*(input$tid+2*input$tt)*(Tcond[(i+j+k)]-Temp[(i+j+k+stepsPerSection*sections*(2*nrows))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
                hd = 0.943*((k_ld(Tcond[(i+j+k)])^3*dens_ld(Tcond[(i+j+k)])^2*9.81*hvap_vd(Temp[i+j+k]))/(mu_ld(Tcond[(i+j+k)])*(input$tid+2*input$tt)*(Tcond[(i+j+k)]-Temp[(i+j+k+stepsPerSection*sections*(2*nrows))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
                SAfracw = (x[i+j+k]/dens_lw(Temp[(i+j+k)]))/(x[i+j+k]/dens_lw(Temp[(i+j+k)])+(1-x[i+j+k])/dens_ld(Temp[(i+j+k)]))
                h_s[(i+j+k)] = SAfracw*hw + (1-SAfracw)*hd
              }
            }
          }
        }
        
        # Set up a new ht vector and to match hs
        h_t_extended = c(h_t, rep(NA,sections*stepsPerSection))
        for (i in seq(2,sections, by = 2)) {
          for (j in (nrows+1):1) {
            if (j == 1){
              h_t_extended[((j-1)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-1)*sections*stepsPerSection + (i)*stepsPerSection)] = NA
            } else {
              h_t_extended[((j-1)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-1)*sections*stepsPerSection + (i)*stepsPerSection)] = 
                h_t[((j-2)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-2)*sections*stepsPerSection + (i)*stepsPerSection)]
            }
          }
        }
        
        tesths <<- h_s
        testhtext <<- h_t_extended
        
        
        # Sensible heat coefficents
        shellCoefficient1 = 1 - h_s*Tube_vector*SA_ot/(stepsPerSection*sections*2*Fs*Cp[1:(stepsPerSection*sections*(nrows+1))])
        shellCoefficient2 = -1 - h_s*Tube_vector*SA_ot/(stepsPerSection*sections*2*Fs*Cp[1:(stepsPerSection*sections*(nrows+1))])
        shellCoefficient3 = h_s*Tube_vector*SA_ot/(stepsPerSection*sections*Fs*Cp[1:(stepsPerSection*sections*(nrows+1))])
        tubeCoefficient1 = -1
        tubeCoefficient2 = 1 - h_t*sum(Ntubes_row)*SA_it/(sections*stepsPerSection*input$Ft*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))])
        tubeCoefficient3 = h_t*sum(Ntubes_row)*SA_it/(sections*stepsPerSection*input$Ft*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))])
        wallCoefficient1 = h_s*SA_ot/(2*(h_s*SA_ot + h_t_extended*SA_it))
        wallCoefficient2 = wallCoefficient1
        wallCoefficient3 = h_t_extended*SA_it/(h_s*SA_ot + h_t_extended*SA_it)
        wallCoefficient4 = -1
        
        
        # Condensation coefficent adjustments
        if (count %%2 == 0 && (input$fluid == "vw" || input$fluid == "vwe" || input$fluid == "vwd")) {
          for (i in 1:length(shellCoefficient1)){
            if (Temp[i] <= Tcond[i]) {
              shellCoefficient1[i] = 0
              shellCoefficient2[i] = 1
              shellCoefficient3[i] = 0
              
              rhs_vector[i] = Tcond[i] 
            }
          }
          for (j in seq(sections*stepsPerSection*nrows+stepsPerSection+1,sections*stepsPerSection*(nrows+1),by = stepsPerSection*2)) {
            if (Temp[j] <= Tcond[j]) {
              T_matrix[j:(stepsPerSection+j-1),(j-stepsPerSection):(j-1)] = matrix(data = 0,ncol = stepsPerSection, nrow = stepsPerSection)
            }
          }
          if (sections > 2){
            for (k in seq(1+stepsPerSection*2,sections*stepsPerSection,by = stepsPerSection*2)) {
              if (Temp[k] <= Tcond[k]) {
                T_matrix[k:(stepsPerSection+k-1),(k-stepsPerSection):(k-1)] = matrix(data = 0,ncol = stepsPerSection, nrow = stepsPerSection)
              }
            }
          }
        } else if (count %%2 == 1 && (input$fluid == "vw" || input$fluid == "vwe" || input$fluid == "vwd") && count != 1) {
          for (i in 1:length(tubeCoefficient2)){
            tubeCoefficient1 = 1
            tubeCoefficient2[i] = 0
            tubeCoefficient3[i] = 0
            rhs_vector[i+stepsPerSection*sections*(nrows+1)] = Temp[i+stepsPerSection*sections*(nrows+1)] 
          }
        }
        
        
        # Coefficent adjustments for total condensation
        # Downward condensation
        if (count %%2 == 0 && (input$fluid == "vw" || input$fluid == "vwe" || input$fluid == "vwd")) {
          for (i in seq(1,sections*stepsPerSection,by = stepsPerSection*2)) {
            for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
              for (k in seq(0,stepsPerSection-1)) { 
                if (Fs[i+j+k] == 0) {
                  tubeCoefficient2[i+j+k] = 1
                  tubeCoefficient3[i+j+k] = 0
                  wallCoefficient1[i+j+k] = 1/4
                  wallCoefficient2[i+j+k] = 1/4
                  wallCoefficient3[i+j+k] = 1/2
                }
              }
            }
          }
          # Upward condensation
          for (i in seq(1+stepsPerSection+sections*stepsPerSection,sections*stepsPerSection*2,by = stepsPerSection*2)) {
            for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
              for (k in seq(0,stepsPerSection-1)) {
                if (Fs[i+j+k] == 0) {
                  tubeCoefficient2[i+j+k-sections*stepsPerSection] = 1
                  tubeCoefficient3[i+j+k-sections*stepsPerSection] = 0
                  wallCoefficient1[i+j+k] = 1/4
                  wallCoefficient2[i+j+k] = 1/4
                  wallCoefficient3[i+j+k] = 1/2
                }
              }
            }
          }
          # Downward sensible
        } else if (count %%2 == 1 && (input$fluid == "vw" || input$fluid == "vwe" || input$fluid == "vwd")) {
          for (i in seq(1+sections*stepsPerSection,sections*stepsPerSection*2,by = stepsPerSection*2)) {
            for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
              for (k in seq(0,stepsPerSection-1)) { 
                if (Fs[i+j+k] == 0) {
                  shellCoefficient1[i+j+k-sections*stepsPerSection] = 0
                  shellCoefficient2[i+j+k-sections*stepsPerSection] = 1
                  shellCoefficient3[i+j+k-sections*stepsPerSection] = 0
                  rhs_vector[i+j+k] = Tcond[i+j+k]
                  wallCoefficient1[i+j+k-sections*stepsPerSection] = 1/4
                  wallCoefficient2[i+j+k-sections*stepsPerSection] = 1/4
                  wallCoefficient3[i+j+k-sections*stepsPerSection] = 1/2
                }
              }
            }
          }
          # Upward sensible
          for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
            for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
              for (k in seq(0,stepsPerSection-1)) {
                if (Fs[i+j+k] == 0) {
                  shellCoefficient1[i+j+k+sections*stepsPerSection] = 0
                  shellCoefficient2[i+j+k+sections*stepsPerSection] = 1
                  shellCoefficient3[i+j+k+sections*stepsPerSection] = 0
                  rhs_vector[i+j+k] = Tcond[i+j+k]
                  wallCoefficient1[i+j+k+sections*stepsPerSection] = 1/4
                  wallCoefficient2[i+j+k+sections*stepsPerSection] = 1/4
                  wallCoefficient3[i+j+k+sections*stepsPerSection] = 1/2
                }
              }
            }
          }
        }
        
        
        
        # Shell for the shell equation (downward flow)
        for (i in seq(1+stepsPerSection*sections,sections*stepsPerSection*2,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            for (k in seq(0,stepsPerSection-1,by = 1)) {
              if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && (input$fluid == "vw" || input$fluid == "vwe" || input$fluid == "vwd")) {
                T_matrix[(i+j+k),(i+j+k)] = shellCoefficient2[i+j+k]
              } else {
                T_matrix[(i+j+k),(i+j+k)] = shellCoefficient2[(i+j+k-sections*stepsPerSection)]
                T_matrix[(i+j+k),(i+j+k-sections*stepsPerSection)] = shellCoefficient1[(i+j+k-sections*stepsPerSection)]
              }
            }
          }
        }
        
        
        # Shell for the shell equation (upward flow)
        for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            for (k in seq(0,stepsPerSection-1,by = 1)) {
              if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && (input$fluid == "vw" || input$fluid == "vwe" || input$fluid == "vwd")) {
                T_matrix[(i+j+k),(i+j+k)] = shellCoefficient2[(i+j+k)]
              } else {
                T_matrix[(i+j+k),(i+j+k)] = shellCoefficient2[(i+j+k+sections*stepsPerSection)]
                T_matrix[(i+j+k),(i+j+k+sections*stepsPerSection)] = shellCoefficient1[(i+j+k+sections*stepsPerSection)]
              }
            }
          }
        }
        
        
        # Wall part for the shell equation (downward flow)
        for (i in seq(1+stepsPerSection*sections,sections*stepsPerSection*2,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            for (k in seq(0,stepsPerSection-1,by = 1)) {
              if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && (input$fluid == "vw" || input$fluid == "vwe" || input$fluid == "vwd")) {
                T_matrix[(i+j+k),(i+j+k+stepsPerSection*sections*(2*nrows))] = shellCoefficient3[i+j+k]
              } else {
                T_matrix[(i+j+k),(sections*stepsPerSection*2*length(Ntubes_row)+i+j+k)] = shellCoefficient3[(i+j+k-sections*stepsPerSection)]
              }
            }
          }
        }
        
        
        # Wall part for the shell equation (upward flow)
        for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            for (k in seq(0,stepsPerSection-1,by = 1)) {
              if (count %%2 == 0 && Temp[i+j+k] <= Tcond[i+j+k] && (input$fluid == "vw" || input$fluid == "vwe" || input$fluid == "vwd")) {
                T_matrix[(i+j+k),(i+j+k+stepsPerSection*sections*(2*nrows+1))] = shellCoefficient3[i+j+k]
              } else {
                T_matrix[(i+j+k),(sections*stepsPerSection*(2*length(Ntubes_row)+1)+i+j+k)] = shellCoefficient3[(i+j+k+sections*stepsPerSection)]  
              }
            }
          }
        }
        
        # Tube section of matrix
        for (i in seq(1+stepsPerSection*sections*(nrows+1),stepsPerSection*sections*(nrows+2)-1)) {
          for (j in seq(0,(nrows-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            T_matrix[(i+j),(i+j)] = tubeCoefficient1
            T_matrix[(i+j),(i+j+1)] = tubeCoefficient2[(i+j-sections*stepsPerSection*(nrows+1))]
            T_matrix[(i+j),(i+j+1+nrows*sections*stepsPerSection)] = tubeCoefficient3[(i+j-sections*stepsPerSection*(nrows+1))]
          }
        }
        
        
        # Wall for the wall equation (downward flow)
        for (i in seq(1+stepsPerSection*sections*(2*nrows+1),sections*stepsPerSection*(2*nrows+2),by = stepsPerSection*2)) {
          for (j in seq(0,(nrows-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j-stepsPerSection*sections*(2*nrows+1)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows+1))] = 
              diag(stepsPerSection)*wallCoefficient1[(i+j-stepsPerSection*sections*(2*nrows+1)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows+1))]
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j-stepsPerSection*sections*(2*nrows)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows))] = 
              diag(stepsPerSection)*wallCoefficient2[(i+j-stepsPerSection*sections*(2*nrows+1)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows+1))]
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j-stepsPerSection*sections*(nrows)):(stepsPerSection+i+j-1-stepsPerSection*sections*(nrows))] = 
              diag(stepsPerSection)*wallCoefficient3[(i+j-stepsPerSection*sections*(2*nrows+1)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows+1))]
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = 
              diag(stepsPerSection)*wallCoefficient4
            
          }
        }
        
        # Wall for the wall equation (upward flow)
        for (i in seq(1+stepsPerSection*sections*(2*nrows+1)+stepsPerSection,sections*stepsPerSection*(2*nrows+2),by = stepsPerSection*2)) {
          for (j in seq(0,(nrows-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j-stepsPerSection*sections*(2*nrows+1)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows+1))] = 
              diag(stepsPerSection)*wallCoefficient1[(i+j-stepsPerSection*sections*(2*nrows)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows))]
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j-stepsPerSection*sections*(2*nrows)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows))] = 
              diag(stepsPerSection)*wallCoefficient2[(i+j-stepsPerSection*sections*(2*nrows)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows))]
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j-stepsPerSection*sections*(nrows)):(stepsPerSection+i+j-1-stepsPerSection*sections*(nrows))] = 
              diag(stepsPerSection)*wallCoefficient3[(i+j-stepsPerSection*sections*(2*nrows)):(stepsPerSection+i+j-1-stepsPerSection*sections*(2*nrows))]
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = 
              diag(stepsPerSection)*wallCoefficient4
            
          }
        }
        
        testmatrix <<- T_matrix
        testb <<- rhs_vector
        testshell1 <<-shellCoefficient1
        testshell2 <<-shellCoefficient2
        testwall1 <<-wallCoefficient1
        
        # Solve matrix problem using LU decomposition
        x_vector = lusys(T_matrix,rhs_vector)
        
        testxvector <<- x_vector
        
        for (i in 1:(stepsPerSection*sections*(nrows+1))) {
          if (input$fluid == "vw") {
            if (x_vector[i] < Tcond[i]) {
              x_vector[i] = Tcond[i]
            }
          } else if (input$fluid == "vwe" || input$fluid == "vwd") {
            if (x_vector[i] < Tcond[i]) {
              x_vector[i] = Tcond[i]
            }
          }
        }
        Temp = x_vector
        
        testtemp <<- Temp
        
        # Error calc
        if (count %%2 == 1) {
          
          err = sqrt(mean(((x_vector - SensTempOld)/SensTempOld)^2))*100
          print(err)
          err = err+100
          
          SensTempOld = x_vector
        } else {
          
          err = sqrt(mean(((x_vector - CondTempOld)/CondTempOld)^2))*100
          print(err)
          
          CondTempOld = x_vector
        }
        
        if (count %%2 == 0){
          Q = (input$Ft/sum(Ntubes_row))*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*(
            Temp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))] - 
              Temp[(stepsPerSection*sections*(nrows+1)+2):(stepsPerSection*sections*(2*nrows+1)+1)])
          for (i in 1:nrows) {
            Q[stepsPerSection*sections*i] = Q[stepsPerSection*sections*i-1]
          }
          
          
          
          if (input$fluid == "vw" || input$fluid == "vwe" || input$fluid == "vwd") {
            
            # Set up a new Q vector
            Q_extended = c(Q, rep(NA,sections*stepsPerSection))
            for (i in seq(2,sections, by = 2)) {
              for (j in (nrows+1):1) {
                if (j == 1){
                  Q_extended[((j-1)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-1)*sections*stepsPerSection + (i)*stepsPerSection)] = NA
                } else {
                  Q_extended[((j-1)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-1)*sections*stepsPerSection + (i)*stepsPerSection)] = 
                    Q[((j-2)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-2)*sections*stepsPerSection + (i)*stepsPerSection)]
                }
              }
            }
            Q_extended = Q_extended*Tube_vector
            if (input$fluid == "vw") {
              for (i in 1:(stepsPerSection*sections)) {
                for (j in seq(0,(nrows-1)*sections*stepsPerSection, by = sections*stepsPerSection)) {
                  if (ceiling(i/stepsPerSection) %%2 == 1 && Temp[i+j] <= Tcond[i+j]) {
                    Mcond[i+j] = Q_extended[i+j]/(Hvap[i+j]+cp_lw((Tcond[i+j]+Temp[i+j+(2*nrows+1)*sections*stepsPerSection])/2)*(Tcond[i+j]-Temp[i+j+(2*nrows+1)*sections*stepsPerSection]))
                  } else if (ceiling(i/stepsPerSection) %%2 == 0 && Temp[i+j+sections*stepsPerSection] <= Tcond[i+j]) {
                    Mcond[i+j+sections*stepsPerSection] = Q_extended[i+j+sections*stepsPerSection]/(Hvap[i+j+sections*stepsPerSection]+cp_lw((Tcond[i+j]+Temp[i+j+(2*nrows+1)*sections*stepsPerSection])/2)*(Tcond[i+j]-Temp[i+j+(2*nrows+1)*sections*stepsPerSection]))
                  } else if (ceiling(i/stepsPerSection) %%2 == 1) {
                    Mcond[i+j] = 0
                  } else {
                    Mcond[i+j+sections*stepsPerSection] = 0
                  }
                }
              }
            } else if (input$fluid == "vwe") {
              for (i in 1:(stepsPerSection*sections)) {
                for (j in seq(0,(nrows-1)*sections*stepsPerSection, by = sections*stepsPerSection)) {
                  if (ceiling(i/stepsPerSection) %%2 == 1 && Temp[i+j] <= Tcond[i+j]) {
                    Mcond[i+j] = Q_extended[i+j]/(Hvap[i+j]+cp_lwe(x[i+j],(Tcond[i+j]+Temp[i+j+(2*nrows+1)*sections*stepsPerSection])/2)*(Tcond[i+j]-Temp[i+j+(2*nrows+1)*sections*stepsPerSection]))
                  } else if (ceiling(i/stepsPerSection) %%2 == 0 && Temp[i+j+sections*stepsPerSection] <= Tcond[i+j]) {
                    Mcond[i+j+sections*stepsPerSection] = Q_extended[i+j+sections*stepsPerSection]/(Hvap[i+j+sections*stepsPerSection]+cp_lwe(x[i+j],(Tcond[i+j]+Temp[i+j+(2*nrows+1)*sections*stepsPerSection])/2)*(Tcond[i+j]-Temp[i+j+(2*nrows+1)*sections*stepsPerSection]))
                  } else if (ceiling(i/stepsPerSection) %%2 == 1) {
                    Mcond[i+j] = 0
                  } else {
                    Mcond[i+j+sections*stepsPerSection] = 0
                  }
                }
              }
            } else if (input$fluid == "vwd") {
              for (i in 1:(stepsPerSection*sections)) {
                for (j in seq(0,(nrows-1)*sections*stepsPerSection, by = sections*stepsPerSection)) {
                  if (ceiling(i/stepsPerSection) %%2 == 1 && Temp[i+j] <= Tcond[i+j]) {
                    Mcond[i+j] = Q_extended[i+j]/(Hvap[i+j]+cp_lwd(x[i+j],(Tcond[i+j]+Temp[i+j+(2*nrows+1)*sections*stepsPerSection])/2)*(Tcond[i+j]-Temp[i+j+(2*nrows+1)*sections*stepsPerSection]))
                  } else if (ceiling(i/stepsPerSection) %%2 == 0 && Temp[i+j+sections*stepsPerSection] <= Tcond[i+j]) {
                    Mcond[i+j+sections*stepsPerSection] = Q_extended[i+j+sections*stepsPerSection]/(Hvap[i+j+sections*stepsPerSection]+cp_lwd(x[i+j],(Tcond[i+j]+Temp[i+j+(2*nrows+1)*sections*stepsPerSection])/2)*(Tcond[i+j]-Temp[i+j+(2*nrows+1)*sections*stepsPerSection]))
                  } else if (ceiling(i/stepsPerSection) %%2 == 1) {
                    Mcond[i+j] = 0
                  } else {
                    Mcond[i+j+sections*stepsPerSection] = 0
                  }
                }
              }
            }
            
            for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
              
              Mcond[i:(i+stepsPerSection-1)] = 0
              
            }
            for (i in seq(sections*stepsPerSection*nrows+1,sections*stepsPerSection*(nrows+1),by = stepsPerSection*2)) {
              
              Mcond[i:(i+stepsPerSection-1)] = 0
              
            }
            
            # Calculating cumulative mass loss due to cond.
            p = c(rep(1,stepsPerSection))
            checkFs = c(rep(1,stepsPerSection))
            for (i in 1:sections) {
              for (j in 1:(nrows+1)) {
                if (i %% 2 == 0) {
                  k = nrows + 2 - j
                  if (k == nrows + 1){
                    Fs[((i-1)*stepsPerSection +1+nrows*stepsPerSection*sections):(i*stepsPerSection+nrows*stepsPerSection*sections)] = 
                      mean(Fs[((i-2)*stepsPerSection +1+nrows*stepsPerSection*sections):((i-1)*stepsPerSection+nrows*stepsPerSection*sections)])
                    if (Fs[((i-1)*stepsPerSection +1+nrows*stepsPerSection*sections)] == 0) {
                      Mcond[((i-1)*stepsPerSection +1+nrows*stepsPerSection*sections):(i*stepsPerSection+nrows*stepsPerSection*sections)] = 0
                    }
                  } else {
                    Fs[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = 
                      Fs[((i-1)*stepsPerSection +1+k*stepsPerSection*sections):(i*stepsPerSection+k*stepsPerSection*sections)] -
                      Mcond[((i-1)*stepsPerSection +1+k*stepsPerSection*sections):(i*stepsPerSection+k*stepsPerSection*sections)]
                    for (m in 0:(stepsPerSection-1)) {
                      if (Fs[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections)+m] < 0) {
                        Mcond[((i-1)*stepsPerSection +1+k*stepsPerSection*sections)+m] = Fs[((i-1)*stepsPerSection +1+k*stepsPerSection*sections)+m]
                        Fs[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections)+m] = 0
                      }
                    }
                  }
                  
                  for (m in 0:(stepsPerSection-1)) {
                    if (Fs_old[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections)+m] == 0) {
                      Fs[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections)+m] = 0
                      if (p[m+1] == 1) {
                        Mcond[((i-1)*stepsPerSection+1+k*stepsPerSection*sections)+m] = Fs[((i-1)*stepsPerSection +1+k*stepsPerSection*sections)+m]
                      }
                      p[m+1] = p[m+1]+1
                    }
                    
                    if (Fs[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections)+m] == 0) {
                      checkFs[m+1] = checkFs[m+1]+1
                    }
                    if (checkFs[m+1] > 1 && k != nrows + 1) {
                      Fs[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections)+m] = 0
                    } else if (checkFs[m+1] > 1 && k == nrows + 1 && mean(Fs[((i-1)*stepsPerSection +1+nrows*stepsPerSection*sections):(i*stepsPerSection+nrows*stepsPerSection*sections)]) > 0) {
                      checkFs = c(rep(1,stepsPerSection))
                    }
                    
                  }
                  
                } else {
                  if (j == 1 && i != 1){
                    Fs[((i-1)*stepsPerSection +1):(i*stepsPerSection)] = mean(Fs[((i-2)*stepsPerSection +1):((i-1)*stepsPerSection)])
                    if (Fs[((i-1)*stepsPerSection +1)] == 0) {
                      Mcond[((i-1)*stepsPerSection +1):(i*stepsPerSection)] = 0
                    }
                  } else if (j == 1 && i == 1) {
                    Fs[1:stepsPerSection] = input$Fs/stepsPerSection
                  } else {
                    Fs[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = 
                      Fs[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections):(i*stepsPerSection+(j-2)*stepsPerSection*sections)] -
                      Mcond[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections):(i*stepsPerSection+(j-2)*stepsPerSection*sections)] 
                    
                    for (m in 0:(stepsPerSection-1)) {
                      if (Fs[((i-1)*stepsPerSection+1+(j-1)*stepsPerSection*sections)+m] <0) {
                        Mcond[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections)+m] = Fs[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections)+m]
                        Fs[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections)+m] = 0
                      }
                    }
                  }
                  for (m in 0:(stepsPerSection-1)) {
                    if (Fs_old[((i-1)*stepsPerSection+1+(j-1)*stepsPerSection*sections)+m] == 0) {
                      Fs[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections)+m] = 0
                      if (p[m+1] == 1) {
                        Mcond[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections)+m] = Fs[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections)+m]
                      }
                      p[m+1] = p[m+1]+1
                    }
                    
                    if (Fs[((i-1)*stepsPerSection+1+(j-1)*stepsPerSection*sections)+m] == 0) {
                      checkFs[m+1] = checkFs[m+1]+1
                    }
                    if (checkFs[m+1] > 1 && j != 1) {
                      Fs[((i-1)*stepsPerSection+1+(j-1)*stepsPerSection*sections)+m] = 0
                    } else if (checkFs[m+1] > 1 && j == 1 && mean(Fs[((i-2)*stepsPerSection +1):((i-1)*stepsPerSection)]) > 0) {
                      checkFs = c(rep(1,stepsPerSection))
                    }
                    
                    
                    
                  }
                }
              }
            }
            
            
            # # Replace negative Fs with 0
            # Fs[Fs<0] = 0
            # 
            # Fs_old = Fs
            # 
            # # Put in Fs = 0 where the flow rate stops changing
            # for (i in seq(1+sections*stepsPerSection,sections*stepsPerSection*2,by = stepsPerSection*2)) {
            #   for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            #     for (k in seq(0,stepsPerSection-1)) {
            #       if (Temp[i+j+k] == Temp[i+j+k-sections*stepsPerSection] && Fs_old[i+j+k] == Fs_old[i+j+k-sections*stepsPerSection]) {
            #         Fs[i+j+k] = 0
            #       }
            #     }
            #   }
            # }
            # for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
            #   for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            #     for (k in seq(0,stepsPerSection-1)) {
            #       if (Temp[i+j+k] == Temp[i+j+k+sections*stepsPerSection] && Fs_old[i+j+k] == Fs_old[i+j+k+sections*stepsPerSection]) {
            #         Fs[i+j+k] = 0
            #       }
            #     }
            #   }
            # }
            # for (i in seq((stepsPerSection+1),sections*stepsPerSection,by = stepsPerSection)) {
            #   if (ceiling(i/stepsPerSection) %% 2 == 1) {
            #     Fs[i:(i+stepsPerSection-1)] = mean(Fs[(i-stepsPerSection):(i-1)])
            #   } else {
            #     Fs[(i+sections*stepsPerSection*nrows):(i+stepsPerSection-1+sections*stepsPerSection*nrows)] = mean(Fs[(i-stepsPerSection+sections*stepsPerSection*nrows):(i-1+sections*stepsPerSection*nrows)])
            #   }
            # }
            tesFsold <<- Fs_old
            tesFs <<- Fs
            testQ <<- Q
            
            testQEXT <<- Q_extended
            testcond <<- Mcond
          }
        }
        Fs_old = Fs
        
        # Updating x and y for mixtures
        checky = c(rep(1,stepsPerSection))
        if (input$fluid == "vwe" || input$fluid == "vwd") {
          for (i in 1:sections) {
            for (j in 1:(nrows+1)) {
              if (i %% 2 == 0) {
                k = nrows + 2 - j
                if (k == nrows + 1){
                  y[((i-1)*stepsPerSection +1+nrows*stepsPerSection*sections):(i*stepsPerSection+nrows*stepsPerSection*sections)] = 
                    mean(y[((i-2)*stepsPerSection +1+nrows*stepsPerSection*sections):((i-1)*stepsPerSection+nrows*stepsPerSection*sections)])
                } else if (any(Fs[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)]<= 0)) {
                 y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = 
                    yfinal
                } else {
                 
                  y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = 
                    (Fs[((i-1)*stepsPerSection +1+k*stepsPerSection*sections):(i*stepsPerSection+k*stepsPerSection*sections)]*
                       y[((i-1)*stepsPerSection +1+k*stepsPerSection*sections):(i*stepsPerSection+k*stepsPerSection*sections)]-
                       Mcond[((i-1)*stepsPerSection +1+k*stepsPerSection*sections):(i*stepsPerSection+k*stepsPerSection*sections)]*
                       x[((i-1)*stepsPerSection +1+k*stepsPerSection*sections):(i*stepsPerSection+k*stepsPerSection*sections)])/
                    Fs[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)]
                  
                  if (input$fluid == "vwe") {
                    if (any(y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] <= yfinal)) {
                      y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = 
                        yfinal
                    }
                  } else {
                    if (input$fracw < 0.587583) {
                      if (any(y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] >= yfinal)) {
                        y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = 
                          yfinal
                      } 
                    } else {
                      if (any(y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] <= yfinal)) {
                        y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = 
                          yfinal
                      } 
                    }
                  }
                }
                
                
                
                if (input$fluid == "vwe") {
                  if (any(y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] >= input$fracw)) {
                    y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = input$fracw
                  }
                } else {
                  if (input$fracw < 0.587583) {
                    if (any(y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] <= input$fracw)) {
                      y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = 
                        input$fracw
                    } 
                  } else {
                    if (any(y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] >= input$fracw)) {
                      y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = 
                        input$fracw
                    } 
                  }
                }
                
                #Fix all compositions past the azeotrope to the azeotrope
                for (m in 0:(stepsPerSection-1)) {
                  if (y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections)+m] == yfinal) {
                    checky[m+1]=checky[m+1]+1
                  }
                  if (checky[m+1] > 1) {
                    y[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections)+m] = yfinal
                  }
                }
                
                
                
              } else {
                if (j == 1 && i != 1){
                  y[((i-1)*stepsPerSection +1):(i*stepsPerSection)] = mean(y[((i-2)*stepsPerSection +1):((i-1)*stepsPerSection)])
                } else if (j == 1 && i == 1) {
                  y[1:stepsPerSection] = input$fracw
                } else if (any(Fs[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)]<= 0)) {
                  y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = 
                    yfinal
               } else {
                 y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = 
                    (Fs[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections):(i*stepsPerSection+(j-2)*stepsPerSection*sections)]*
                       y[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections):(i*stepsPerSection+(j-2)*stepsPerSection*sections)]-
                       Mcond[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections):(i*stepsPerSection+(j-2)*stepsPerSection*sections)]*
                       x[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections):(i*stepsPerSection+(j-2)*stepsPerSection*sections)])/
                    Fs[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)]
                  

                 if (input$fluid == "vwe") {
                   if (any(y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] <= yfinal)) {
                     y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = 
                       yfinal
                   }
                 } else {
                   if (input$fracw < 0.587583) {
                     if (any(y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] >= yfinal)) {
                       y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = 
                         yfinal
                     } 
                   } else {
                     if (any(y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] <= yfinal)) {
                       y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = 
                         yfinal
                     } 
                   }
                 }
                }
                
                
                
                
                if (input$fluid == "vwe") {
                  if (any(y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] >= input$fracw)) {
                    y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = input$fracw
                  }
                } else {
                  if (input$fracw < 0.587583) {
                    if (any(y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] <= input$fracw)) {
                      y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = input$fracw
                    } 
                  } else {
                    if (any(y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] >= input$fracw)) {
                      y[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = input$fracw
                    } 
                  }
                }
                
                #Fix all compositions past the azeotrope to the azeotrope
                for (m in 0:(stepsPerSection-1)) {
                  if (y[((i-1)*stepsPerSection+1+(j-1)*stepsPerSection*sections)+m] == yfinal) {
                    checky[m+1]=checky[m+1]+1
                  }
                  if (checky[m+1] > 1) {
                    y[((i-1)*stepsPerSection+1+(j-1)*stepsPerSection*sections)+m] = yfinal
                  }
                }
                
              }
            }
          }
        
        Tcond = Tcondy(y)
        
        if (input$fluid == "vwd") {
            if (input$fracw < 0.587583) {
              x = nleqslv(rep(0.005,times = length(x)),Tcondx)$x
            }
            else {
              x = nleqslv(rep(0.995,times = length(x)),Tcondx)$x
            }

        } else {
            x = nleqslv(x,Tcondx)$x            
        }
        
          
        if (input$fluid == "vwe") {
          for (i in 1:length(x)) {
            if (y[i] == yfinal || x[i] < yfinal) {
              x[i] = yfinal
            }
          }
        } else {
          
          if (input$fracw < 0.587583) {
            for (i in 1:length(x)) {
              if (y[i] == yfinal || x[i] > yfinal) {
                x[i] = yfinal
              }
            }
          } else {
            for (i in 1:length(x)) {
              if (y[i] == yfinal || x[i] < yfinal) {
                x[i] = yfinal
              }
            }
          }
        }
          
        for (i in 1:(sections*stepsPerSection*(nrows+1))) {
          if (Temp[i] <= Tcond[i]) {
            Temp[i] = Tcondy(y[i])
          }
        }
      }
        
        count = count + 1
      }
      
      
      # Initialize heatmap data matrix without baffle columns
      heatmapdata = matrix(ncol = stepsPerSection*sections, nrow = 2*length(Ntubes_row)+1)
      
      # Insert tube temperatures
      for (i in 1:(length(Ntubes_row)+1)) {
        
        
        heatmapdata[2*i - 1,1:(stepsPerSection*sections)] = x_vector[((i-1)*stepsPerSection*sections + 1):(i*stepsPerSection*sections)]
        
      }
      
      for (i in 1:length(Ntubes_row)) {
        heatmapdata[2*i,1:(stepsPerSection*sections)] = x_vector[(i*stepsPerSection*sections + sections*stepsPerSection*length(Ntubes_row) + 1):(i*stepsPerSection*sections + sections*stepsPerSection*(length(Ntubes_row) + 1))]
        
      }
      
      
      heatmapdata2 <<- matrix(ncol = stepsPerSection*sections,nrow = 0)
      for (i in 1:nrow(heatmapdata)) {
        if (i %%2 == 0) {
          heatmapdata2 = rbind(heatmapdata2,heatmapdata[i,],heatmapdata[i,],heatmapdata[i,])
          
          
        } else {
          heatmapdata2 = rbind(heatmapdata2,heatmapdata[i,],heatmapdata[i,])
          
          
        }
      }
      
      
      heatmapdata2 <<- heatmapdata2 # save as global variable for debugging
      
      
      tubelines = list()
      counter = 0
      # Create border lines for tubes in heatmap
      for (i in 1:nrows) {
        counter = counter + 2
        tubelines = c(tubelines,
                      list(list(type = 'line', xref = "x", yref = "y",
                                x0 = -0.5, x1 = sections*stepsPerSection-0.5,
                                y0 = counter-0.5, y1 = counter-0.5,
                                line = list(color = "black", width = 2)
                      )
                      )
        )
        counter = counter + 3
        tubelines = c(tubelines,
                      list(list(type = 'line', xref = "x", yref = "y",
                                x0 = -0.5, x1 = sections*stepsPerSection-0.5,
                                y0 = counter-0.5, y1 = counter-0.5,
                                line = list(color = "black", width = 2)
                      )
                      )
        )
      
      }
      
      # Create lines for baffles in heatmap
      baffleline = list(
        type = "line",
        line = list(color = "black", width = 5),
        xref = "x",
        yref = "y"
      )
      
      # Store all line objects in one variable
      lines = list()
      for (i in 1:(sections-1)) {
        baffleline[["x0"]] = i*stepsPerSection-0.5
        baffleline[["x1"]] = i*stepsPerSection-0.5
        if ((i %% 2) == 0) {
          baffleline[c("y0", "y1")] = c(0.5,nrow(heatmapdata2)-0.5)
        } else {
          baffleline[c("y0", "y1")] = c(-0.5,nrow(heatmapdata2)-1.5)
        }
        lines = c(list(baffleline),lines )
      }
      
      lines = c(tubelines,lines)
      
      heatmaptext = list()
      for (i in 1:nrow(heatmapdata2)) {
        heatmaptext[[i]] = paste(round(heatmapdata2[i,1:(sections*stepsPerSection)],digits = 2), "K", sep = " ")
      }
      
      
      # Create arrows
      arrows = list(
        x = rep(ncol(heatmapdata2)-stepsPerSection,nrows),
        y = seq(3,3+(nrows-1)*5,5),
        text = "",
        xref = "x",
        yref = "y",
        arrowcolor = "white",
        arrowwidth = 2,
        showarrow = TRUE,
        arrowhead = 1,
        axref = "x",
        ayref = "y",
        ax = rep(ncol(heatmapdata2) - 1,nrows),
        ay = seq(3,3+(nrows-1)*5,5)
      )
      
      
      # Create plot and output
      output$plot1 = renderPlotly({
        
        # dotted cross-section line
        csline = list(
          type = "line",
          line = list(width = 5,dash="dashdot", color = "#FFFFFF"),
          xref = "x",
          yref = "y",
          x0 = (input$slider/(0.5*input$L/(sections*stepsPerSection))-1)/2,
          x1 = (input$slider/(0.5*input$L/(sections*stepsPerSection))-1)/2,
          y0 = -1,
          y1 = nrow(heatmapdata2)
        )
        
        
        plot_ly(z = heatmapdata2, type = "heatmap", hoverinfo = 'text', 
                text = heatmaptext,
                colorscale = list(c(0, "rgb(63,76,217)"), list(1, "#ed4a45")), 
                colorbar = list(len = 1, title = list(text = "Temperature (K)",side = "right"))) %>%
          add_annotations(
            x = stepsPerSection/2-0.5,
            y = nrow(heatmapdata2)/3,
            text = "",
            xref = "x",
            yref = "y",
            arrowcolor = "white",
            arrowwidth = 4,
            showarrow = TRUE,
            arrowhead = 1,
            axref = "x",
            ayref = "y",
            ax = stepsPerSection/2-0.5,
            ay = 0
          )  %>%
          config(displayModeBar = F) %>%
          layout(
            shapes = c(lines,list(csline)),
            annotations = arrows,
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            xaxis = list(
              ticks = "",
              showticklabels = F,
              showgrid = F,
              zeroline = F,
              fixedrange = T
            ),
            yaxis = list(
              ticks = "",
              showticklabels = F,
              showgrid = F,
              zeroline = F,
              fixedrange = T
            )
          )
      })
      
      
      # Set-up and output circle plot
      output$plot2 =  renderPlotly({
        x = seq(-1,1,1/100)
        circlefun = function(x) {
          
          y = sqrt(1-x^2)
          return(y)
        }
        y = circlefun(x)
        
        x = c(x,rev(x))
        y = c(y,-y)
        
        
        normpitchx = pitch_x/(input$sid/2)
        normpitchy = pitch_y/(input$sid/2)
        
        
        ypos = c(seq(-1+0.5*normpitchy,-1+(nrows-0.5)*normpitchy,normpitchy),1)
        
        
        plot2 = plot_ly(type = "scatter", mode = "lines", showlegend = F) %>%
          add_trace(x = c(-1,1), y = c(-1,-1), line = list(color = "rgba(0,0,0,0)")) %>%
          config(displayModeBar = F) %>%
          layout(
            xaxis = list(
              zeroline = FALSE,
              showline = FALSE,
              showticklabels = FALSE,
              showgrid = FALSE,
              range = c(-1.1,1.1)
            ),
            yaxis = list(
              zeroline = FALSE,
              showline = FALSE,
              showticklabels = FALSE,
              showgrid = FALSE,
              range = c(-1.1,1.1)
            ),
            paper_bgcolor = "rgba(0,0,0,0)",
            plot_bgcolor = "rgba(0,0,0,0)"
          )
        
        
        shellcolorvalues = round(heatmapdata[seq(1,2*nrows+1,2),]) - floor(min(heatmapdata)) + 1
        tubecolorvalues = round(heatmapdata[seq(2,2*nrows,2),]) - floor(min(heatmapdata)) + 1
        colorrange = colorRampPalette(c("#3f4cd9","#ed4a45")) (ceiling(max(heatmapdata))-floor(min(heatmapdata)) + 1)
        
        for (i in 1:(nrows+1)) {
          plot2 = add_trace(plot2,type = "scatter", mode = "lines",x = c(-1,1),y = c(ypos[i],ypos[i]), 
                            fill = "tonextx", fillcolor = colorrange[shellcolorvalues[i,input$slider/(0.5*input$L/(sections*stepsPerSection))/2+0.5]], line = list(color = "rgba(0,0,0,0)"))
        }
        
        plot2 = add_trace(plot2,x = x*2, y = y*2)
        plot2 = add_trace(plot2,x = x, y = y, fill = "tonextx", fillcolor = "#E2E2E2", line = list(color = "black"))
        
        normchord = chord_length/(input$sid/2)
        normtid = input$tid/(input$sid)*454 # 454px is width of the shell in our plot 
        
        if (input$config == "square"){
          if (sum(Ntubes_row %%2) == 0 || sum(Ntubes_row %%2) == nrows) {
            xpos = seq(-min(normchord[Ntubes_row %in% max(Ntubes_row)])/2 + 0.5*normpitchx,min(normchord[Ntubes_row %in% max(Ntubes_row)])/2 - 0.5*normpitchx,normpitchx)
            for (i in 1:nrows) {
              if (length(xpos) == Ntubes_row[i]) {
                xposi = xpos
              } else{
                xposi = xpos[-c((1:(max(Ntubes_row)-Ntubes_row[i])/2),length(xpos):(length(xpos)-(max(Ntubes_row)-Ntubes_row[i])/2 + 1))]
              }
              plot2 = add_trace(plot2, mode = "markers", x = xposi, y = rep(ypos[i],Ntubes_row[i]),
                                marker = list(size = normtid, line = list(color = "black",width = 2),
                                              color = colorrange[rep(tubecolorvalues[i,input$slider/(0.5*input$L/(sections*stepsPerSection))/2+0.5],Ntubes_row[i])]))
            }
            
          } else {
            
          }
          
          
        } else {
          if (sum(Ntubes_row %%2) == (nrows/2 + 0.5) || sum(Ntubes_row %%2) == (nrows/2 - 0.5)) {
            # x positions of rows alligned with largest row
            xpos1 = seq(-min(normchord[Ntubes_row %in% max(Ntubes_row)])/2 + 0.5*normpitchx,min(normchord[Ntubes_row %in% max(Ntubes_row)])/2 - 0.5*normpitchx,normpitchx)
            
            # x positions of rows not alligned with largest row
            xpos2 = xpos1[-1] - normpitchx/2
            
            for (i in 1:nrows) {
              if (length(xpos1) == Ntubes_row[i]) {
                xposi = xpos1
              } else if (abs(i-length(Ntubes_row)/2 + 0.5) %%2) {
                xposi = xpos1[-c((1:(max(Ntubes_row)-Ntubes_row[i])/2),length(xpos1):(length(xpos1)-(max(Ntubes_row)-Ntubes_row[i])/2 + 1))]
              } else if (length(xpos2) == Ntubes_row[i]) {
                xposi = xpos2
              } else {
                xposi = xpos2[-c((1:(max(Ntubes_row)-1-Ntubes_row[i])/2),length(xpos2):(length(xpos2)-(max(Ntubes_row)-1-Ntubes_row[i])/2 + 1))]
              }
              plot2 = add_trace(plot2, mode = "markers", x = xposi, y = rep(ypos[i],Ntubes_row[i]),
                                marker = list(size = normtid, line = list(color = "black",width = 2),
                                              color = colorrange[rep(tubecolorvalues[i,input$slider/(0.5*input$L/(sections*stepsPerSection))/2+0.5],Ntubes_row[i])]))
            }
            
          } else {
            
          }
          
        }
        plot2
        
        
        
      })
      
      # Update slider bar
      updateSliderInput(session, "slider", min = 0.5*input$L/(sections*stepsPerSection) , max = input$L-0.5*input$L/(sections*stepsPerSection), label = "X-Position (m)",
                        step = 0.5*input$L/(sections*stepsPerSection), value = 0.5*input$L/(sections*stepsPerSection))
      showElement("slider")
      showTab("tabs","Output")
      
      # Create and output results table
      if (input$fluid == "lw") {
        column1 = c("Output Shell Temperature (K)","Output Tube Temperature (K)", HTML("&Delta;T Shell (K)"), 
                    HTML("&Delta;T Tube (K)"), "Q (W)")
        
        if (sections %%2) {
          outputshelltemp = mean(heatmapdata[1,(stepsPerSection*(sections-1)+1):(stepsPerSection*sections)])
          outputtubetemp = sum(heatmapdata[seq(2,2*nrows,2),1]*Ntubes_row)/sum(Ntubes_row)
          
        } else {
          outputshelltemp = mean(heatmapdata[2*nrows+1,(stepsPerSection*(sections-1)+1):(stepsPerSection*sections)])
          outputtubetemp = sum(heatmapdata[seq(2,2*nrows,2),1]*Ntubes_row)/sum(Ntubes_row)
        }
        deltas = outputshelltemp - input$Tsi
        deltat = outputtubetemp - input$Tti  
        
        column2 = c(outputshelltemp,outputtubetemp,deltas,deltat, sum(Q_extended,na.rm =T))
        
        
      } else if (input$fluid == "vw") {
        
        column1 = c("Output Shell Temperature (K)","Output Tube Temperature (K)", HTML("&Delta;T Shell (K)"), 
                    HTML("&Delta;T Tube (K)"), "Mass Condensed (kg/s)", "% Condensed", "Q (W)")
        
        if (sections %%2) {
          outputshelltemp = mean(heatmapdata[1,(stepsPerSection*(sections-1)+1):(stepsPerSection*sections)])
          outputtubetemp = sum(heatmapdata[seq(2,2*nrows,2),1]*Ntubes_row)/sum(Ntubes_row)
          
        } else {
          outputshelltemp = mean(heatmapdata[2*nrows+1,(stepsPerSection*(sections-1)+1):(stepsPerSection*sections)])
          outputtubetemp = sum(heatmapdata[seq(2,2*nrows,2),1]*Ntubes_row)/sum(Ntubes_row)
        }
        deltas = outputshelltemp - input$Tsi
        deltat = outputtubetemp - input$Tti  
        
        column2 = c(outputshelltemp,outputtubetemp,deltas,deltat,sum(Mcond),sum(Mcond)/input$Fs*100, sum(Q_extended,na.rm =T))
        
      } else if (input$fluid == "vwe" || input$fluid == "vwd") {
        
        column1 = c("Output Shell Temperature (K)","Output Tube Temperature (K)", HTML("&Delta;T Shell (K)"), 
                    HTML("&Delta;T Tube (K)"), "Mass Condensed (kg/s)", "% Condensed", "Q (W)")
        
        if (sections %%2) {
          outputshelltemp = mean(heatmapdata[1,(stepsPerSection*(sections-1)+1):(stepsPerSection*sections)])
          outputtubetemp = sum(heatmapdata[seq(2,2*nrows,2),1]*Ntubes_row)/sum(Ntubes_row)
          
        } else {
          outputshelltemp = mean(heatmapdata[2*nrows+1,(stepsPerSection*(sections-1)+1):(stepsPerSection*sections)])
          outputtubetemp = sum(heatmapdata[seq(2,2*nrows,2),1]*Ntubes_row)/sum(Ntubes_row)
        }
        deltas = outputshelltemp - input$Tsi
        deltat = outputtubetemp - input$Tti  
        
        column2 = c(outputshelltemp,outputtubetemp,deltas,deltat,sum(Mcond),sum(Mcond)/input$Fs*100, sum(Q_extended,na.rm =T))
        
      }
      
      Results = t(data.frame(column2))
      colnames(Results) = column1
      
      output$resulttable = renderTable(Results,digits = 6, sanitize.text.function = function(x) x, align = 'c')
      
      updateTabsetPanel(session, "tabs", selected = "Output")
      
    })
    }
  })
  
  
}