server = function(input,output,session){
  hideElement("slider")
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
    dens_lw = function(x) {
      Tau = 1 - x/647.096
      y = 17.863 + 58.606*Tau^0.35 - 95.396*Tau^(2/3) + 213.89*Tau - 141.26*Tau^(4/3) # mol/dm^3
      y = y*18.015 # kg/m3
      return(y)
    }
    
    cp_lw = function(x) {
      y = 276370 - 2090.1*x + 8.125*x^2 - 0.014116*x^3 + 9.3701*10^(-6)*x^4 # J/(kmol K)
      y = y/18.015 # J/(kg K)
      return(y)
    }
    
    mu_lw = function(x) {
      y = exp(-52.843 + 3703.6/x + 5.866*log(x) -5.879*10^(-29)*T^10) # Pa s
      return(y)
    }
    
    k_lw = function(x) {
      y = -0.432 + 0.0057255*x - 0.000008078*x^2 + 1.861*10^(-9)*x^3 # W/(m K)
      return(y)
    }
    
    # Vapour water properties
    hvap_vw = function(x) {
      Tau = 1 - x/647.096
      y = 5.2053*10^7*(1-Tau)^(0.3199 - 0.212*Tau + 0.25795*Tau^2) # J/kmol
      y = y/18.015 # J/kg
      return(y)
    }
    
    cp_vw = function(x) {
      y = 0.33363*10^5 + 0.2697*10^5*((2.6105*10^3/x)/(sinh((2.6105*10^3)/x)))^2 + 0.08896*10^5*((1169/x)/(sinh((1169)/x)))^2 # J/(kmol K)
      y = y/18.015 # J/(kg K)
      return(y)
    }
    
    mu_vw = function(x) {
      y = 1.7096*10^-8*x^1.1146 # Pa s
      return(y)
    }
    
    k_vw = function(x) {
      y = 6.2041*10^-6*x^1.3973 # W/m K
      return(y)
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
      
      
      err = 1
      count = 1
      while (err > 0.001) { 
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
      
      
      err = 1
      count = 1
      while (err > 0.001) { 
        
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
        
        Tcond = 1730.63/(8.07131-log10(input$P/101325*760)) - 233.426 + 273.15
        for (i in 1:(stepsPerSection*sections*(nrows+1))) {
          if (Temp[i] < Tcond && input$fluid == "vw") {
            Temp[i] = Tcond
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
        Pr_w = Mu[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))]*Cp[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))]/K[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))]
        
        # Create Pr_w vector that matches the shell vector
        Pr_w_aligned = c(Pr_w,rep(NA,stepsPerSection*sections))
        for (i in seq(2,sections, by = 2)) {
          for (j in (nrows+1):1) {
            if (j == 1){
              Pr_w_aligned[((j-1)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-1)*sections*stepsPerSection + (i)*stepsPerSection)] = -999
            } else {
              Pr_w_aligned[((j-1)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-1)*sections*stepsPerSection + (i)*stepsPerSection)] = 
                Pr_w[((j-2)*sections*stepsPerSection + (i-1)*stepsPerSection + 1):((j-2)*sections*stepsPerSection + (i)*stepsPerSection)]
            }
          }
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
              if (count %%2 == 0 && Temp[i+j+k] <= Tcond && input$fluid == "vw") {
                h_s[(i+j+k)] = 0.943*((k_lw(Tcond)^3*dens_lw(Tcond)^2*9.81*Hvap[i+j+k])/(mu_lw(Tcond)*(input$tid+2*input$tt)*(Tcond-Temp[(i+j+k+stepsPerSection*sections*(2*nrows+1))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
              }
            }
          }
        }
        for (i in seq(1+stepsPerSection+sections*stepsPerSection,sections*stepsPerSection*2,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            for (k in seq(0,stepsPerSection-1)) {
              if (count %%2 == 0 && Temp[i+j+k] <= Tcond && input$fluid == "vw") {
                h_s[(i+j+k)] = 0.943*((k_lw(Tcond)^3*dens_lw(Tcond)^2*9.81*Hvap[i+j+k])/(mu_lw(Tcond)*(input$tid+2*input$tt)*(Tcond-Temp[(i+j+k+stepsPerSection*sections*(2*nrows))])))^0.25*(0.6+0.42*Tube_vector[i+j+k]^-0.25)
              }
            }
          }
        }
        
        
        tesths <<- h_s

        
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
        if (count %%2 == 0 && input$fluid == "vw") {
          for (i in 1:length(shellCoefficient1)){
            if (Temp[i] <= Tcond) {
              shellCoefficient1[i] = 0
              shellCoefficient2[i] = 1
              shellCoefficient3[i] = 0
              
              rhs_vector[i] = Tcond  
            }
          }
          for (j in seq(sections*stepsPerSection*nrows+stepsPerSection+1,sections*stepsPerSection*(nrows+1),by = stepsPerSection*2)) {
            if (Temp[j] <= Tcond) {
              T_matrix[j:(stepsPerSection+j-1),(j-stepsPerSection):(j-1)] = matrix(data = 0,ncol = stepsPerSection, nrow = stepsPerSection)
            }
          }
          if (sections > 2){
            for (k in seq(1+stepsPerSection*2,sections*stepsPerSection,by = stepsPerSection*2)) {
              if (Temp[k] <= Tcond) {
                T_matrix[k:(stepsPerSection+k-1),(k-stepsPerSection):(k-1)] = matrix(data = 0,ncol = stepsPerSection, nrow = stepsPerSection)
              }
            }
          }
        } else if (count %%2 == 1 && input$fluid == "vw" && count != 1) {
          for (i in 1:length(tubeCoefficient2)){
            tubeCoefficient1 = 1
            tubeCoefficient2[i] = 0
            tubeCoefficient3[i] = 0
            rhs_vector[i+stepsPerSection*sections*(nrows+1)] = Temp[i+stepsPerSection*sections*(nrows+1)] 
          }
        }
        
        
        # Coefficent adjustments for total condensation
        # Downward condensation
        if (count %%2 == 0 && input$fluid == "vw") {
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
        } else if (count %%2 == 1 && input$fluid == "vw") {
          for (i in seq(1+sections*stepsPerSection,sections*stepsPerSection*2,by = stepsPerSection*2)) {
            for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
              for (k in seq(0,stepsPerSection-1)) { 
                if (Fs[i+j+k] == 0) {
                  shellCoefficient1[i+j+k-sections*stepsPerSection] = 0
                  shellCoefficient2[i+j+k-sections*stepsPerSection] = 1
                  shellCoefficient3[i+j+k-sections*stepsPerSection] = 0
                  rhs_vector[i+j+k] = Tcond
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
                  rhs_vector[i+j+k] = Tcond
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
            if (count %%2 == 0 && Temp[i+j] <= Tcond && input$fluid == "vw") {
              T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = diag(stepsPerSection)*shellCoefficient2[(i+j):(stepsPerSection+i+j-1)]
            } else {
              T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = diag(stepsPerSection)*shellCoefficient2[(i+j-sections*stepsPerSection):(stepsPerSection+i+j-1-sections*stepsPerSection)]
              T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j-sections*stepsPerSection):(i+j-sections*stepsPerSection+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient1[(i+j-sections*stepsPerSection):(stepsPerSection+i+j-1-sections*stepsPerSection)]
            }
          }
        }
        
        
        # Shell for the shell equation (upward flow)
        for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            if (count %%2 == 0 && Temp[i+j] <= Tcond && input$fluid == "vw") {
              T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = diag(stepsPerSection)*shellCoefficient2[(i+j):(stepsPerSection+i+j-1)]
            } else {
              T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = diag(stepsPerSection)*shellCoefficient2[(i+j+sections*stepsPerSection):(stepsPerSection+i+j-1+sections*stepsPerSection)]
              T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j+sections*stepsPerSection):(i+j+sections*stepsPerSection+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient1[(i+j+sections*stepsPerSection):(stepsPerSection+i+j-1+sections*stepsPerSection)]
            }
          }
        }
        
        
        # Wall part for the shell equation (downward flow)
        for (i in seq(1+stepsPerSection*sections,sections*stepsPerSection*2,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            if (count %%2 == 0 && Temp[i+j] <= Tcond && input$fluid == "vw") {
              T_matrix[(i+j),(i+j+stepsPerSection*sections*(2*nrows))] = shellCoefficient3[i+j]
            } else {
              T_matrix[(i+j):(stepsPerSection+i+j-1),(sections*stepsPerSection*2*length(Ntubes_row)+i+j):(sections*stepsPerSection*2*length(Ntubes_row)+i+j+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient3[(i+j-sections*stepsPerSection):(stepsPerSection+i+j-1-sections*stepsPerSection)]
            }
          }
        }
        
        
        # Wall part for the shell equation (upward flow)
        for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            if (count %%2 == 0 && Temp[i+j] <= Tcond && input$fluid == "vw") {
              T_matrix[(i+j),(i+j+stepsPerSection*sections*(2*nrows+1))] = shellCoefficient3[i+j]
            } else {
              T_matrix[(i+j):(stepsPerSection+i+j-1),(sections*stepsPerSection*(2*length(Ntubes_row)+1)+i+j):(sections*stepsPerSection*(2*length(Ntubes_row)+1)+i+j+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient3[(i+j+sections*stepsPerSection):(stepsPerSection+i+j-1+sections*stepsPerSection)]  
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
        testwall1 <<-wallCoefficient1
        
        # Solve matrix problem using LU decomposition
        x_vector = lusys(T_matrix,rhs_vector)
        
        testxvector <<- x_vector
        
        for (i in 1:(stepsPerSection*sections*(nrows+1))) {
          if (x_vector[i] < Tcond && input$fluid == "vw") {
            x_vector[i] = Tcond
          }
        }
        Temp = x_vector
        
        testtemp <<- Temp
        
        # Error calc
        if (count %%2 == 1) {
          
          err = sqrt(mean(((x_vector - SensTempOld)/SensTempOld)^2))
          print(err)
          err = err+1
          
          SensTempOld = x_vector
        } else {
          
          err = sqrt(mean(((x_vector - CondTempOld)/CondTempOld)^2))
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
          
          
          
          if (input$fluid == "vw"){
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
            
            for (i in 1:(stepsPerSection*sections*(nrows+1))) {
              if (Temp[i] <= Tcond && input$fluid == "vw") {
                Mcond[i] = Q_extended[i]/Hvap[i]
              } else {
                Mcond[i] = 0
              }
            }
            
            for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
              
              Mcond[i:(i+stepsPerSection-1)] = 0
              
            }
            for (i in seq(sections*stepsPerSection*nrows+1,sections*stepsPerSection*(nrows+1),by = stepsPerSection*2)) {
              
              Mcond[i:(i+stepsPerSection-1)] = 0
              
            }
            
            
            # Calculating cumulative mass loss due to cond.
            for (i in 1:sections) {
              for (j in 1:(nrows+1)) {
                if (i %% 2 == 0) {
                  k = nrows + 2 - j
                  if (k == nrows + 1){
                    Fs[((i-1)*stepsPerSection +1+nrows*stepsPerSection*sections):(i*stepsPerSection+nrows*stepsPerSection*sections)] = 
                      mean(Fs[((i-2)*stepsPerSection +1+nrows*stepsPerSection*sections):((i-1)*stepsPerSection+nrows*stepsPerSection*sections)])
                  } else {
                    Fs[((i-1)*stepsPerSection+1+(k-1)*stepsPerSection*sections):(i*stepsPerSection+(k-1)*stepsPerSection*sections)] = 
                      Fs[((i-1)*stepsPerSection +1+k*stepsPerSection*sections):(i*stepsPerSection+k*stepsPerSection*sections)] -
                      Mcond[((i-1)*stepsPerSection +1+k*stepsPerSection*sections):(i*stepsPerSection+k*stepsPerSection*sections)] 
                  }
                } else {
                  if (j == 1 && i != 1){
                    Fs[((i-1)*stepsPerSection +1):(i*stepsPerSection)] = mean(Fs[((i-2)*stepsPerSection +1):((i-1)*stepsPerSection)])
                  } else if (j == 1 && i == 1) {
                    Fs[1:stepsPerSection] = input$Fs/stepsPerSection
                  } else {
                    Fs[((i-1)*stepsPerSection +1+(j-1)*stepsPerSection*sections):(i*stepsPerSection+(j-1)*stepsPerSection*sections)] = 
                      Fs[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections):(i*stepsPerSection+(j-2)*stepsPerSection*sections)] -
                      Mcond[((i-1)*stepsPerSection +1+(j-2)*stepsPerSection*sections):(i*stepsPerSection+(j-2)*stepsPerSection*sections)] 
                  }
                }
              }
            }
            
            # Replace negative Fs with 0
            Fs[Fs<0] = 0
            
            Fs_old = Fs
            
            # Put in Fs = 0 where the flow rate stops changing
            for (i in seq(1+sections*stepsPerSection,sections*stepsPerSection*2,by = stepsPerSection*2)) {
              for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
                for (k in seq(0,stepsPerSection-1)) {
                  if (Temp[i+j+k] == Temp[i+j+k-sections*stepsPerSection] && Fs_old[i+j+k] == Fs_old[i+j+k-sections*stepsPerSection]) {
                    Fs[i+j+k] = 0
                  }
                }
              }
            }
            for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
              for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
                for (k in seq(0,stepsPerSection-1)) {
                  if (Temp[i+j+k] == Temp[i+j+k+sections*stepsPerSection] && Fs_old[i+j+k] == Fs_old[i+j+k+sections*stepsPerSection]) {
                    Fs[i+j+k] = 0
                  }
                }
              }
            }
            for (i in seq((stepsPerSection+1),sections*stepsPerSection,by = stepsPerSection)) {
              if (ceiling(i/stepsPerSection) %% 2 == 1) {
                Fs[i:(i+stepsPerSection-1)] = mean(Fs[(i-stepsPerSection):(i-1)])
              } else {
                Fs[(i+sections*stepsPerSection*nrows):(i+stepsPerSection-1+sections*stepsPerSection*nrows)] = mean(Fs[(i-stepsPerSection+sections*stepsPerSection*nrows):(i-1+sections*stepsPerSection*nrows)])
              }
            }
            
            tesFs <<- Fs
            testQ <<- Q_extended
            testcond <<- Mcond
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
      
      
      heatmapdata <<- heatmapdata # save as global variable for debugging
      
      tubelines = list()
      # Create border lines for tubes in heatmap
      for (i in 1:(2*length(Ntubes_row))) {
        tubelines = c(tubelines,
                      list(list(type = 'line', xref = "x", yref = "y",
                                x0 = -0.5, x1 = sections*stepsPerSection-0.5,
                                y0 = i-0.5, y1 = i-0.5,
                                line = list(color = "black", width = 10)
                      )
                      )
        )
      }
      
      # Create lines for baffles in heatmap
      baffleline = list(
        type = "line",
        line = list(color = "black", width = 10),
        xref = "x",
        yref = "y"
      )
      
      # Store all line objects in one variable
      lines = list()
      for (i in 1:(sections-1)) {
        baffleline[["x0"]] = i*stepsPerSection-0.5
        baffleline[["x1"]] = i*stepsPerSection-0.5
        if ((i %% 2) == 0) {
          baffleline[c("y0", "y1")] = c(0,2*length(Ntubes_row)+0.5)
        } else {
          baffleline[c("y0", "y1")] = c(-0.5,2*length(Ntubes_row))
        }
        lines = c(lines, list(baffleline))
      }

      lines = c(lines,tubelines)

      heatmaptext = list()
      for (i in 1:nrow(heatmapdata)) {
        heatmaptext[[i]] = paste(round(heatmapdata[i,1:(sections*stepsPerSection)],digits = 2), "K", sep = " ")
      }
      
      
      # Create plot and output
      output$plot1 = renderPlotly({
        
        
        csline = list(
          type = "line",
          line = list(width = 5,dash="dashdot", color = "white"),
          xref = "x",
          yref = "y",
          x0 = round(input$slider*sections*stepsPerSection/(0.5*input$L))/2-0.5,
          x1 = round(input$slider*sections*stepsPerSection/(0.5*input$L))/2-0.5,
          y0 = -1,
          y1 = 2*nrows +1
        )
        
        
        plot_ly(z = heatmapdata, type = "heatmap", hoverinfo = 'text', 
                text = heatmaptext,
                colorscale = list(c(0, "rgb(0, 0, 255)"), list(1, "rgb(255, 0, 0)")), 
                colorbar = list(len = 1, title = list(text = "Temperature (K)",side = "right"))) %>%
          config(displayModeBar = F) %>%
          layout(
            shapes = c(lines,list(csline)),
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
        colorrange = colorRampPalette(c("#0000ff","#ff0000")) (ceiling(max(heatmapdata))-floor(min(heatmapdata)) + 1)
        
        for (i in 1:(nrows+1)) {
          plot2 = add_trace(plot2,type = "scatter", mode = "lines",x = c(-1,1),y = c(ypos[i],ypos[i]), 
                            fill = "tonextx", fillcolor = colorrange[shellcolorvalues[i,round(input$slider*sections*stepsPerSection/(0.5*input$L))/2+0.5]], line = list(color = "rgba(0,0,0,0)"))
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
                                              color = colorrange[rep(tubecolorvalues[i,round(input$slider*sections*stepsPerSection/(0.5*input$L))/2+0.5],Ntubes_row[i])]))
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
                                              color = colorrange[rep(tubecolorvalues[i,round(input$slider*sections*stepsPerSection/(0.5*input$L))/2+0.5],Ntubes_row[i])]))
            }
            
          } else {
            
          }
          
        }
        plot2
        
        
        
      })
      
      # Update slider bar
      updateSliderInput(session, "slider", min = 0.5*input$L/(sections*stepsPerSection) , max = input$L-0.5*input$L/(sections*stepsPerSection), label = "X-Position (m)",
                        step = input$L/(sections*stepsPerSection), value = 0.5*input$L/(sections*stepsPerSection))
      showElement("slider")


    })
    }
  })
  observeEvent(input$slider,{

  })
  
}