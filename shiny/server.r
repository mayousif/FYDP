server = function(input,output,session){
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
    Ntubes_row <<- as.numeric(unlist(strsplit(input$ntubesrow,split = ",")))
    
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
    
    
    # Double-pipe, no baffles --------------------------------------------------------------------------------------------  
    if (input$baffles == 0) {withProgress(value = 0.5, message = "Calculation in progress...",{
      
      # Cross sectional area of the shell (m2)
      CSA_s = pi/4*(input$sid^2-(input$tid+2*input$tt)^2)
      
      # Cross sectional area of the tube (m2)
      CSA_t = pi/4*(input$tid^2)
      
      # Hydraulic diameter shell/pipe
      dhs = input$sid - (input$tid + input$tt*2)
      dht = input$tid
      
      # Inital temperature estimate vector
      Temp <<- c(rep(input$Tsi,times = input$modelsteps),
                 rep(input$Tti,times = input$modelsteps),
                 rep(mean(c(input$Tti,input$Tsi)),times = input$modelsteps))
      
      err = 1
      while (err > 0.001) { 
        # Create matrix containing all calclated temp values
        T_matrix <<- matrix(data = 0, ncol = input$modelsteps*3,nrow = input$modelsteps*3)
  
        # Properties
        Dens <<- dens_lw(Temp)
        Cp <<- cp_lw(Temp)
        Mu <<- mu_lw(Temp)
        K <<- k_lw(Temp)
        
        # Reynolds #
        vs <<- input$Fs/(Dens[1:input$modelsteps]*CSA_s)
        vt <<- input$Ft/(Dens[(input$modelsteps+1):(2*input$modelsteps)]*CSA_t)
        Re_s <<- Dens[1:input$modelsteps]*dhs*vs/Mu[1:input$modelsteps]
        Re_t <<- Dens[(input$modelsteps+1):(2*input$modelsteps)]*dht*vt/Mu[(input$modelsteps+1):(2*input$modelsteps)]
        
        # Prandtl #
        Pr_s <<- Mu[1:input$modelsteps]*Cp[1:input$modelsteps]/K[1:input$modelsteps]
        Pr_t <<- Mu[(input$modelsteps+1):(2*input$modelsteps)]*Cp[(input$modelsteps+1):(2*input$modelsteps)]/K[(input$modelsteps+1):(2*input$modelsteps)]
        
        # Graetz #
        Gz_s <<- dhs/input$L*Re_s*Pr_s
        Gz_t <<- dht/input$L*Re_t*Pr_t
        
        # Nusselt #
        if (mean(Re_s) > 2100) {
          Nu_s <<- 0.027*Re_s^0.8*Pr_s^(1/3)*(Mu[1:input$modelsteps]/Mu[(2*input$modelsteps+1):(3*input$modelsteps)])^0.14
        } else {
          Nu_s <<- 1.86*Gz_s^(1/3)*(Mu[1:input$modelsteps]/Mu[(2*input$modelsteps+1):(3*input$modelsteps)])^0.14
        }
        
        if (mean(Re_t) > 2100) {
          Nu_t <<- 0.027*Re_t^0.8*Pr_t^(1/3)*(Mu[(input$modelsteps+1):(2*input$modelsteps)]/Mu[(2*input$modelsteps+1):(3*input$modelsteps)])^0.14
        } else {
          Nu_t <<- 1.86*Gz_t^(1/3)*(Mu[(input$modelsteps+1):(2*input$modelsteps)]/Mu[(2*input$modelsteps+1):(3*input$modelsteps)])^0.14
        }
        
        # Heat transfer coeffcient
        h_s <<- Nu_s*K[1:input$modelsteps]/dhs
        h_t <<- Nu_t*K[(input$modelsteps+1):(2*input$modelsteps)]/dht
        
        
        # Insert coefficents for shell
        T_matrix[1,1] = 1
        for (i in 2:input$modelsteps){
          T_matrix[i,i] = 1
          T_matrix[i,i-1] = -1 + h_s[i-1]*SA_ot*dz/(input$L*input$Fs*Cp[i-1])
          T_matrix[i,i-1+2*input$modelsteps] = -h_s[i-1]*SA_ot*dz/(input$L*input$Fs*Cp[i-1])
        }
        
        # Insert coefficents for tube
        T_matrix[2*input$modelsteps,2*input$modelsteps] = 1
        for (i in (input$modelsteps+1):(2*input$modelsteps-1)) {
          T_matrix[i,i] = 1
          T_matrix[i,i+1] = -1 + h_t[i+1-input$modelsteps]*SA_it*dz/(input$L*input$Ft*Cp[i+1])
          T_matrix[i,i+1+input$modelsteps] = -h_t[i+1-input$modelsteps]*SA_it*dz/(input$L*input$Ft*Cp[i+1])
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
        x_vector <<- lusys(T_matrix,rhs_vector)
        
        err = sqrt(mean(((x_vector-Temp)/Temp)^2))*100
        Temp = x_vector
        print(err)
      }
      
      # Initialize heatmap data matrix without baffle columns
      heatmapdata = matrix(ncol = input$modelsteps, nrow = 3)
      
      heatmapdata[1,1:input$modelsteps] = x_vector[1:input$modelsteps]
      heatmapdata[2,1:input$modelsteps] = x_vector[(input$modelsteps+1):(2*input$modelsteps)]
      heatmapdata[3,1:input$modelsteps] = x_vector[1:input$modelsteps]
      
      heatmapdata <<- heatmapdata # save as global variable for debugging
      
      # Create border lines for tubes in heatmap
      tubelines = list(
        list(type = 'line', xref = "x", yref = "y",
             x0 = -0.5, x1 = input$modelsteps-0.5,
             y0 = 1.5, y1 = 1.5,
             line = list(color = "black", width = 10)
        ),
        list(type = 'line', xref = "x", yref = "y",
             x0 = -0.5, x1 = input$modelsteps-0.5,
             y0 = 0.5, y1 = 0.5,
             line = list(color = "black", width = 10)
        )
      )
      
      # Create plot and output
      output$plot1 = renderPlotly({
        plot_ly(z = heatmapdata, type = "heatmap", hoverinfo = 'text', 
                text = list(paste(round(heatmapdata[1,1:(input$modelsteps)],digits = 2), "K", sep = " "),
                            paste(round(heatmapdata[2,1:(input$modelsteps)],digits = 2), "K", sep = " "),
                            paste(round(heatmapdata[3,1:(input$modelsteps)],digits = 2), "K", sep = " ")
                ),
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
    # Shell and tube with baffles --------------------------------------------------------------------------------------------  
    } else if (input$baffles > 24) {withProgress(value = 0.5, message = "Calculation in progress...",{
      
      # Number of sections seperated by baffles
      sections = input$baffles + 1
      
      # Number of model steps per section
      stepsPerSection = ceiling(input$modelsteps/sections)
      
      # Cross sectional area of the tube (m2)
      CSA_t <<- pi/4*(input$tid^2)
      
      # Hydraulic diameter pipe
      dht = input$tid
      
      # Calculate space between tubes
      pitch_x = input$sid/max(Ntubes_row + 1)
      pitch_y = input$sid/(length(Ntubes_row) + 1)
      
      # Chord length and cross-sectional area of shell per row
      chord_length = sqrt(2*pitch_y*(1:length(Ntubes_row))*input$sid/2 - (pitch_y*(1:length(Ntubes_row)))^2)
      CSA_s = c(chord_length*input$L/(sections*stepsPerSection),0)
      
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
      Temp <<- c(rep(input$Tsi,times = (nrows+1)*sections*stepsPerSection),
                 rep(input$Tti,times = nrows*sections*stepsPerSection),
                 rep(mean(c(input$Tti,input$Tsi)),times = nrows*sections*stepsPerSection))
      
      err = 1
      while (err > 0.001) { 
        # Create matrix containing all calclated temp values
        T_matrix = matrix(data = 0, ncol = stepsPerSection*sections*(3*nrows+1),nrow = stepsPerSection*sections*(3*nrows+1))
        
        # Initial values + averaging for next section
        for (i in seq(1,sections*stepsPerSection,by = stepsPerSection*2)) {
          
          T_matrix[i:(stepsPerSection+i-1),i:(stepsPerSection+i-1)] = diag(stepsPerSection)
          
        }
        for (i in seq(sections*stepsPerSection*nrows+stepsPerSection+1,sections*stepsPerSection*(nrows+1),by = stepsPerSection*2)) {
          
          T_matrix[i:(stepsPerSection+i-1),i:(stepsPerSection+i-1)] = diag(stepsPerSection)
          
        }
        for (i in seq(1+stepsPerSection*2,sections*stepsPerSection,by = stepsPerSection*2)) {
          
          T_matrix[i:(stepsPerSection+i-1),(i-stepsPerSection):(i-1)] = matrix(data = -1/stepsPerSection,ncol = stepsPerSection, nrow = stepsPerSection)
          
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
        
        # Properties
        Dens <<- dens_lw(Temp)
        Cp <<- cp_lw(Temp)
        Mu <<- mu_lw(Temp)
        K <<- k_lw(Temp)
         
        # Reynolds #
        vs = input$Fs/(Dens[1:(stepsPerSection*sections)]*CSA_s_vector)
        vt = input$Ft/(Dens[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*CSA_t*sum(Ntubes_row))
        
        vs_max = vs*pitch_x/(pitch_x-(input$tid + 2*input$tt))
        
        if (input$config == "triangle" && pitch_x >= input$tid + 2*input$tt + 2*(sqrt((pitch_x/2)^2 + pitch_y^2) - input$tid - 2*input$tt)) {
          pitch_diag = sqrt((pitch_x/2)^2 + pitch_y^2)
          vs_max = pitch_x/2*vs/(pitch_diag - input$tid - 2*input$tt)
        } 
        
        Re_s <<- Dens[1:(stepsPerSection*sections*(nrows+1))]*dht*vs_max/Mu[1:(stepsPerSection*sections*(nrows+1))]
        Re_t <<- Dens[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*dht*vt/Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]
        Re_d <<- Dens[1:(stepsPerSection*sections*(nrows+1))]*dht*vs/Mu[1:(stepsPerSection*sections*(nrows+1))]
        
        
        # Prandtl #
        Pr_s <<- Mu[1:(stepsPerSection*sections*(nrows+1))]*Cp[1:(stepsPerSection*sections*(nrows+1))]/K[1:(stepsPerSection*sections*(nrows+1))]
        Pr_t <<- Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/K[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]
        Pr_w <<- Mu[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))]*Cp[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))]/K[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))]
        
        
        # Create Pr_w vector that matches the shell vector
        Pr_w_aligned <<- vector(mode = "numeric",length = stepsPerSection*sections*(nrows+1)) - 999
        
        for (i in seq(1,stepsPerSection*sections*(nrows+1), by = 2*stepsPerSection)) {
          Pr_w_aligned[i:(i+stepsPerSection-1)] = Pr_w[i:(i+stepsPerSection-1)]
        }
        
        for (i in seq(stepsPerSection*sections + stepsPerSection + 1,stepsPerSection*sections*(nrows + 1), by = 2*stepsPerSection)) {
          Pr_w_aligned[i:(i+stepsPerSection-1)] = Pr_w[(i-stepsPerSection*sections):(i+stepsPerSection-stepsPerSection*sections-1)]
        }
        
        Pr_w_aligned[is.na(Pr_w_aligned)] = -999
        
        # Graetz #
        Gz_t <<- dht/input$L*Re_t*Pr_t
        
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
            if (mean(Re_s) < 100) {
              Nu_s = Nuss_func(0.8,0.4)
            } else if (mean(Re_s) < 1000) {
              Nu_s = 0.3 + (0.62*(Re_d^0.5)*Pr_s^(1/3))/((1+(0.4/Pr_s)^(2/3))^0.25)*((1+(Re_d/282000)^0.625)^(0.8))
            } else if (mean(Re_s) < 200000) {
              Nu_s = Nuss_func(0.27,0.63)
            } else {
              Nu_s = Nuss_func(0.021,0.84)
            }
            
          } else {
            if (mean(Re_s) < 100) {
              Nu_s = Nuss_func(0.9,0.4)
            } else if (mean(Re_s) < 1000) {
              Nu_s = 0.3 + (0.62*(Re_d^0.5)*Pr_s^(1/3))/((1+(0.4/Pr_s)^(2/3))^0.25)*((1+(Re_d/282000)^0.625)^(0.8))
            } else if (mean(Re_s) < 200000) {
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
          
        if (mean(Re_t) > 2100) {
          Nu_t <<- 0.027*Re_t^0.8*Pr_t^(1/3)*(Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/Mu[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))])^0.14
        } else {
          Nu_t <<- 1.86*Gz_t^(1/3)*(Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/Mu[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))])^0.14
        }
        
        
        # Heat transfer coeffcient
        h_s <<- Nu_s*K[1:(stepsPerSection*sections*(nrows+1))]/dht
        h_t <<- Nu_t*K[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/dht

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
        
        # Set up a new ht vector because stein said so (and to match hs)
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
        
        shellCoefficient1 = 1 - h_s*Tube_vector*SA_ot/(sections*2*input$Fs*Cp[1:(stepsPerSection*sections*(nrows+1))])
        shellCoefficient2 = -1 - h_s*Tube_vector*SA_ot/(sections*2*input$Fs*Cp[1:(stepsPerSection*sections*(nrows+1))])
        shellCoefficient3 = h_s*Tube_vector*SA_ot/(sections*input$Fs*Cp[1:(stepsPerSection*sections*(nrows+1))])
        tubeCoefficient1 = -1
        tubeCoefficient2 = 1 - h_t*sum(Ntubes_row)*SA_it/(sections*stepsPerSection*input$Ft*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))])
        tubeCoefficient3 = h_t*sum(Ntubes_row)*SA_it/(sections*stepsPerSection*input$Ft*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))])
        wallCoefficient1 <<- h_s*SA_ot/(2*(h_s*SA_ot + h_t_extended*SA_it))
        wallCoefficient2 = wallCoefficient1
        wallCoefficient3 <<- h_t_extended*SA_it/(h_s*SA_ot + h_t_extended*SA_it)
        wallCoefficient4 = -1
        
        # Shell for the shell equation (downward flow)
        for (i in seq(1+stepsPerSection*sections,sections*stepsPerSection*2,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = diag(stepsPerSection)*shellCoefficient2[(i+j-sections*stepsPerSection):(stepsPerSection+i+j-1-sections*stepsPerSection)]
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j-sections*stepsPerSection):(i+j-sections*stepsPerSection+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient1[(i+j-sections*stepsPerSection):(stepsPerSection+i+j-1-sections*stepsPerSection)]
          }
        }
        
        # Shell for the shell equation (upward flow)
        for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = diag(stepsPerSection)*shellCoefficient2[(i+j+sections*stepsPerSection):(stepsPerSection+i+j-1+sections*stepsPerSection)]
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j+sections*stepsPerSection):(i+j+sections*stepsPerSection+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient1[(i+j+sections*stepsPerSection):(stepsPerSection+i+j-1+sections*stepsPerSection)]
          }
        }
        
        # Wall part for the shell equation (downward flow)
        for (i in seq(1+stepsPerSection*sections,sections*stepsPerSection*2,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(sections*stepsPerSection*2*length(Ntubes_row)+i+j):(sections*stepsPerSection*2*length(Ntubes_row)+i+j+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient3[(i+j-sections*stepsPerSection):(stepsPerSection+i+j-1-sections*stepsPerSection)]
          }
        }
        
        # Wall part for the shell equation (upward flow)
        for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(sections*stepsPerSection*(2*length(Ntubes_row)+1)+i+j):(sections*stepsPerSection*(2*length(Ntubes_row)+1)+i+j+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient3[(i+j+sections*stepsPerSection):(stepsPerSection+i+j-1+sections*stepsPerSection)]  
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
        
        # Known inlet tube stuff for matrix
        for (i in seq(sections*stepsPerSection*(length(Ntubes_row)+2),sections*stepsPerSection*(2*length(Ntubes_row)+1),by = sections*stepsPerSection)) {
          T_matrix[i,i] = 1
          rhs_vector[i] = input$Tti
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
        T_matrix <<- T_matrix
        
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
        
        
        T_matrix <<- T_matrix
        rhs_vector <<- rhs_vector
        
        # Solve matrix problem using LU decomposition
        x_vector <<- lusys(T_matrix,rhs_vector)
        err = sqrt(mean(((x_vector - Temp)/Temp)^2))
        Temp = x_vector
        print(err)
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
        plot_ly(z = heatmapdata, type = "heatmap", hoverinfo = 'text', 
                text = heatmaptext,
                colorscale = "RdBu", 
                colorbar = list(len = 1, title = list(text = "Temperature (K)",side = "right"))) %>%
          config(displayModeBar = F) %>%
          layout(
            shapes = lines,
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
      
    # Shell and tube with baffles version 2 --------------------------------------------------------------------------------------------  
    } else if (input$baffles > 0) {withProgress(value = 0.5, message = "Calculation in progress...",{
      
      # Number of sections seperated by baffles
      sections = input$baffles + 1
      
      # Number of model steps per section
      stepsPerSection = ceiling(input$modelsteps/sections)
      
      # Cross sectional area of the tube (m2)
      CSA_t = pi/4*(input$tid^2)
      
      # Hydraulic diameter pipe
      dht = input$tid
      
      # Calculate space between tubes
      pitch_x = input$sid/max(Ntubes_row + 1)
      pitch_y = input$sid/(length(Ntubes_row) + 1)
      
      # Chord length and cross-sectional area of shell per row
      chord_length = sqrt(2*pitch_y*(1:length(Ntubes_row))*input$sid/2 - (pitch_y*(1:length(Ntubes_row)))^2)
      CSA_s = c(chord_length*input$L/(sections*stepsPerSection),0)
      
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
      
      # Create mass flowrate vector
      Fs = rep(input$Fs/stepsPerSection,stepsPerSection*sections*(nrows+1))
      
      # Create condensed mass vector
      Mcond = rep(0,stepsPerSection*sections*(nrows+1))
      
      # Create RH vector
      RH = rep(input$RH,stepsPerSection*sections*(nrows+1))
      
      err = 1
      karisux = 1
      while (karisux < 3) { 
        
        # Create matrix containing all calclated temp values
        T_matrix = matrix(data = 0, ncol = stepsPerSection*sections*(3*nrows+1),nrow = stepsPerSection*sections*(3*nrows+1))
        
        test1 <<- sections
        test2 <<-stepsPerSection
        test3 <<- T_matrix
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
        
        if (input$config == "triangle" && pitch_x >= input$tid + 2*input$tt + 2*(sqrt((pitch_x/2)^2 + pitch_y^2) - input$tid - 2*input$tt)) {
          pitch_diag = sqrt((pitch_x/2)^2 + pitch_y^2)
          vs_max = pitch_x/2*vs/(pitch_diag - input$tid - 2*input$tt)
        } 
        
        testvs <<- vs
        
        Re_s = Dens[1:(stepsPerSection*sections*(nrows+1))]*dht*vs_max/Mu[1:(stepsPerSection*sections*(nrows+1))]
        Re_t = Dens[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*dht*vt/Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]
        Re_d = Dens[1:(stepsPerSection*sections*(nrows+1))]*dht*vs/Mu[1:(stepsPerSection*sections*(nrows+1))]
        
        testres <<- Re_s
        testret <<- Re_t
        
        
        
        # Prandtl #
        Pr_s = Mu[1:(stepsPerSection*sections*(nrows+1))]*Cp[1:(stepsPerSection*sections*(nrows+1))]/K[1:(stepsPerSection*sections*(nrows+1))]
        Pr_t = Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/K[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]
        Pr_w = Mu[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))]*Cp[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))]/K[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))]
        
        
        # Create Pr_w vector that matches the shell vector
        Pr_w_aligned = vector(mode = "numeric",length = stepsPerSection*sections*(nrows+1)) - 999
        
        for (i in seq(1,stepsPerSection*sections*(nrows+1), by = 2*stepsPerSection)) {
          Pr_w_aligned[i:(i+stepsPerSection-1)] = Pr_w[i:(i+stepsPerSection-1)]
        }
        
        for (i in seq(stepsPerSection*sections + stepsPerSection + 1,stepsPerSection*sections*(nrows + 1), by = 2*stepsPerSection)) {
          Pr_w_aligned[i:(i+stepsPerSection-1)] = Pr_w[(i-stepsPerSection*sections):(i+stepsPerSection-stepsPerSection*sections-1)]
        }
        
        Pr_w_aligned[is.na(Pr_w_aligned)] = -999
        
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
            if (mean(Re_s) < 100) {
              Nu_s = Nuss_func(0.8,0.4)
            } else if (mean(Re_s) < 1000) {
              Nu_s = 0.3 + (0.62*(Re_d^0.5)*Pr_s^(1/3))/((1+(0.4/Pr_s)^(2/3))^0.25)*((1+(Re_d/282000)^0.625)^(0.8))
            } else if (mean(Re_s) < 200000) {
              Nu_s = Nuss_func(0.27,0.63)
            } else {
              Nu_s = Nuss_func(0.021,0.84)
            }
            
          } else {
            if (mean(Re_s) < 100) {
              Nu_s = Nuss_func(0.9,0.4)
            } else if (mean(Re_s) < 1000) {
              Nu_s = 0.3 + (0.62*(Re_d^0.5)*Pr_s^(1/3))/((1+(0.4/Pr_s)^(2/3))^0.25)*((1+(Re_d/282000)^0.625)^(0.8))
            } else if (mean(Re_s) < 200000) {
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
        
        if (mean(Re_t) > 2100) {
          Nu_t = 0.027*Re_t^0.8*Pr_t^(1/3)*(Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/Mu[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))])^0.14
        } else {
          Nu_t = 1.86*Gz_t^(1/3)*(Mu[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]/Mu[(stepsPerSection*sections*(2*nrows+1)+1):(stepsPerSection*sections*(3*nrows+1))])^0.14
        }
        
        testNUT <<- Nu_t
        testNUS <<- Nu_s
        testKKK <<- K
        
        # Heat transfer coeffcient
        h_s = Nu_s*K[1:(stepsPerSection*sections*(nrows+1))]/dht
        h_s = h_s*300
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
        testht <<- h_t_extended

        shellCoefficient1 = 1 - h_s*Tube_vector*SA_ot/(stepsPerSection*sections*2*Fs*Cp[1:(stepsPerSection*sections*(nrows+1))])
        shellCoefficient2 = -1 - h_s*Tube_vector*SA_ot/(stepsPerSection*sections*2*Fs*Cp[1:(stepsPerSection*sections*(nrows+1))])
        shellCoefficient3 = h_s*Tube_vector*SA_ot/(stepsPerSection*sections*Fs*Cp[1:(stepsPerSection*sections*(nrows+1))])
        # for (i in 1:(stepsPerSection*sections*(nrows+1))) {
        #   if (Temp[i] <= Tcond && input$fluid == "vw") {
        #     shellCoefficient1[i] =  -h_s[i]*Tube_vector[i]*SA_ot/(sections*2*stepsPerSection)
        #     shellCoefficient2[i] = shellCoefficient1[i]
        #     shellCoefficient3[i] = h_s[i]*Tube_vector[i]*SA_ot/(sections*stepsPerSection)
        #     
        #     rhs_vector[i] = -Mcond[i]*Hvap[i]
        #   }
        # }
        for (i in 1:(sections*stepsPerSection)) {
          for (j in 0:nrows){
            if (Temp[i+j*stepsPerSection*sections] <= Tcond && input$fluid == "vw") {
              shellCoefficient1[i+j*stepsPerSection*sections] =  -h_s[i+j*stepsPerSection*sections]*Tube_vector[i+j*stepsPerSection*sections]*SA_ot/(sections*2*stepsPerSection)
              shellCoefficient2[i+j*stepsPerSection*sections] = shellCoefficient1[i+j*stepsPerSection*sections]
              shellCoefficient3[i+j*stepsPerSection*sections] = h_s[i+j*stepsPerSection*sections]*Tube_vector[i+j*stepsPerSection*sections]*SA_ot/(sections*stepsPerSection)
              if (ceiling(i/stepsPerSection) %%2 == 0 && j != nrows) {
                rhs_vector[i+j*stepsPerSection*sections] = -Mcond[i+(j+1)*stepsPerSection*sections]*Hvap[i+(j+1)*stepsPerSection*sections]
              } else if (ceiling(i/stepsPerSection) %%2 == 1 && j != 0) {
                rhs_vector[i+j*stepsPerSection*sections] = -Mcond[i+(j-1)*stepsPerSection*sections]*Hvap[i+(j-1)*stepsPerSection*sections]
              }
            }
          }
        }
        
        
        tubeCoefficient1 = -1
        tubeCoefficient2 = 1 - h_t*sum(Ntubes_row)*SA_it/(sections*stepsPerSection*input$Ft*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))])
        tubeCoefficient3 = h_t*sum(Ntubes_row)*SA_it/(sections*stepsPerSection*input$Ft*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))])
        wallCoefficient1 = h_s*SA_ot/(2*(h_s*SA_ot + h_t_extended*SA_it))
        wallCoefficient2 = wallCoefficient1
        wallCoefficient3 = h_t_extended*SA_it/(h_s*SA_ot + h_t_extended*SA_it)
        wallCoefficient4 = -1
        
        # Shell for the shell equation (downward flow)
        for (i in seq(1+stepsPerSection*sections,sections*stepsPerSection*2,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = diag(stepsPerSection)*shellCoefficient2[(i+j-sections*stepsPerSection):(stepsPerSection+i+j-1-sections*stepsPerSection)]
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j-sections*stepsPerSection):(i+j-sections*stepsPerSection+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient1[(i+j-sections*stepsPerSection):(stepsPerSection+i+j-1-sections*stepsPerSection)]
          }
        }
        
        # Shell for the shell equation (upward flow)
        for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j):(stepsPerSection+i+j-1)] = diag(stepsPerSection)*shellCoefficient2[(i+j+sections*stepsPerSection):(stepsPerSection+i+j-1+sections*stepsPerSection)]
            T_matrix[(i+j):(stepsPerSection+i+j-1),(i+j+sections*stepsPerSection):(i+j+sections*stepsPerSection+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient1[(i+j+sections*stepsPerSection):(stepsPerSection+i+j-1+sections*stepsPerSection)]
          }
        }
        
        # Wall part for the shell equation (downward flow)
        for (i in seq(1+stepsPerSection*sections,sections*stepsPerSection*2,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(sections*stepsPerSection*2*length(Ntubes_row)+i+j):(sections*stepsPerSection*2*length(Ntubes_row)+i+j+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient3[(i+j-sections*stepsPerSection):(stepsPerSection+i+j-1-sections*stepsPerSection)]
          }
        }
        
        # Wall part for the shell equation (upward flow)
        for (i in seq(1+stepsPerSection,sections*stepsPerSection,by = stepsPerSection*2)) {
          for (j in seq(0,(length(Ntubes_row)-1)*sections*stepsPerSection,by = sections*stepsPerSection)) {
            
            T_matrix[(i+j):(stepsPerSection+i+j-1),(sections*stepsPerSection*(2*length(Ntubes_row)+1)+i+j):(sections*stepsPerSection*(2*length(Ntubes_row)+1)+i+j+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient3[(i+j+sections*stepsPerSection):(stepsPerSection+i+j-1+sections*stepsPerSection)]  
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
        
        
        T_matrix = T_matrix
        rhs_vector = rhs_vector
        
        #Fs = Fs
        testmatrix <<- T_matrix
        testb <<- rhs_vector
        
        
        # Solve matrix problem using LU decomposition
        x_vector = lusys(T_matrix,rhs_vector)
        testxvector <<- x_vector
        for (i in 1:(stepsPerSection*sections*(nrows+1))) {
          if (x_vector[i] < Tcond && input$fluid == "vw") {
            x_vector[i] = Tcond
          }
        }
        
        
        err = sqrt(mean(((x_vector - Temp)/Temp)^2))
        Temp = x_vector
        print(err)
        
        Q = input$Ft*Cp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))]*(
          Temp[(stepsPerSection*sections*(nrows+1)+1):(stepsPerSection*sections*(2*nrows+1))] - 
          Temp[(stepsPerSection*sections*(nrows+1)+2):(stepsPerSection*sections*(2*nrows+1)+1)])
        for (i in 1:nrows) {
          Q[stepsPerSection*sections*i] = Q[stepsPerSection*sections*i-1]
        }
        

        
        if (input$fluid == "vw"){
          # Set up a new ht vector and to match hs
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
          
          
          
          
          
          
        }
        
         testcond <<- Mcond
         #testhvap <<- Hvap
         testdens <<- Dens
         
         
        karisux = karisux +1
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
      
      
      heatmapdata = heatmapdata # save as global variable for debugging
      
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
        plot_ly(z = heatmapdata, type = "heatmap", hoverinfo = 'text', 
                text = heatmaptext,
                colorscale = "RdBu", 
                colorbar = list(len = 1, title = list(text = "Temperature (K)",side = "right"))) %>%
          config(displayModeBar = F) %>%
          layout(
            shapes = lines,
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
    }
  })
}





