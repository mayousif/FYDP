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
    
    if (input$baffles == 0 && input$ntubes == 1){
      
      parms = c(Ft = input$Ft, Fs = input$Fs, Cpt = input$Cpt, Cps = input$Cps, ht = input$ht, 
                hs = input$hs, SA_it = SA_it, SA_ot = SA_ot, dz = dz, L = input$L)
      
      doublePipeNoBaffles = function(z,Temp,parms){
        
        Tw = (parms[["hs"]]*Temp[1]*parms[["SA_ot"]] + parms[["ht"]]*Temp[2]*parms[["SA_it"]])/(parms[["ht"]]*parms[["SA_it"]] + parms[["hs"]]*parms[["SA_ot"]])
        dTs_dz= parms[["hs"]]*(Tw-Temp[1])*parms[["SA_ot"]]/(parms[["Fs"]]*parms[["Cps"]]*parms[["L"]])
        dTt_dz = -parms[["ht"]]*(Tw-Temp[2])*parms[["SA_it"]]/(parms[["Ft"]]*parms[["Cpt"]]*parms[["L"]])
  
        return(list(c(dTs_dz,dTt_dz)))
      }
      
      init = c(input$Tsi,NA)
      end = c(NA,input$Tti)
      
      sol = as.data.frame(bvptwp(yini = init, x = seq(0, input$L, by = dz), func = doublePipeNoBaffles, yend = end, parms = parms))
      
      output$plot1 = renderPlotly({
       plot_ly(x = sol[,1], y = sol[,2], type = "scatter", mode = "lines", name = "Shell-side") %>%
          add_trace(y = sol[,3], name = "Tube-side") %>%
          layout(
            xaxis = list(
              title = "x (m)",
              zeroline = F
            ),
            yaxis = list(
              title = "Temperature (K)",
              zeroline = F
            ),
            legend = list(
              orientation = "h",  
              xanchor = "center",
              x = 0.5,
              y = -0.2
            )
          )
      })
      
    } else if (input$baffles > 0 && input$ntubes == 1) {withProgress(value = 0.5, message = "Calculation in progress...",{
      
      # Number of sections seperated b baffles
      sections <<- input$baffles + 1
      
      # Number of model steps per section
      stepsPerSection <<- ceiling(input$modelsteps/sections)
      
      # Create matrix containing all calclated temp values
      T_matrix = matrix(data = 0, ncol = stepsPerSection*sections*4,nrow = stepsPerSection*sections*4)
      
      # Initial values + averaging for next section
      for (i in seq(1,sections*stepsPerSection*2,by = stepsPerSection*2)) {
        
        T_matrix[i:(stepsPerSection+i-1),i:(stepsPerSection+i-1)] = diag(stepsPerSection)
        
      }
      for (i in seq(1+stepsPerSection*2,sections*stepsPerSection*2,by = stepsPerSection*2)) {
        
        T_matrix[i:(stepsPerSection+i-1),(i-stepsPerSection):(i-1)] = matrix(data = -1/stepsPerSection,ncol = stepsPerSection, nrow = stepsPerSection)
        
      }
      
      # Right hand side solution to matrix equation
      rhs_vector = rep(0, stepsPerSection*sections*4)
      rhs_vector[1:stepsPerSection] = input$Tsi
      
      # Matrix coefficients
      shellCoefficient1 = 1 - input$hs*SA_ot/(sections*2*input$Fs*input$Cps)
      shellCoefficient2 = -1 - input$hs*SA_ot/(sections*2*input$Fs*input$Cps)
      shellCoefficient3 = input$hs*SA_ot/(sections*input$Fs*input$Cps)
      tubeCoefficient1 = -1
      tubeCoefficient2 = 1 - input$ht*SA_it/(sections*stepsPerSection*input$Ft*input$Cpt)
      tubeCoefficient3 = input$ht*SA_it/(sections*stepsPerSection*input$Ft*input$Cpt)
      wallCoefficient1 = input$hs*SA_ot/(2*(input$hs*SA_ot + input$ht*SA_it))
      wallCoefficient2 = wallCoefficient1
      wallCoefficient3 = input$ht*SA_it/(input$hs*SA_ot + input$ht*SA_it)
      wallCoefficient4 = -1
      
      # Shell top and bottom for the shell equation
      for (i in seq(1+stepsPerSection,sections*stepsPerSection*2,by = stepsPerSection*2)) {
        
        T_matrix[i:(stepsPerSection+i-1),i:(stepsPerSection+i-1)] = diag(stepsPerSection)*shellCoefficient2
        T_matrix[i:(stepsPerSection+i-1),(i-stepsPerSection):(i-1)] = diag(stepsPerSection)*shellCoefficient1
        
      }
      
      # Wall part for the shell equation
      for (i in seq(1+stepsPerSection,sections*stepsPerSection*2,by = stepsPerSection*2)) {
        j = stepsPerSection*sections*3+1 + (i-stepsPerSection-1)/2
        T_matrix[i:(stepsPerSection+i-1),j:(j+stepsPerSection-1)] = diag(stepsPerSection)*shellCoefficient3
      }
      
      # Tube section of matrix
      for (i in seq(1+stepsPerSection*sections*2,stepsPerSection*sections*3-1)) {
        T_matrix[i,i] = tubeCoefficient1
        T_matrix[i,i+1] = tubeCoefficient2
      }
      
      # Wall part of tube equation
      for (i in seq(2+stepsPerSection*sections*3,stepsPerSection*sections*4)) {
        T_matrix[i-stepsPerSection*sections-1,i] = tubeCoefficient3
      }
      
      # Known inlet tube stuff for matrix
      T_matrix[stepsPerSection*sections*3,stepsPerSection*sections*3] = 1
      rhs_vector[stepsPerSection*sections*3] = input$Tti
      
      #
      for (i in seq(1, stepsPerSection*sections,by = stepsPerSection)) {
        rows = (i+stepsPerSection*sections*3):(i+stepsPerSection*sections*3+stepsPerSection-1)
        col1 = (2*i-1):(2*i+stepsPerSection-2)
        col2 = (2*i+stepsPerSection-1):(2*i+stepsPerSection*2-2)
        col3 = rows - stepsPerSection*sections
        
        T_matrix[rows,col1] = diag(stepsPerSection)*wallCoefficient1
        T_matrix[rows,col2] = diag(stepsPerSection)*wallCoefficient2
        T_matrix[rows,col3] = diag(stepsPerSection)*wallCoefficient3
        T_matrix[rows,rows] = diag(stepsPerSection)*wallCoefficient4
        
      }
      
      # Solve matrix problem using LU decomposition
      x_vector <<- lusys(T_matrix,rhs_vector)
      
    })
      
      # Initialize heatmap data matrix without baffle columns
      heatmapdata = matrix(ncol = stepsPerSection*sections, nrow = 3)
      
      # Insert tube temperatures
      heatmapdata[2,] = x_vector[(stepsPerSection*sections*2 + 1):(stepsPerSection*sections*3)]
      
      # Insert shell temperatures
      hmsection = 1
      counter = 2
      for (i in 1:(sections*2)) {
        
        block = ceiling(i/2)
        if (counter == 2 && hmsection == 1) {
          heatmapdata[hmsection,(1 + (block-1)*stepsPerSection):(block*stepsPerSection)] = x_vector[(1 + stepsPerSection*(i-1)):(stepsPerSection*i)]
          hmsection = 3
          counter = 1
        } else if (counter == 2 && hmsection == 3) {
          heatmapdata[hmsection,(1 + (block-1)*stepsPerSection):(block*stepsPerSection)] = x_vector[(1 + stepsPerSection*(i-1)):(stepsPerSection*i)]
          hmsection = 1
          counter = 1
        } else if (counter < 2 && hmsection == 1) {
          heatmapdata[hmsection,(1 + (block-1)*stepsPerSection):(block*stepsPerSection)] = x_vector[(1 + stepsPerSection*(i-1)):(stepsPerSection*i)]
          counter = counter + 1
        } else if (counter < 2 && hmsection == 3) {
          heatmapdata[hmsection,(1 + (block-1)*stepsPerSection):(block*stepsPerSection)] = x_vector[(1 + stepsPerSection*(i-1)):(stepsPerSection*i)]
          counter = counter + 1
        }
        
      }
      
      heatmapdata <<- heatmapdata # save as global variable for debugging
      
      # Create border lines for tubes in heatmap
      tubelines = list(
        list(type = 'line', xref = "x", yref = "y",
             x0 = -0.5, x1 = sections*stepsPerSection-0.5,
             y0 = 1.5, y1 = 1.5,
             line = list(color = "black", width = 10)
        ),
        list(type = 'line', xref = "x", yref = "y",
             x0 = -0.5, x1 = sections*stepsPerSection-0.5,
             y0 = 0.5, y1 = 0.5,
             line = list(color = "black", width = 10)
        )
      )
      
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
          baffleline[c("y0", "y1")] = c(0,2.5)
        } else {
          baffleline[c("y0", "y1")] = c(-0.5,2)
        }
        lines = c(lines, list(baffleline))
      }
      
      lines = c(lines,tubelines)
      
      # Create plot and output
      output$plot1 = renderPlotly({
        plot_ly(z = heatmapdata, type = "heatmap", hoverinfo = 'text', 
                text = list(paste(round(heatmapdata[1,1:(sections*stepsPerSection)],digits = 2), "K", sep = " "),
                            paste(round(heatmapdata[2,1:(sections*stepsPerSection)],digits = 2), "K", sep = " "),
                            paste(round(heatmapdata[3,1:(sections*stepsPerSection)],digits = 2), "K", sep = " ")
                ),
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
      
    } else if (input$baffles > 0 && input$ntubes > 1) {withProgress(value = 0.5, message = "Calculation in progress...",{
      
      
      })
    }
    
  })
}





