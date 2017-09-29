wt <-
function(x, start = 1, dt = 1, dj = 1/20, 
         lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3),
         make.pval = T, method = "white.noise", params = NULL, 
         n.sim = 100, save.sim = F, outfile) {
 
                                 
  ###############################################################################
  ## Call function WaveletTransform
  ## Retrieve the wavelet transform, power, phases, amplitudes
  ###############################################################################
  
  # wavelet transform
  WT = WaveletTransform(x, dt = dt, dj = dj, 
                        lowerPeriod = lowerPeriod, upperPeriod = upperPeriod)

  print("Saving main run results...")
  save(list = c("WT"), file = outfile, compress = T)
  
  # keep these for COI and p-value calculations
  Power = WT$Power
  Power.avg = WT$Power.avg
  
  Period = WT$Period
  nr  = WT$nr
  nc  = WT$nc
  
  rm(WT)
  gc()
  
  ##################################################################################################
  ## Compute the power ridge
  ##################################################################################################
  
  Ridge = ridge(Power)
  gc()
  
  ###############################################################################
  ## Compute the cone of influence COI
  ###############################################################################
  
  coi = COI(start = start, dt = dt, nc = nc, nr = nr, Period = Period)
  coi = list(Ridge, coi.1=coi$coi.1, coi.2=coi$coi.2,
  		   axis.1=coi$axis.1, axis.2=coi$axis.2, dt=dt, dj=dj)
  
  print("Saving auxiliary info on the main run...")
  save(list = c("coi"), file = paste0(outfile, "_aux"), compress=T)
  rm(Period, Ridge, coi)
  gc()
  

  ###############################################################################
  ## Compute p values for significance check
  ###############################################################################
    
  Power.pval = NULL
  Power.avg.pval = NULL
  series.sim = NULL
  
  if (make.pval == T) {
       
      Power.pval = matrix(0, nrow = nr, ncol = nc)
      Power.avg.pval = rep(0, nr)
      
      if (save.sim == T) { series.sim = matrix(NA, nrow=nc, ncol=n.sim) }
       
      pb = txtProgressBar(min = 0, max = n.sim, style = 3) # create a progress bar
      for(ind.sim in 1:n.sim){
      
          x.sim = SurrogateData(x, method = method)
          
          if (save.sim == T) { series.sim[,ind.sim] = x.sim }
          
          WT.sim = WaveletTransform(x.sim, dt = dt, dj = dj, 
                                    lowerPeriod = lowerPeriod, upperPeriod = upperPeriod)
                                   
          Power.sim = WT.sim$Power
          Power.avg.sim = rowMeans(Power.sim)  
          
          rm(WT.sim)
          gc()
          
          Power.pval[Power.sim >= Power] = Power.pval[Power.sim >= Power] + 1
          Power.avg.pval[Power.avg.sim >= Power.avg] = Power.avg.pval[Power.avg.sim >= Power.avg] + 1
          setTxtProgressBar(pb, ind.sim)
      }
      close(pb)

      # p-values
      
      Power.pval = Power.pval / n.sim 
      Power.avg.pval = Power.avg.pval / n.sim  
      
      print("Saving simulation results...")
      save(list = c("Power.pval"), file = paste0(outfile, "_pval"), compress = T)
      save(list = c("Power.avg.pval"), file = paste0(outfile, "_pvalavg"), compress = T)
  }  
  

  ###############################################################################
  ## Prepare the output
  ###############################################################################

  output = list(Power = Power, Power.avg = Power.avg,
                Power.pval = Power.pval, Power.avg.pval = Power.avg.pval, 
                nc = nc, nr = nr,    
                series.sim = series.sim)
                
  return(invisible(output))
}
