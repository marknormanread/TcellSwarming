# Contains functionality common to both time-variant and -invariant forms of the swarming metric. 
#
# Jack Hywood and Mark N. Read, 2018
# 
# 

# Functionality common to both time-variant and -invariant forms of the metric. 
swarm_metric_preprocessing <- function(dtb_mat, time_step_indices, 
                                       rB, rR, r_resolution=1,
                                       time_step_sec,
                                       nsplines=NULL,
                                       lambda=NULL,
                                       nsplines_as_prop_rR=0.1,
                                       select_nsplines_prop=TRUE,
                                       min_loglam=-3, max_loglam=10,
                                       working_dir,  # Where to write graphs to.
                                       Norder=6, Lfdval=4,  # Used for smoothing with B splines. 
                                       # Limits for plotting. Expressed relative to the bolus=0
                                       r_plot_min=-rB, r_plot_max=(rR-rB), 
                                       kymograph_zlim=NULL  # re: surface plots. Set to c(lower, upper). 
                                       )
{
  rA <- rR-rB  # Defines the radius of the annulus from the bolus to the edge of the observed area.
  time_step_min = time_step_sec / 60
  time_step_hr = time_step_min / 60
  result <- initialise_Kraw(time_step_indices, rA, rB, rR, r_resolution=r_resolution, dtb_mat, r_plot_max=r_plot_max)
  
  r_full          <- result$r_full
  Kraw_mat        <- result$Kraw_mat
  n_agents_total  <- result$n_agents_total
  n_agents_R      <- result$n_agents_R
  n_agents_rmax   <- result$n_agents_rmax
  n_agents_B      <- result$n_agents_B
  n_agents_B_adj  <- result$n_agents_B_adj
  R_vol           <- result$R_vol
  
  n_time_steps = length(time_step_indices)
  # Write csv file for K functions:
  write.table(Kraw_mat, file=paste(working_dir, "/Kraw.txt", sep=''))  # Keeps the rownames
  
  print("Calculating M and L")
  result <- well_L_M(r_full=r_full, Kraw_mat=Kraw_mat, rB=rB, rR=rR)
  L_mat <- result$L_mat
  M     <- result$M
  r_M_indices <- result$r_M_indices  # Corresponds to indices of r >= 0. 

  r_plot_range_indices = which(r_plot_min <= r_full & r_full <= r_plot_max)  # Includes r within the bolus. 
  r_plot_values = r_full[r_plot_range_indices] 
  
  # Save M to filesystem.
  M_df <- data.frame(time_step_indices, M)
  names(M_df)[1] <- "Time_step"
  M_df$Time_s = (M_df$Time_step - 1) * time_step_sec  # -1 because R starts counting from 1, not zero. 
  write.csv(M_df, file=paste(working_dir, "/M", ".csv", sep=''), row.names=FALSE)

  pdf(paste(working_dir, '/M.pdf', sep=''))
  plot(time_step_indices * time_step_hr, M, type="l", lwd=2, ylim=c(-1, 1), xlab="Time (h)", ylab="M")
  dev.off()
  
  fds_L_mat <- fds(x=r_plot_values[r_M_indices], y=t(L_mat))  # (raw) L function values
  pdf(paste(working_dir, '/L.pdf', sep=''))
  plot(fds_L_mat, xlab="r", ylab="L", ylim=c(0, max(L_mat)))
  dev.off()
  
  # -----------------------------------
  # Set up for smoothing time series data with b-splines.
  Norder=6
  Lfdval=4

  # The following sets the number of splines if they have not been set by the user.
  #
  # From our experience analysing simulations and experimental data, the most appropriate
  # number of splines to use is dependent on the data being analysed.
  #
  # In our experience, we find that setting the number of splines ~0.05 x rA works well
  # for analysis, and that setting the number of splines <20 is problematic.
  #
  # However, we suggest that users investigate their data for a number of choices of splines.
  if (is.null(nsplines)) {
    print("User has not selected number of splines to use in smoothing (nsplines).")
    if (select_nsplines_prop) {
      print("Choosing nsplines as a proportion of r.")
      nsplines <- min(200, max(20, ceiling(rR * nsplines_as_prop_rR)))
    } else {
      # The following chooses the number of splines used based on minimising the GCV measure.
      # While this approach is commonly employed in this manner, we find that it typically significantly overestimates the
      # appropriate number of splines required.
      # We believe that this is likely due to a) commonly analysing over a large range of r values, and b) the fact that K
      # functions are monotonically increasing step functions.
      print("Choosing nsplines based on best fit as assessed through generalised cross validation.")
      nsplines <- choose_nsplines_gcv(rA=rA, r_full=r_full, Kraw_mat=Kraw_mat, Lfdval=Lfdval, Norder=Norder)
    }
  }
  print(paste('Number of splines to be used in smoothing =', nsplines))

  # Generate smoothed K functions, and choose the spline basis for other smoothing operations.
  result <- smooth_Kraw(r_full=r_full, Kraw_mat=Kraw_mat, 
                        nsplines=nsplines, lambda=lambda,
                        min_loglam=min_loglam, max_loglam=max_loglam, 
                        Lfdval=Lfdval, Norder=Norder)
  Ksmooth           <- result$Ksmooth
  gamma_splinebasis <- result$gamma_splinebasis
  
  
  # First and second derivatives of K (smoothed) w.r.t. r.
  Ksm_dr <- deriv.fd(Ksmooth$fd, 1)
  Ksm_dr_mat <- eval.fd(r_full, Ksm_dr) 
  Ksm_dr2 <- deriv.fd(Ksmooth$fd, 2)
  Ksm_dr2_mat <- eval.fd(r_full, Ksm_dr2)
  
  
  # Find differences in smoothed K between r and r+1, for all time points (and r). Has shape [r-1, n_time_steps-1]
  delta_Ksmooth <- (eval.fd(r_full, Ksmooth$fd)[,-1]  # From all Ksmooth t=2:end ...
                    - eval.fd(r_full, Ksmooth$fd)[, -n_time_steps]  # ... subract all Ksmooth t=1:end-1. IE, find dKsmooth.
                    ) / (time_step_min)
  deltaK_fd <- Data2fd(argvals=r_full, y=delta_Ksmooth, basisobj=gamma_splinebasis)
  
  # Repeats n_agents_R x [number of time steps]
  n_agents_R_mat <- t(matrix(rep(n_agents_R, length(r_full)), nrow=n_time_steps))

  # lambda (densities) calculated from the r=0=origin. 
  lambda_r <- eval.fd(r_full, Ksm_dr) * n_agents_R_mat / ( (rB+r_full) * 2 * pi * R_vol )
  lambda_dr <- ( (eval.fd(r_full, Ksm_dr2) * (rB + r_full)) - eval.fd(r_full, Ksm_dr) )  *
               ( n_agents_R_mat / ( (rB+r_full)^2 * 2 * pi * R_vol ) )  # First derivative

  # "fds" variables arrange data for rainbow plots. 
  fds_Kraw <- fds(x=r_full, y=t(Kraw_mat))  # Raw K function values
  pdf(paste(working_dir, '/K_raw.pdf', sep=''))
  plot(fds_Kraw, xlab="r", ylab="K (raw)", ylim=c(0, max(Kraw_mat)), xlim=range(r_plot_values))
  dev.off()
  
  pdf(paste(working_dir, '/nTotalhr.pdf', sep=''))
  plot(time_step_indices * time_step_min / 60, n_agents_total, ylim=c(0, max(n_agents_total)), col="black", type="l", xlab="Time", ylab="Agents")
  lines(time_step_indices * time_step_min / 60, n_agents_R, col="red")
  lines(time_step_indices * time_step_min / 60, n_agents_rmax, col="blue")
  legend("bottomleft", legend=c("total observed", "<= rA", "<= r_plot_max"), col=c("black", "red", "blue"), lty=1, cex=0.8)
  dev.off()
  
  pdf(paste(working_dir, '/agents_in_bolus.pdf', sep=''))
  plot(time_step_indices * time_step_min / 60 , n_agents_B, ylim=c(0, max(n_agents_B)), type="l", xlab="Time", ylab="Agents")
  dev.off()
  
  fds_Ksm <- fds(x=r_full, y=eval.fd(r_full, Ksmooth$fd))  # Smoothed K function values
  pdf(paste(working_dir, '/K_smoothed.pdf', sep=''))
  plot(fds_Ksm, xlab="r", ylab="K (smoothed)", ylim=c(min(Kraw_mat), max(Kraw_mat)), xlim=range(r_plot_values))
  dev.off()
  
  fds_Kdr <- fds(x=r_full, y=eval.fd(r_full, Ksm_dr))
  pdf(paste(working_dir, '/K_dr.pdf', sep=''))
  plot(fds_Kdr, xlab="r", ylab="K'", xlim=range(r_plot_values))
  dev.off()
  
  fds_lambda_r <- fds(x=r_full, y=lambda_r)
  pdf(paste(working_dir, '/lambda_r.pdf', sep=''))
  plot(fds_lambda_r, xlab="r", ylab="lambda", xlim=range(r_plot_values))
  dev.off()
  
  fds_lambda_dr <- fds(x=r_full, y=lambda_dr)
  pdf(paste(working_dir, '/lambda_dr.pdf', sep=''))
  plot(fds_lambda_dr, xlab="r", ylab="lambda'", 
       xlim=range(r_plot_values),
       ylim=range(lambda_dr[floor(dim(lambda_dr)[1]/5),], lambda_dr[ceiling(dim(lambda_dr)[1]*4/5),]))
  dev.off()
  

  # Plot data as surfaces to separate out time and r.
  # source("diverge_anchored.R")
  plot_surface(x=r_full, y=time_step_indices * time_step_hr , z=Ksm_dr_mat, centred_val=0,
               path=paste(working_dir, '/K_dr_surface.pdf', sep=''))

  # Converts density units from cells/um2 to cells/mm2.
  lambda_r_plot = lambda_r * 1e6
  # Specify a zlim common to multiple experiments
  # Lambda can exhibit some very severe large-magnitude values at the extremes of r. 
  # This is one reason it can be worth subsetting the data across R. 
  # Subsetting performed on the input to these surface plots such that zval ranges aren't crazy (the plot will calculate
  # the colour bar ranges for everything input, not just the subsets requested with xlim and ylim).
  plot_surface(x=r_full[r_plot_range_indices[-1]], y=time_step_indices * time_step_hr, z=lambda_r_plot[r_plot_range_indices[-1], ],
               path=paste(working_dir, '/lambda_r_surface_commonZlim.pdf', sep=''),
               xlim=range(r_plot_values),
               zlim=kymograph_zlim, 
               plot_contours=FALSE)
  # Unrestricted zlim
  plot_surface(x=r_full[r_plot_range_indices[-1]], y=time_step_indices * time_step_hr, z=lambda_r_plot[r_plot_range_indices[-1], ],
               xlim=range(r_plot_values),
               path=paste(working_dir, '/lambda_r_surface.pdf', sep=''),
               plot_contours=FALSE)
  print(paste('range of lambda values = ', range(lambda_r_plot), sep=''))

  result <- list(r_full=r_full,
                 Ksmooth=Ksmooth,
                 gamma_splinebasis=gamma_splinebasis,
                 delta_Ksmooth=delta_Ksmooth,
                 deltaK_fd=deltaK_fd,
                 Kraw_mat=Kraw_mat,
                 n_agents_R=n_agents_R,
                 n_agents_B=n_agents_B,
                 R_vol=R_vol,
                 n_agents_R_mat=n_agents_R_mat,
                 lambda_r=lambda_r,
                 lambda_dr=lambda_dr,
                 Norder=Norder,
                 Lfdval=Lfdval,
                 nsplines=nsplines,
                 lambda_r_plot=lambda_r_plot,
                 r_plot_range_indices=r_plot_range_indices
                )
  print("Completed common pre-processing stage.")
  return(result)
}


plot_surface <- function(x, y, z, xlab='r', ylab='Time', centred_val=NULL, path, 
                         xlim=range(x), ylim=range(y), zlim=NULL,  # Set to c(lower, upper) to specify limits. 
                         plot_contours=TRUE)
{
  if (is.null(zlim)) {
    zlim = range(z)
  }
  min_val <- min(z, na.rm=TRUE)
  max_val <- max(z, na.rm=TRUE)
  colmap <- NULL
  if (! is.null(centred_val)) {
    # MNR commented out @ 20181126; incompatible with zlim use to specify a common max range of colours plotted
    # across several plots. designer.colors works, but be warned, the help for that function is utterly terrible. 
    # colmap <- diverge_anchored(min_val=min_val, max_val=max_val, whiteCentred=centred_val, num_colours=256)
    colmap = designer.colors(250, c("blue","white", "red"))
  }
  pdf(path)
  # Automatically includes a colour bar
  image.plot(x=x, y=y, z=z, xlab=xlab,  ylab=ylab, col=colmap,
             xlim=xlim, ylim=ylim, zlim=zlim)
  
  if (plot_contours) {
    contour(x=x, y=y, z=z, add=TRUE, nlevels=1,
            levels=pretty(range(min_val, max_val), 8), col='dark grey')
  
    # Heavier weight for z=centred_val, if selected.
    if (! is.null(centred_val)) {
      contour(x=x, y=y, z=z, add=TRUE, nlevels=1, levels=c(centred_val), lwd=1.0)
    }
  }
  dev.off()
}


# Function inputs:
# time_step_indices - numeric vector containing the absolute-time stamps for the discrete time points; no repetition of values.
# rA - numeric, the radius of annulus around (and not including) the bolus. 
# r_resolution - Sets resolution of r over which analysis is performend.
# dtb_mat - list, contains items "distance_bolus" and "Time". Distance_bolus values <= 0 are considered within the bolus.
#
# Return in list:
# r_full - vector of "r" indices to process, which can be a periodic subset of all available. 
#     The indices refer to values stored in dtb_mat. 
# Kraw_mat - 2D matrix of numeric values, [time_step, r index]. r_resolution maps "r index" values onto the range 0:rA.
# n_agents_B - numeric vector[time_step], number of agents within bolus at given time step
# n_agents_R - numeric vector[time_step], number of agents within R (including bolus) at given time step. 
initialise_Kraw <- function(time_step_indices, rA, rB, rR, r_resolution=1, dtb_mat, 
                            r_plot_max)  # Not used in calculating K, only the number of agents captured at this "r".
{
  n_time_steps <- length(time_step_indices)  # Number of time steps in data set
  
  # r values used to produce K functions; i.e. rB + r. Length out can be adjusted as required.
  # Can contain negative values for regions WITHIN the bolus. 
  # 'r' here must match how the input data is spatially defined, which is 0=bolus edge.  
  r_full <- seq(ceiling(-rB), floor(rA), length.out=floor((rR) / r_resolution))  # Stick to whole numbers.
  print(paste('range of r values to be analysed (r_full) = ', min(r_full), ' to ', max(r_full), sep=''))
  # When simulaing 3D, make sure you increase chemokine secretion rates in proportion to 2D vs 3D volume. 
  R_vol <- pi * rR^2  # Volume of region R; 2D.
  
  # Following used to calculate Well K functions, and number of agents in A and B for each time point.
  # This is raw data, it is smoothed later in the script.
  Kraw_mat <- matrix(0, n_time_steps, length(r_full))
  n_agents_B <- rep(0, n_time_steps)  # Number of agents in B, i.e. inside bolus, at time step j
  n_agents_A <- rep(0, n_time_steps)  
  # Number of agents within region R at time step j. R may be selected as subset of total environment 
  # (e.g. if bolus is non-centric in well).
  n_agents_R <- rep(0, n_time_steps)  
  n_agents_rmax <- rep(0, n_time_steps)  # Num of agents within the user-selected maximum range of r to examine
  n_agents_total <- rep(0, n_time_steps)  # Number of agents within entire simulation/experiment.
  
  # Populate Kraw_mat, iterate over time steps. 
  for (j in 1:n_time_steps) 
  {
    time_sample <- time_step_indices[j]  # Extract the absolute time stamp of this time step.
    # Collect cell distances from bolus, for cells observed at time j
    distances_j <- dtb_mat$distance_bolus[which(dtb_mat$Time == time_sample)]
    n_agents_B[j] <- sum(distances_j <= 0)  # Number of agents within bolus. Value of <=0 indicates inside bolus.
    # Number of agents within R, which may be a subset of actual environment
    n_agents_R[j] <- sum(distances_j <= rA)  # rA not defined from origin, is intead defined from bolus edge (r=0)
    n_agents_A[j] <- n_agents_R[j] - n_agents_B[j]
    n_agents_rmax[j] <- sum(distances_j <= r_plot_max)  # r_plot_max may != rA, if subset of r subject to analysis. 
    n_agents_total[j] <- length(distances_j)  # Number of agents in the total environment (not subsetted)

    # Produce a vector counting the number of agents within "r" of the bolus (rB + r of centre) at current time step.
    # r_full can contain negative values if taking consideration of regions WITHIN the bolus. Bolus edge = 0. 
    k <- sapply(r_full, function(r) sum(distances_j <= r))  # Well K function for this point in time. 

    k <- k * R_vol / n_agents_R[j]

    Kraw_mat[j, ] <- k
    print(paste('Building Well K functions, processing time index ', j, sep=''))
  }

  result <- list(r_full=r_full, Kraw_mat=Kraw_mat, 
                 n_agents_B=n_agents_B,
                 n_agents_R=n_agents_R, 
                 n_agents_total=n_agents_total,
                 n_agents_rmax=n_agents_rmax,
                 R_vol=R_vol)
  return(result)
}


# Create L functions and calculate M. 
# Inputs:
# r_full - vector of "r" values to analyse. Takes values of -rB (the origin) to rA: corresponds to r=0=bolus edge. 
# Kraw_mat - 2D matrix of numeric values, [time_step, r index]. r_resolution maps "r index" values onto the range 0:rA.
# rB - radius of bolus
# rR - radius of entire region R (including bolus)   
well_L_M <- function(r_full, Kraw_mat, rB, rR)
{
  # Produce L_{i}(r) functions. 2D matrix, L_mat = [time step, r values]
  R_vol <- pi * rR^2  # Volume of region R; 2D.
  K = Kraw_mat
  
  if (! (min(r_full) <= 0 & max(r_full) > 0)) { stop("r_full does not contain 0.") }

  r_M_indices <- which(0 <= r_full)
  L_mat <- sqrt((K[, r_M_indices] / R_vol*(pi*(rR^2-rB^2)) + pi * rB^2) /pi) - rB
  
  # Produce metric M_{i}.
  dr <- r_full[2] - r_full[1]  # Step distnace in r, for integration below.
  
  # Complicated operation. We've tried to separate it across lines for clarity.
  M <- (2 / (rR - rB)^2) *  # Normalisation coefficient.
    (rowSums(t(  # M value calculated for each point in time, hence summing over r.
      t(L_mat) - dr*(0:(length(r_M_indices)-1)))  # Subtract r values from L values, for given r. Processes all times at once.
    ) * dr  # Integration.
    )
  
  return(list(L_mat=L_mat,  # L_mat - 2D numeric matrix [time step, r value]. 
              M=M,  # M - numeric vector[time step], contains values of swarming metric M. 
              r_M_indices=r_M_indices))  # indices of r_full used in defining L and M (r >= 0)
}


# The following can be used IF you want to choose the number of splines used based on minimising the GCV measure.
# While this approach is commonly employed in this manner, we find that it typically significantly overestimates the 
# appropriate number of splines required.
# We believe that this is likely due to a) commonly analysing over a large range of r values, and b) the fact that K 
# functions are monotonically increasing step functions.
choose_nsplines_gcv <- function(rA, r_full, Kraw_mat, Lfdval, Norder=6)
{
  nsplines_vec <- seq(10, 10 * floor(rA * 0.25/10), 10)  # Possible numbers of splines to trail. 
  n_vec <- length(nsplines_vec)
  gcvsave <- rep(NA, n_vec)
  for (ispl in 1:n_vec) {  # Cycle through each (log) lambda value
    nsplines_i <- nsplines_vec[ispl]
    gamma_splinebasis <- create.bspline.basis(c(0, rA), nbasis=nsplines_i, norder=Norder)
    fdParobj <- fdPar(gamma_splinebasis, Lfdobj=Lfdval)
    smoothlist <- smooth.basis(argvals=r_full, y=t(Kraw_mat), fdParobj)
    gcvsave[ispl] <- sum(smoothlist$gcv)
  }
  print(paste('nsplines selected through gcv procedure =', nsplines_vec[which.min(gcvsave)]))
  nsplines <- nsplines_vec[which.min(gcvsave)]  # The optimal lambda value
  return(nsplines)
}


smooth_Kraw <- function(r_full, Kraw_mat, nsplines, lambda, min_loglam, max_loglam, Lfdval, Norder)
{
  # Choose the the optimal smoothing parameter lambda by minimising GCV.
  # We find that this performs well, resulting in accurate estimates by the linear model.
  gamma_splinebasis <- create.bspline.basis(range(r_full), nsplines, norder=Norder)
  if (is.null(lambda)) {
    print("User has not selected the smoothing penalty to use in smoothing (lambda); determine using cross validation.")
    # Gives ranges of log(lambda) values searched. 
    # May need to be adjusted for individual data sets, though this range has worked well for all of ours.
    # Warning messages are generated if the max/min value specified by user is selected; the optimal value may lie beyond 
    # the range explored. 
    loglam <- seq(min_loglam, max_loglam, 0.25)
    nlam <- length(loglam)
    gcvsave <- rep(NA, nlam)
    for (ilam in 1:nlam) {  # cycle through each (log) lambda value
      lambda_i <- 10^loglam[ilam]
      fdParobj <- fdPar(gamma_splinebasis, Lfdobj=Lfdval, lambda=lambda_i)
      smoothlist <- smooth.basis(argvals=r_full, y=t(Kraw_mat), fdParobj)
      gcvsave[ilam] <- sum(smoothlist$gcv)
    }
    optimal_loglam <- loglam[which.min(gcvsave)]
    print(paste('Optimal log10 lambda =', optimal_loglam))
    lambda <- 10^(optimal_loglam)  # The optimal lambda value
    if (min_loglam == optimal_loglam) { 
      print(paste("WARNING! GCV loglam procedure selected the minimum specified possibility, 10^", min_loglam,
                  "\n  You may wish to select a LOWER min_loglam value to ensure there isn't a better possibility", sep=''))
    }
    if (max_loglam == optimal_loglam) { 
      print(paste("WARNING! GCV loglam procedure selected the maximum specified possibility, 10^", min_loglam,
                  "\n  You may wish to select a HIGHER min_loglam value to ensure there isn't a better possibility", sep=''))
    }
  }
  print(paste('Smoothing parameter =', lambda))
  # Using optimal lambda as established above.
  fdParobj <- fdPar(gamma_splinebasis, Lfdobj=Lfdval, lambda=lambda)
  Ksmooth <- smooth.basis(argvals=r_full, y=t(Kraw_mat), fdParobj)  # Smoothed K function values
  
  return(list(Ksmooth=Ksmooth, gamma_splinebasis=gamma_splinebasis))
}
