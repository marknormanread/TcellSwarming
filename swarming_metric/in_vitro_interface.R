#!/usr/bin/env Rscript
#
# Can be run from within R (I use Rstudio) or from the command line:
# $> Rstudio script_name.R <dataset name>
#
# Script is used to make publication quality-level graphs. 
#
# Scipt decomposes spatio temporal positional data of agents around a region (bolus) in chemoattractive (or repulsive)
# and diffusive components.
#
# The input csv file must have a column named "Time", and another named "distance_bolus".
#
# Can be launched from command line with with Rscript, but only from this working directory (to load dependent scripts 
# contained herein).
#
# lsos.R can be a helpful debuging tool.
# $> source("lsos.R")
# $> lsos()
#
# Mark N. Read, 2017 
rm(list=ls())

require(fda)
require(rainbow)
require(fields)

user = "Mark"
if (user == "Mark") {
  setwd("~/dropbox_usyd/projects/biro/neutroswarm/imaris/swarming_metric_biro/")  
}


# Representative experiments. Comment out as needed. 
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  #-----
  # dataset = "Cog-preembed31817"  # See if statements below for text to put here, to run a specific dataset.
  # dataset = "Cog6917"  # See if statements below for text to put here, to run a specific dataset.
  # dataset = "Non-cognate8917"  # See if statements below for text to put here, to run a specific dataset.
  #-----
  # dataset = "human_control"
  # dataset = "human_aCD3CD28"
  #-----
  # dataset = "EL4-CCL3"
  # dataset = "EL4-CCL3-SNFKL"
  # dataset = "EL4-CCL3-CCL4"
  # dataset = "EL4-CCL3-CCL4-SNFKL"
  #-----
  dataset = "CCR5KO_OT1"
  # dataset = "OT1"
} else {
  dataset <- args[1] 
}

# Several variables must be set for each specific process being analysed. These variables they are:

# Values defining the geometry of the observed region and well:
# rB; radius of the well, B, from the origin (or what is assumed to be the centre of the bolus.
# rR; radius of the observed region, R, from the origin (or what is assumed to be the centre of the bolus.

# User defined maximum and minimum times to be analysed:
# time_index_min = -1; Mimimum time iteration to examine. -1 for all available data.
# time_index_max = -1; Maximum time iteration to examine. -1 for all available data.

# The resolution at which r is discretised for use in analysis.
# r_resolution = 1

# The number of splines used to transform empirical step K functions to smoothed functional data objects:
# nsplines = NULL; The user can choose to set a specific number of splines to use in smoothing. Set to NULL for automatic slection of splines.

# select_nsplines_prop = T; We fined from our analysis of simple simulations that choosing the number of splines based on the length of A works well. Set this to T to enable this..
# If the above is set to F, nsplines is estimated via minimising the GCV measure, however, we find this leads to very inacurate results.

# nsplines_as_prop_rR = 0.1;  We fine 0.1 x rR works well if selecting number of splines based on rR.

# The order of the smoothing parameter used to smooth K functions:
# lambda = NULL; If user does not specific, lambda is selected by minimising gcv.

# The "fda" package can be quite fragile with respect the input data.
# If you receive errors relating to "negative eigenvalues", then you are running into problems.
# Selecting a smaller/finer resolution can prevent these errors.
# Failing this, we have tried to detect what features of input data case these issues, and it appears that regions "r'
# of very low (or even zero) agent density are responsible.
# These regions (at the extremes of r, not in the centre of it) can be excluded using the variables abs_rmin and abs_rmax.
# However, this excludes those regions from ALL time points.
# Given the dynamic processes we analyse, with agent repulsion and attraction, it would be preferable for the time-variant
# mode of operation to exclude these regions ONLY for the time windows where they cause problems.
# We have been unable to predict in advance when fda will cause errors.
# We've looked at truncating r based on agent
# density, but to no avail. We have left the functionality to perform this online filtering for a given window in place,
# though currently disabled, in case other users wish to experiment.


# Finding common ranges (time, r, lambda) across the three representative experiments. 
# Cog-preembed31817
# rB <- 1099 # (1082.5+1115.5)/2  # Defines the radius of the bolus
# rR <- 2985 # Shortest distance to the edge of the well
# time samples labelled in range [1, 163]
# time step = 363.3066s
# max experimental time = 163 * 363.3066 = 59218.98s ~ 16:26:59 
# Range of values on lambda_r_surface = [-300, 1000]
# 
# Cog6917
# rB <- 1067.25 # (1024.5+1110)/2 
# rR <- 2578 # Shortest distance to the edge of the well
# time samples labelled in range [1, 155]
# time step = 373.381s
# max experimental time = 155 * 373.381 = 57874.05s ~ 16:04:34 
# Range of values on lambda_r_surface = [-154.1129  951.5368]
#
# Non-cognate8917
# rB <- 1076.75 # = (1043+1110.5)/2 
# rR <- 2855 # Shortest distance to the edge of the well
# time samples labelled in range [1, 155]
# time step = 373.47s
# max experimental time = 155 * 373.47s = 57887.85s ~ 16:04:48 
# Range of values on lambda_r_surface = [-300, 500]
# 
# Simulations
# Sim 20181211-well-chemotaxisSA/#57 Range of values on lambda_r_surface = [-100, 1600]
# Sim 20181211-well-chemotaxisSA/#57 Range of values on lambda_r_surface = [-100, 2000]

lambda_n = 10^4 #4 Smoothing parameter used to smooth the number of agents in R. Used to analyse birth death processes.
splines_n = 100

r_resolution = 1  # The resolution at which r is discretised for use in analysis. Increasing this increases the resolution used to define K functions.

nsplines = NULL
if(is.null(nsplines)) {
  select_nsplines_prop = T
  nsplines_as_prop_rR = 0.1
}

lambda = NULL


if (dataset == "human_control" || dataset == "human_aCD3CD28") {
  if (dataset == "human_control")  {
    rB <- (1199.5 + 1202) / 2  # Defines the radius of the bolus. Shortest and longest distances, take mean. 
    rR <- 2981 # Shortest distance to the edge of the well from the origin.
    # Timestep is 5:00.020074 ; 195 time observations in this dataset. 
    time_step_sec = (5 * 60) + 0.020074  # = 300.0201 Length of time step in seconds.
    
    working_dir <- paste('~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/human_T_cells/Control', sep="")
    
  }  else if (dataset == "human_aCD3CD28")  {
    rB <- (1299.5+1302.5)/2  # Defines the radius of the bolus
    rR <- 2517
    # Timestep is 5:00.020128 ; 195 time observations in this dataset. 
    time_step_sec = (5 * 60) + 0.020128  # = 300.0201 Length of time step in seconds.
    
    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/human_T_cells/aCD3CD28'
  }
  r_plot_min = -1000  # These are interpreted in relation to how "r" is to be plotted - with r=0 at the bolus edge.
  r_plot_max = 1200  
  time_index_min = 1
  time_index_max = 168  # Smallest experimental time point over 16h, given the time step size. Calculations above. 
  kymograph_zlim = c(-300, 2300)  # Override generic case above to be consistent across the two "human" experiments. 
}

if (dataset == "Cog-preembed31817" || dataset == "Cog6917" || dataset == "Non-cognate8917") { 
  if (dataset == "Cog-preembed31817") {
    rB <- (1082.5 +1115.5) / 2  # Defines the radius of the bolus
    rR <- 2985 # Shortest distance to the edge of the well
    # = 363.3066 Length of time step in seconds.
    time_step_sec = (6 * 60) + 03.306607  
    # 163 time observations in this dataset. 
    time_index_min = 1
    time_index_max = 158  # Smallest experimental time point over 16h, given the time step size. Calculations above. 
   
    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/2017-08-31 cognate pre-embedded - for Mark'
    
  } else if (dataset == "Cog6917")  {
    rB <- (1024.5+1110)/2
    rR <- 2578

    # 155 time observations in this dataset (labelled from 1 to 155)
    time_step_sec = (6 * 60) + 13.380961
    time_index_min = 1
    time_index_max = 155  # Smallest experimental time point over 16h, given the time step size. Calculations above. 
  
    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/2017-09-06 Cognate - Data for Mark'
    
  } else  if (dataset == "Non-cognate8917") {
    rB <- (1043+1110.5)/2
    rR <- 2855
    # 155 time observations in this dataset (labelled from 1 to 155)
    time_step_sec = (6 * 60) + 13.470
    
    time_index_min = 1
    time_index_max = 155  # Smallest experimental time point over 16h, given the time step size. Calculations above.   
    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/2017-09-08 Non-cognate - for Mark'
  }
  
  r_plot_min = -1000  # These are interpreted in relation to how "r" is to be plotted - with r=0 at the bolus edge.
  r_plot_max = 1500
  kymograph_zlim = c(-300, 1550)
}

if (dataset == "EL4-CCL3" || dataset == "EL4-CCL3-SNFKL") { 
  if (dataset == "EL4-CCL3")  
  {
    rB <- (851+856.5)/2
    rR <- 2956
    # 4:58.783751 ; 177 time observations in this dataset (labelled from 1 to 177)
    time_step_sec = (4 * 60) + 58.783751
    
    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/EL4-CCL3 with and without SNFKL/EL4-CCL3'
  } else if (dataset == "EL4-CCL3-SNFKL")  # REPRESENTATIVE OF EXPERIMENT
  {
    rB <- (1199.5 + 1202) / 2
    rR <- 3001
    # 4:58.832250
    time_step_sec = (4 * 60) + 58.832250
  
    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/EL4-CCL3 with and without SNFKL/EL4-CCL3-SNFKL'
  }
  
  r_plot_min = -750  # r_plot_min and r_plot_max values are expressed w.r.t. the bolus edge (= 0).
  r_plot_max = 1550
  kymograph_zlim = c(-100, 1250)  # Override generic case above to be consistent across the two "human" experiments. 
  
  time_index_max = -1
  time_index_min = -1  # Smallest experimental time point over 16h, given the time step size. Calculations above. 
}

if (dataset == "EL4-CCL3-CCL4" || dataset == "EL4-CCL3-CCL4-SNFKL") {
  if (dataset == "EL4-CCL3-CCL4") 
  {
    rB <- (1018.5 + 1112.5)/2
    rR <- 2558
    # 4:58.844494 ; 177 time observations in this dataset (labelled from 1 to 177
    time_step_sec = (4 * 60) + 58.844494
    
    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/EL4-CCL3-CCL4 with and without SNFKL/EL4-CCL3-CCL4'
    
  } else if (dataset == "EL4-CCL3-CCL4-SNFKL")
  {
    rB <- (1257.5 + 1328) / 2
    rR <- 2806
    # 5:00:020585 ; 177 time observations in this dataset (labelled from 1 to 177
    time_step_sec = (5 * 60) + 00.020585

    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/EL4-CCL3-CCL4 with and without SNFKL/EL4-CCL3-CCL4-SNFKL'
  }
  
  r_plot_min = -900  # These are interpreted in relation to how "r" is to be plotted - with r=0 at the bolus edge.
  r_plot_max = 1300
  kymograph_zlim = c(-100, 1000)
  
  time_index_min = -1
  time_index_max = -1 
}


if (dataset == 'CCR5KO_OT1' || dataset == 'OT1')
{
  # OT1 and CCR5KO were imaged within the same well experiment, hence they have the same values for timestep, rR and rB.
  if (dataset == 'CCR5KO_OT1') {
    rB <- (1021.5 + 1135.5) / 2  # Defines the radius of the bolus
    rR <- 2465  # Defines the radius of the observed area. Shortest distance between bolus centre and well perimeter.  
    time_step_sec = (5 * 60) + 0.02  # Length of time step in seconds.
    
    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/2018-11-09_cognate_OT1_CCR5KO-OT1/CCR5KO_OT1'
  } else if (dataset == 'OT1') {
    rB <- (1021.5 + 1135.5) / 2  # Defines the radius of the bolus
    rR <- 2465  # Defines the radius of the observed area. Shortest distance between bolus centre and well perimeter.  
    time_step_sec = (5 * 60) + 0.02  # Length of time step in seconds.
    
    working_dir <- '~/dropbox_usyd/projects/biro/neutroswarm/imaris/whole_well/2018-11-09_cognate_OT1_CCR5KO-OT1/OT1'
  }
  
  r_plot_min = -1000  # These are interpreted in relation to how "r" is to be plotted - with r=0 at the bolus edge.
  r_plot_max = 1300
  kymograph_zlim = c(-100, 700)
  
  time_index_min = -1
  time_index_max = -1   
}


# --------------
# Read and define the "distance to bolus" files:
dtb_file <- paste(working_dir, "/distances_to_bolus_R.csv", sep='')

if (!file.exists(dtb_file)) {
  print(dtb_file)
  stop('specified distance to bolus file does not exist')
}

# Distance to bolus. Table must contain two columns, distance and then times at which those distances were observed.
# Distance to bolus = 0 is the bolus edge. Negative values are within the bolus, from the edge. 
dtb_mat <- read.csv(dtb_file,
                    skip=0,
                    colClasses=c(NA,NA))
dtb_mat <- dtb_mat[order(dtb_mat$Time), ]  # Sort by monotonically increasing time.

# Truncate time as required by user (there are some sections of experiments that may be unrepresentative)
time_step_indices <- unique(dtb_mat$Time)
if (! exists('time_index_max') || time_index_max == -1) { time_index_max <- length(time_step_indices)}
if (! exists('time_index_min') || time_index_min == -1) { time_index_min <- 1}
dtb_mat <- dtb_mat[which((dtb_mat$Time >= time_step_indices[time_index_min]) & (dtb_mat$Time <= time_step_indices[time_index_max])),]
time_step_indices <- unique(dtb_mat$Time)

# If not otherwise specified, draw full range of available data. 
if (! exists('r_plot_min') || r_plot_min == -1) { r_plot_min = -rB }
if (! exists('r_plot_max') || r_plot_max == -1) { r_plot_max = rR - rB}

source('~/dropbox_usyd/projects/biro/neutroswarm/imaris/swarming_metric_biro/swarm_metric_commons.R')
result = swarm_metric_preprocessing(dtb_mat, time_step_indices=time_step_indices, 
                                    rB, rR, r_resolution, 
                                    time_step_sec=time_step_sec, 
                                    nsplines=NULL,
                                    lambda=NULL,
                                    nsplines_as_prop_rR=0.1,
                                    select_nsplines_prop=TRUE,
                                    min_loglam=-3, max_loglam=10,
                                    working_dir=working_dir,
                                    Norder=6, Lfdval=4,  # Used for smoothing with B splines. 
                                    r_plot_min=r_plot_min, r_plot_max=r_plot_max,  # Limits for plotting r 
                                    kymograph_zlim=kymograph_zlim)  
print('range of lambda values = ')
print(range(result$lambda_r_plot))
