# Applies the swarm metric to motilisim-generated "well" data. 
#
# 20181128: Note that the brth_dth=T process calculates f on the correct ranges. However, it screws up M. 
# For M you need to run brth_dth=F. 
#
# Mark N. Read, 2018

# Clean the R environment. 
rm(list=ls())

require(fda)
require(rainbow)
require(fields)

setwd("/Users/markread/dropbox_usyd/projects/biro/neutroswarm/imaris/swarming_metric_biro/")

######### Set up the analysis ##########

# Point R towards the data to be analysed. User supplies directory through command line argument, or typing it in below.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  working_dir = "/Users/markread/dropbox_usyd/projects/biro/neutroswarm/imaris/simulation/20181212-well-unbiased/rep2/"
  
  # working_dir <- "~/dropbox_usyd/projects/biro/neutroswarm/imaris/simulation/20181211-well-chemotaxisSA/data/57/rep0"
} else {
  # Directory where data is contained, must be the firt argument. 
  # This is best supplied as an absolute path - it will work with relative paths, but needs to be relative to where this
  # script is located on the filesystem, which can be complicated to figure out. 
  working_dir <- args[1] 
}

rB = 1200  # Radius of bolus
rR = 3400  # Radius of well
rA = rR-rB  # Distance from bolus to well perimeter

# These are values used for plotting in the analysis; set to be common for all three representative experiments. 
# These are interpreted in relation to how "r" is to be plotted - with r=0 at the bolus edge, or at the origin?
r_plot_min = -1000
r_plot_max = 1500
kymograph_zlim = c(-300, 1550) #c(-200, 1600)

time_step_sec = 300
time_step_min = time_step_sec / 60

time_index_min = -1  # Mimimum time iteration (not absolute time stamp) to examine. -1 for all available data.
time_index_max = -1  # Maximum time iteration (not absolute time stamp) to examine. -1 for all available data.

r_resolution = 1  # Resolution for analysis across radius
nsplines = NULL

lambda_n = 10^4  # Smoothing parameter used to smooth the number of agents in R. Used to analyse birth death processes.
splines_n = 100


######### Start the analysis ##########  
  
dtb_file = paste(working_dir, "/distances_to_bolus.csv", sep='')
if (! file.exists(dtb_file)) {
  print(dtb_file)
  stop('Specified distance to bolus file does not exist')
}


# Distance to bolus. Table must contain two columns, distance and then times at which those distances were observed
dtb_mat = read.csv(dtb_file,
                   skip=0,
                   colClasses=c(NA, NA))
dtb_mat = dtb_mat[order(dtb_mat$Time),]  # Sort by monotonically increasing time.
# Downstream analysis expects time_step_indices to be time steps. However upstream scripts express time in minutes. Convert. 
dtb_mat$Time = dtb_mat$Time / time_step_min  
# NOTE imaris starts counting time from t=1, but t=0 is more standard. Convert here so first time = t=1
dtb_mat$Time = dtb_mat$Time + 1

# Vector of unique times observed. Should be in seconds.
# Represents the raw data, user-specified subsetting thereof occurs below. 
time_step_indices = unique(dtb_mat$Time)

# If not otherwise specified, draw full range of available data. 
if (time_index_max == -1)
{ time_index_max = length(time_step_indices) }
if (time_index_min == -1)
{ time_index_min = 1 }

# If not otherwise specified, draw full range of available data. 
if (r_plot_min == -1)  
{ r_plot_min = -rB}
if (r_plot_max == -1)
{ r_plot_max = rR - rB}

dtb_mat = dtb_mat[which(   (dtb_mat$Time>=time_step_indices[time_index_min])    # Extract user-selected subset of times
                         & (dtb_mat$Time<=time_step_indices[time_index_max]) )
                  , ]  # All columns
time_step_indices <- unique(dtb_mat$Time)  # Update to reflect possible user-selected restrictions 


source('swarm_metric_commons.R')
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

# if (FALSE) {
#   if (TVAR) {
#     print("Time-varying analysis")
#     # Typical in functional data anlaysis to choose this as sqrt(# time points)
#     half_time_window_width = floor(sqrt(TMAX)) 
#     if(brth_dth){
#       source('swarm_metric_tvar_brthdth.R')
#     }else{
#       source('swarm_metric_tvar.R')
#     }
#   } else {
#     print("Time-INVARIANT analysis")
#     if(brth_dth){
#       source('swarm_metric_core_Bsplines_brthdth.R')
#     }else{
#       source('swarm_metric_core_Bsplines.R')
#     }
#   }
# }
