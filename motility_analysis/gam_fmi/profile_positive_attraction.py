"""
Execute from parent directory, as:
$> python3 -m imaris.profile_positive_attraction.py [-ac 0.8]

"""
import motility_analysis.profile as profile
import motility_analysis.graphing as graphing
from imaris.profile_attraction_helper import process
from imaris.profile_attraction_helper import build_profile_handler
from imaris.profile_attraction_helper import chemotactic_response_handler
import os
import sys


outdir = 'imaris/positive_attraction/'

d1 = outdir + '/2016_03_16_Positive_control_2'
d2 = outdir + '/2016_05_03_Positive_control_1'
d3 = outdir + '/2016_05_26_Positive_control_3'
d4 = outdir + '/2017-03-30_IgG_3'
d5 = outdir + '/2017-05-05_IgG_3'
d6 = outdir + '/2017-07-06_Positive_1'

arrest_coefficient_filter = 0.8

if '-ac' in sys.argv:
    arrest_coefficient_filter = float(sys.argv[sys.argv.index('-ac') + 1])

outdir = outdir + '/arrest_coeff_filter_{:.2f}'.format(arrest_coefficient_filter)
os.makedirs(outdir, exist_ok=True)


cog_pre_1 = build_profile_handler(d1, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
chemotactic_response_handler(d1, cog_pre_1, arrest_coefficient_filter=arrest_coefficient_filter)

cog_pre_2 = build_profile_handler(d2, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
chemotactic_response_handler(d2, cog_pre_2, arrest_coefficient_filter=arrest_coefficient_filter)

cog_pre_3 = build_profile_handler(d3, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
chemotactic_response_handler(d3, cog_pre_3, arrest_coefficient_filter=arrest_coefficient_filter)

cog_pre_4 = build_profile_handler(d4, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
chemotactic_response_handler(d4, cog_pre_4, arrest_coefficient_filter=arrest_coefficient_filter)

cog_pre_5 = build_profile_handler(d5, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
chemotactic_response_handler(d5, cog_pre_5, arrest_coefficient_filter=arrest_coefficient_filter)

cog_pre_6 = build_profile_handler(d6, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
chemotactic_response_handler(d6, cog_pre_6, arrest_coefficient_filter=arrest_coefficient_filter)


all_cog_pre = [cog_pre_1, cog_pre_2, cog_pre_3, cog_pre_4, cog_pre_5, cog_pre_6]
chemotactic_response_handler(outdir+'/summary_chemotaxis', all_cog_pre, arrest_coefficient_filter=arrest_coefficient_filter)


# read in each of the data files
# graphs = False
# trim_displacement = None
# fmi_direction = (-1, 0, 0)
# p1 = profile.build_profile(directory=d1, graphs=graphs, trim_displacement=trim_displacement, trim_observations=None,
#                            trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
#                            analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
# p2 = profile.build_profile(directory=d2, graphs=graphs, trim_displacement=trim_displacement, trim_observations=None,
#                            trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
#                            analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
# p3 = profile.build_profile(directory=d3, graphs=graphs, trim_displacement=trim_displacement, trim_observations=None,
#                            trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
#                            analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
# p4 = profile.build_profile(directory=d4, graphs=graphs, trim_displacement=trim_displacement, trim_observations=None,
#                            trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
#                            analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
# p5 = profile.build_profile(directory=d5, graphs=graphs, trim_displacement=trim_displacement, trim_observations=None,
#                            trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
#                            analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
# p6 = profile.build_profile(directory=d6, graphs=graphs, trim_displacement=trim_displacement, trim_observations=None,
#                            trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
#                            analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
# all_prof = [p1, p2, p3, p4, p5, p6]

# calculate forward migration index w.r.t. negative x axis
# for prof in all_prof:
#     for tra in prof.tracks:
#         tra._calculate_forward_migration_index(fmi_vector=[-1., 0., 0.])

# graph_dir = outdir+'/summary'
# os.makedirs(graph_dir, exist_ok=True)
# graphing.plot_profile_graphs(all_prof,  directory=graph_dir)
# graphing.plot_msd_vs_max_dt(all_prof, directory=graph_dir)

# graphing.plot_track_fmi_biases(all_prof, directory=graph_dir)
# graphing.plot_track_fmi(all_prof, directory=graph_dir)
# graphing.plot_instantaneous_fmi(all_prof, directory=graph_dir)

# process(all_prof, outdir)
