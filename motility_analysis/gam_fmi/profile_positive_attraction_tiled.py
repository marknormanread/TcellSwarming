"""
Execute from parent directory, as:
$> python3 -m imaris.profile_positive_attraction_tiled.py [-ac 0.8]

(this is needed to import motility_analysis, which sits in the root project directory)

Script characterises the positive attraction forces attracting cells to the bolus. This is done by drawing 2D density
scatter graphs of x-position, track age, direction of approach to bolus, translational speed and turn speed with
one another.

Mark N. Read, 2017
"""
import motility_analysis.contrast_profiles as contrast
import motility_analysis.profile as profile
import motility_analysis.graphing as graphing
import matplotlib
matplotlib.rcParams['figure.autolayout'] = True
import matplotlib.pyplot as plt
import numpy as np
import os
# import seaborn as sns
import sys
import pandas
from imaris.profile_attraction_helper import process
from imaris.profile_attraction_helper import build_profile_handler


d1 = 'imaris/positive_attraction_tiled1_test/'
dt1 = 26.

d2 = 'imaris/positive_attraction_tiled2/'
dt2 = 30.12

cog1_fp = 'imaris/side_bolus_tiled/cognate1_positive_attraction_tiled1'
cog2_fp = 'imaris/side_bolus_tiled/cognate2_positive_attraction_tiled2'
cog3_fp = 'imaris/side_bolus_tiled/cognate3'
cog4_fp = 'imaris/side_bolus_tiled/cognate4'
cog5_fp = 'imaris/side_bolus_tiled/cognate5'

non_cog1_fp = 'imaris/side_bolus_tiled/non_cognate1'
non_cog2_fp = 'imaris/side_bolus_tiled/non_cognate2'
non_cog3_fp = 'imaris/side_bolus_tiled/non_cognate3'


# Read argument from command line. 
arrest_coefficient_filter = 1.0
outdir = '{:s}/'.format(d1)
if '-ac' in sys.argv:
    arrest_coefficient_filter = float(sys.argv[sys.argv.index('-ac') + 1])


# def chemotactic_response_handler(direct, profiles, arrest_coefficient_filter):
#     """ 
#     Python-based heatmaps of FMI, speeds, directionality vs distance and time. 
#     Outputs measures in a csv file for downstream analysis in R. 
#     Profiles can be either a motility_analyis.profile.Profile object, or a list of such objects. 
#     """
#     # Analyse the chemotactic response of these data sets. Can also draw graphs using python infrastructure. 
#     outdir = '{:s}/arrest_coeff_filter_{:.2f}'.format(direct, arrest_coefficient_filter)
#     os.makedirs(outdir, exist_ok=True)
#     chemotactic_response = process(profiles, outdir)
#     chemotactic_response.to_csv(outdir + '/chemotactic_response.csv', index=False)


# def motility_graphing_handler(direct, profiles, arrest_coefficient_filter):
#     outdir = '{:s}/arrest_coeff_filter_{:.2f}'.format(direct, arrest_coefficient_filter)
#     graph_dir = outdir + '/motility_summary'
#     os.makedirs(graph_dir, exist_ok=True)
#     graphing.plot_profile_graphs(profiles, directory=graph_dir)


# def representative_tracks_handler(direct, profiles):
#     xlim = None  # e.g. [-125, 125]
#     ylim = None
#     os.makedirs(direct, exist_ok=True)
#     for i in list(range(10)):
#         graphing.plot_representative_tracks(Ps=profiles, graph_path=direct + '/representative_tracks_' + str(i),
#                                             n=40, xlim=None, ylim=None,
#                                             sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
#                                             seed=i)


# def build_profile_handler(data_fp, graphs=False, trim_displacement=False, fmi_direction=(-1, 0, 0)):
#     prof = profile.build_profile(directory=data_fp, graphs=False,
#                                  trim_displacement=trim_displacement,
#                                  trim_arrest_coefficient=arrest_coefficient_filter, tracks=None, fmi_direction=fmi_direction,
#                                  analyse_teleports=False, interpolate=True, style="whitegrid")
#     if graphs:
#         motility_graphing_handler(data_fp, prof)
#     return prof

# Build a profile from data in the given directory.                                                                             
cog1 = build_profile_handler(cog1_fp, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
# graphing.output_track_start_end_points([cog1], cog1_fp)
# chemotactic_response_handler(cog1_fp, cog1)
# cog2 = build_profile_handler(cog2_fp, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
# chemotactic_response_handler(cog2_fp, cog2)
cog3 = build_profile_handler(cog3_fp, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
# graphing.output_track_start_end_points([cog3], cog3_fp)
# chemotactic_response_handler(cog3_fp, cog3)
cog4 = build_profile_handler(cog4_fp, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
# graphing.output_track_start_end_points([cog4], cog4_fp)
# chemotactic_response_handler(cog4_fp, cog4)
# cog5 = build_profile_handler(cog5_fp, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
# chemotactic_response_handler(cog5_fp, cog5)

# all_cog = [cog1, cog2, cog3, cog4, cog5]
all_cog = [cog1, cog3, cog4]
# graphing.output_track_start_end_points(all_cog, 'imaris/side_bolus_tiled/cognate134/')
# chemotactic_response_handler('imaris/side_bolus_tiled/cognate134', all_cog)
# motility_graphing_handler('imaris/side_bolus_tiled/cognate134/motility_summary', all_cog)
# representative_tracks_handler('imaris/side_bolus_tiled/cognate134/motility_summary/representative_tracks', all_cog)

non_cog1 = build_profile_handler(non_cog1_fp, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
# graphing.output_track_start_end_points([non_cog1], non_cog1_fp)
# chemotactic_response_handler(non_cog1_fp, non_cog1)
non_cog2 = build_profile_handler(non_cog2_fp, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
# graphing.output_track_start_end_points([non_cog2], non_cog2_fp)
# chemotactic_response_handler(non_cog2_fp, non_cog2)
non_cog3 = build_profile_handler(non_cog3_fp, graphs=False, arrest_coefficient_filter=arrest_coefficient_filter)
# graphing.output_track_start_end_points([non_cog3], non_cog3_fp)
# chemotactic_response_handler(non_cog3_fp, non_cog3)

all_non_cog = [non_cog1, non_cog2, non_cog3]
# graphing.output_track_start_end_points(all_non_cog, 'imaris/side_bolus_tiled/non_cognate123')
# chemotactic_response_handler('imaris/side_bolus_tiled/non_cognate123', all_non_cog)
# motility_graphing_handler('imaris/side_bolus_tiled/non_cognate123/motility_summary', all_non_cog)
# representative_tracks_handler('imaris/side_bolus_tiled/non_cognate123/motility_summary/representative_tracks', all_non_cog)

# contrast.contrast(profile1=all_cog, label1='cog', 
#                   profile2=all_non_cog, label2='non cog',
#                   out_dir='imaris/side_bolus_tiled/cog_vs_non_cog', draw_graphs=True,
#                   # User specified graphs here. Can't do dt_disp if the time step varies. 
#                   settings=contrast.ContrastSettings(dt_disp=False)
#                   )


# graph_dir = 'imaris/side_bolus_tiled/non_cognate123/motility_summary/representative_tracks'
# for i in list(range(3)):
#     graphing.plot_representative_tracks(Ps=all_non_cog, graph_path=graph_dir + '/representative_tracks_' + str(i),
#                                         n=40, xlim=[-125, 125], ylim=[-100, 120],
#                                         sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
#                                         seed=i)

if False:
    # Plot three graphs of representative tracks (randomly selected, not sorted by displacement). 
    sns.set_style("whitegrid")  # Otherwise "dark"
    sns.set_context("talk")
    
    outdir = '{:s}/arrest_coeff_filter_{:.2f}'.format(d1, arrest_coefficient_filter)
    graph_dir = outdir + '/summary'
    
    graphing.plot_representative_tracks(Ps=[p1], graph_path=graph_dir + '/representative_tracks_1',
                                        n=40, xlim=None, ylim=None,
                                        sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True)
    graphing.plot_representative_tracks(Ps=[p1], graph_path=graph_dir + '/representative_tracks_2',
                                        n=40, xlim=None, ylim=None,
                                        sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True)
    graphing.plot_representative_tracks(Ps=[p1], graph_path=graph_dir + '/representative_tracks_3',
                                        n=40, xlim=None, ylim=None,
                                        sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True)                                        
    graphing.plot_profile_graphs([p1], directory=graph_dir)


    outdir = '{:s}/arrest_coeff_filter_{:.2f}'.format(d2, arrest_coefficient_filter)
    graph_dir = outdir + '/summary'
    graphing.plot_representative_tracks(Ps=[p2], graph_path=graph_dir + '/representative_tracks_1',
                                        n=40, xlim=None, ylim=None,
                                        sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True)
    graphing.plot_representative_tracks(Ps=[p2], graph_path=graph_dir + '/representative_tracks_2',
                                        n=40, xlim=None, ylim=None,
                                        sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True)
    graphing.plot_representative_tracks(Ps=[p2], graph_path=graph_dir + '/representative_tracks_3',
                                        n=40, xlim=None, ylim=None,
                                        sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True)
    graphing.plot_profile_graphs([p2], directory=graph_dir)


    outdir = '{:s}/arrest_coeff_filter_{:.2f}'.format('imaris/positive_attraction_tiled_1+2', arrest_coefficient_filter)
    graph_dir = outdir + '/summary'
    graphing.plot_representative_tracks(Ps=cog_prof, graph_path=graph_dir + '/representative_tracks_1',
                                        n=40, xlim=None, ylim=None,
                                        sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True)
    graphing.plot_representative_tracks(Ps=cog_prof, graph_path=graph_dir + '/representative_tracks_2',
                                        n=40, xlim=None, ylim=None,
                                        sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True)
    graphing.plot_representative_tracks(Ps=cog_prof, graph_path=graph_dir + '/representative_tracks_3',
                                        n=40, xlim=None, ylim=None,
                                        sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True)
    graphing.plot_profile_graphs(cog_prof, directory=graph_dir)


if False:
    outdir = '{:s}/arrest_coeff_filter_{:.2f}'.format('imaris/positive_attraction_tiled_1+2', arrest_coefficient_filter)
    graph_dir = outdir + '/summary'
    os.makedirs(graph_dir, exist_ok=True)
    graphing.plot_profile_graphs(cog_prof, directory=graph_dir)

if False:
    chemotactic_response = process(cog_prof, outdir)
    chemotactic_response.to_csv(outdir + '/chemotactic_response.csv', index=False)
    empirical_x_max = np.max(chemotactic_response['dist_bolus'])
    empirical_x_min = np.min(chemotactic_response['dist_bolus'])
    empirical_time_max = np.max(chemotactic_response['dist_bolus'])


def fov_analysis(fmi, fmi_x, fmi_t):
    """ Perform analyses on each field of view individually. Sanity check against analysis that Jorge had run. """
    fov1_h1 = []
    fov2_h1 = []
    fov3_h1 = []
    fov4_h1 = []

    fov1_h2 = []
    fov2_h2 = []
    fov3_h2 = []
    fov4_h2 = []
    x_quarter = (empirical_x_max - empirical_x_min) / 4.
    print('Position x-width = {:.2f}'.format(x_quarter))
    pos1_max = empirical_x_min + x_quarter
    pos2_max = pos1_max + x_quarter
    pos3_max = pos2_max + x_quarter
    pos4_max = pos3_max + x_quarter
    print('Position 1 start = {:.2f}'.format(empirical_x_min))
    print('Position 1 end = {:.2f}'.format(pos1_max))
    print('Position 2 end = {:.2f}'.format(pos2_max))
    print('Position 3 end = {:.2f}'.format(pos3_max))
    print('Position 4 end = {:.2f}'.format(pos4_max))

    for f, x, t in zip(fmi, fmi_x, fmi_t):
        if x < pos1_max:
            if t < 60.:
                fov1_h1.append(f)
            else:
                fov1_h2.append(f)
        elif x < pos2_max:
            if t < 60.:
                fov2_h1.append(f)
            else:
                fov2_h2.append(f)
        elif x < pos3_max:
            if t < 60.:
                fov3_h1.append(f)
            else:
                fov3_h2.append(f)
        elif x < pos4_max:
            if t < 60.:
                fov4_h1.append(f)
            else:
                fov4_h2.append(f)

    print('pos1 h1 = ' + str(len(fov1_h1)))
    print('pos1 h2 = ' + str(len(fov1_h2)))

    print('pos2 h1 = ' + str(len(fov2_h1)))
    print('pos2 h2 = ' + str(len(fov2_h2)))

    print('pos3 h1 = ' + str(len(fov3_h1)))
    print('pos3 h2 = ' + str(len(fov3_h2)))

    print('pos4 h1 = ' + str(len(fov4_h1)))
    print('pos4 h2 = ' + str(len(fov4_h2)))


    print('median: pos1, h1. Mean = {:.2f}; Median = {:.2f} '.format(np.mean(fov1_h1), np.median(fov1_h1)))
    print('median: pos1, h2. Mean = {:.2f}; Median = {:.2f} '.format(np.mean(fov1_h2), np.median(fov1_h2)))
    print('')
    print('median: pos2, h1. Mean = {:.2f}; Median = {:.2f} '.format(np.mean(fov2_h1), np.median(fov2_h1)))
    print('median: pos2, h2. Mean = {:.2f}; Median = {:.2f} '.format(np.mean(fov2_h2), np.median(fov2_h2)))
    print('')
    print('median: pos3, h1. Mean = {:.2f}; Median = {:.2f} '.format(np.mean(fov3_h1), np.median(fov3_h1)))
    print('median: pos3, h2. Mean = {:.2f}; Median = {:.2f} '.format(np.mean(fov3_h2), np.median(fov3_h2)))
    print('')
    print('median: pos4, h1. Mean = {:.2f}; Median = {:.2f} '.format(np.mean(fov4_h1), np.median(fov4_h1)))
    print('median: pos4, h2. Mean = {:.2f}; Median = {:.2f} '.format(np.mean(fov4_h2), np.median(fov4_h2)))

    # prepare data in long-form, as pandas dataframe
    pos = ['pos1'] * (len(fov1_h1) + len(fov1_h2)) + \
          ['pos2'] * (len(fov2_h1) + len(fov2_h2)) + \
          ['pos3'] * (len(fov3_h1) + len(fov3_h2)) + \
          ['pos4'] * (len(fov4_h1) + len(fov4_h2))

    f = fov1_h1 + fov1_h2 + fov2_h1 + fov2_h2 + fov3_h1 + fov3_h2 + fov4_h1 + fov4_h2

    hrs = ['H1'] * len(fov1_h1) + ['H2'] * len(fov1_h2) + \
          ['H1'] * len(fov2_h1) + ['H2'] * len(fov2_h2) + \
          ['H1'] * len(fov3_h1) + ['H2'] * len(fov3_h2) + \
          ['H1'] * len(fov4_h1) + ['H2'] * len(fov4_h2)

    d = pandas.DataFrame.from_dict({'Position': pos, 'Time': hrs, 'FMI': f})
    plt.clf()
    sns.set_style("whitegrid")
    ax = sns.swarmplot(x="Position", y="FMI", hue="Time", data=d, palette="Set2", dodge=True)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2])
    ax = sns.boxplot(x="Position", y="FMI", hue="Time", data=d, palette="Set2")
    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .0))

    plt.savefig(outdir + '/positional_fmi_swarmplot.png', transparent=True, dpi=300)


# fov_analysis(fmi, fmi_x, fmi_t)
