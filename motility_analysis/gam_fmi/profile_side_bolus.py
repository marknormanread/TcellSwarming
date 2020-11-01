"""
Run as:
$> python3 -m imaris.profile_side_bolus_tracks -ac 0.8

Mark N. Read, 2018
"""
import imaris.profile_attraction_helper as profile_attraction_helper
import motility_analysis.contrast_profiles as contrast
import motility_analysis.graphing as graphing
import motility_analysis.profile as profile
import os
import sys
import seaborn as sns


def generate_profile(direc):
    fmi_direction = (-1, 0, 0)

    arrest_coefficient_filter = 1.0
    outdir = '{:s}/'.format(direc)    
    if '-ac' in sys.argv:
        arrest_coefficient_filter = float(sys.argv[sys.argv.index('-ac') + 1])
        outdir = '{:s}/arrest_coeff_filter_{:.2f}'.format(direc, arrest_coefficient_filter)
    os.makedirs(outdir, exist_ok=True)

    p1 = profile.build_profile \
        (directory=direc, graphs=False,
         trim_arrest_coefficient=arrest_coefficient_filter, tracks=None, fmi_direction=fmi_direction,
         analyse_teleports=False, interpolate=True, style="whitegrid")

    return p1, outdir


def generate_graphs(prof, graph_dir):
    if False:
        # os.makedirs(graph_dir, exist_ok=True)
        # graphing.plot_msd_errorbars(prof, directory=graph_dir)
        for i in list(range(3)):
            graphing.plot_representative_tracks(Ps=[prof], graph_path=graph_dir + '/representative_tracks_' + str(i),
                                                n=40, xlim=[-125, 125], ylim=[-100, 120],
                                                sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
                                                seed=i)
        graphing.plot_msd_allT([prof], directory=graph_dir)
        graphing.plot_dt_displacements([prof], directory=graph_dir)
        graphing.plot_displacement_autocorrelation([prof], directory=graph_dir)
        graphing.plot_speed_histo([prof], directory=graph_dir)
        graphing.plot_speed_ecdf([prof], directory=graph_dir)

        graphing.plot_turns_hist([prof], directory=graph_dir)
        graphing.plot_turns_ecdf([prof], directory=graph_dir)
        graphing.plot_instantaneous_fmi([prof], directory=graph_dir)
        graphing.plot_track_fmi([prof], directory=graph_dir)
        graphing.plot_track_fmi_biases([prof], directory=graph_dir)


def chemotactic_analysis(prof, outdir):
    if False:
        chemotactic_response = profile_attraction_helper.process(prof, outdir)
        chemotactic_response.to_csv(outdir + '/chemotactic_response.csv', index=False)


outdir = 'imaris/Side_bolus_tracks/cognate1'
cog1_prof, cog1_fp = generate_profile(outdir)
generate_graphs(cog1_prof, cog1_fp)
chemotactic_analysis(cog1_prof, outdir)
    
outdir = 'imaris/Side_bolus_tracks/cognate2'
cog2_prof, cog2_fp = generate_profile(outdir)
generate_graphs(cog2_prof, cog2_fp)
chemotactic_analysis(cog2_prof, outdir)

outdir = 'imaris/Side_bolus_tracks/cognate3'
cog3_prof, cog3_fp = generate_profile(outdir)
generate_graphs(cog3_prof, cog3_fp)
chemotactic_analysis(cog3_prof, outdir)

# All together
outdir = 'imaris/Side_bolus_tracks/cognate_summary'
cog_prof = profile.build_profile(tracks=cog1_prof.tracks + cog2_prof.tracks + cog3_prof.tracks)
generate_graphs(cog_prof, graph_dir=outdir)
chemotactic_analysis(cog_prof, outdir)


outdir = 'imaris/Side_bolus_tracks/non_cognate1'
non_cog1_prof, non_cog1_fp = generate_profile(outdir)
generate_graphs(non_cog1_prof, non_cog1_fp)
chemotactic_analysis(non_cog1_prof, outdir)

outdir = 'imaris/Side_bolus_tracks/non_cognate2'
non_cog2_prof, non_cog2_fp = generate_profile(outdir)
generate_graphs(non_cog2_prof, non_cog2_fp)
chemotactic_analysis(non_cog2_prof, outdir)

outdir = 'imaris/Side_bolus_tracks/non_cognate3'
non_cog3_prof, non_cog3_fp = generate_profile(outdir)
generate_graphs(non_cog3_prof, non_cog3_fp)
chemotactic_analysis(non_cog3_prof, outdir)
# All together
outdir = 'imaris/Side_bolus_tracks/non_cognate_summary'
non_cog_prof = profile.build_profile(tracks=non_cog1_prof.tracks + non_cog2_prof.tracks + non_cog3_prof.tracks)
generate_graphs(non_cog_prof, graph_dir=outdir)
chemotactic_analysis(non_cog_prof, outdir)


outdir = 'imaris/Side_bolus_tracks/random1'
rand1_prof, rand1_fp = generate_profile(outdir)
generate_graphs(rand1_prof, rand1_fp)
chemotactic_analysis(rand1_prof, outdir)

outdir = 'imaris/Side_bolus_tracks/random2'
rand2_prof, rand2_fp = generate_profile(outdir)
generate_graphs(rand2_prof, rand2_fp)
chemotactic_analysis(rand2_prof, outdir)

outdir = 'imaris/Side_bolus_tracks/random3'
rand3_prof, rand3_fp = generate_profile(outdir)
generate_graphs(rand3_prof, rand3_fp)
chemotactic_analysis(rand3_prof, outdir)
# All together
outdir = 'imaris/Side_bolus_tracks/random_summary'
rand_prof = profile.build_profile(tracks=rand1_prof.tracks + rand2_prof.tracks + rand3_prof.tracks)
generate_graphs(rand_prof, graph_dir=outdir)
chemotactic_analysis(rand_prof, outdir)

# Plot representative tracks
# outdir = 'imaris/Side_bolus_tracks/representative_tracks'
# os.makedirs(outdir, exist_ok=True)
# for i in list(range(10)):
#     graphing.plot_representative_tracks(Ps=[rand_prof], graph_path=outdir + '/random_representative_tracks_' + str(i),
#                                         n=40, xlim=None, ylim=None,
#                                         sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
#                                         seed=i)

# for i in list(range(10)):
#     graphing.plot_representative_tracks(Ps=[non_cog_prof], graph_path=outdir + '/non_cognate_representative_tracks_' + str(i),
#                                         n=40, xlim=None, ylim=None,
#                                         sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
#                                         seed=i)

# for i in list(range(10)):
#     graphing.plot_representative_tracks(Ps=[cog_prof], graph_path=outdir + '/cognate_representative_tracks_' + str(i),
#                                         n=40, xlim=None, ylim=None,
#                                         sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
#                                         seed=i)

# Plot representative tracks
outdir = 'imaris/Side_bolus_tracks/representative_tracks_selected'
os.makedirs(outdir, exist_ok=True)
xlim = [-140, 130]
ylim = [-130, 110]
i = 5
sns.set_style("whitegrid")  # Seems to override matplotlib settings, which are giving inconsistent results anyway. 
sns.set_context("paper")
graphing.plot_representative_tracks(Ps=[rand_prof], graph_path=outdir + '/random_representative_tracks_' + str(i),
                                    n=40, xlim=xlim, ylim=ylim,
                                    sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
                                    fontsize=18.0,
                                    seed=i)

i = 3
graphing.plot_representative_tracks(Ps=[non_cog_prof], graph_path=outdir + '/non_cognate_representative_tracks_' + str(i),
                                    n=40, xlim=xlim, ylim=ylim,
                                    sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
                                    fontsize=18.0,
                                    seed=i)

i = 7
graphing.plot_representative_tracks(Ps=[cog_prof], graph_path=outdir + '/cognate_representative_tracks_' + str(i),
                                    n=40, xlim=xlim, ylim=ylim,
                                    sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
                                    fontsize=18.0,
                                    seed=i)

# MSDs of all three conditions
# graphing.plot_msd_errorbars(profiles=[cog_prof, non_cog_prof, rand_prof], 
#                             labels=['cognate', 'non cognate', 'random'], directory='imaris/Side_bolus_tracks/contrast')


# contrast.contrast(profile1=[cog1_prof, cog2_prof, cog3_prof], 
#                   profile2=[non_cog1_prof, non_cog2_prof, non_cog3_prof], 
#                   profile3=[rand1_prof, rand2_prof, rand3_prof], 
#                   label1='cognate', label2='non cognate', label3='random', 
#                   out_dir='imaris/Side_bolus_tracks/contrast', draw_graphs=True,
#                   p1_colour='b', p2_colour='g', p3_colour='r',
#                   p1_p2_contrast_name='cognate_vs_non_cognate.xml',
#                   p1_p3_contrast_name='cognate_vs_random.xml', 
#                   settings=contrast.ContrastSettings(dt=cog1_prof.tracks[0].timestep_s))