"""
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Mark N. Read, 2017.


Utility functions that graph data captured in one or more Profiles.

Mark N. Read, 2017
"""
import math
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
# If you are generating graphs on a system with no window manager (for instance, the Sun grid engine if using a cluster)
# then you may need to uncomment this line. If you do, Profile.plot_tracks isn't going to work.
# matplotlib.use('Agg')   # Force matplotlib to not use any Xwindows backend. Must be called before any pyplot import.
matplotlib.rcParams['figure.autolayout'] = True
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import random
import scipy.stats
from statsmodels.nonparametric.smoothers_lowess import lowess
from . import geometry
from . import plot_utils
from . import profile
import pandas
from operator import itemgetter
from itertools import *


def plot_profile_graphs(ps, directory, style="dark", histograms=False, bias_testing_graphs=False):
    """
    Draws graphs of statistical analyses based on the supplied Profile objects.
    :param ps: either a single Profile, or a list of them
    :param directory: where to write graphs to on file system
    """
    sns.set_style(style)
    sns.set_context("talk")
    if isinstance(ps, profile.Profile):
        ps = [ps]
    if not os.path.exists(directory):
        os.mkdir(directory)

    print('plot_speed_ecdf')
    plot_speed_ecdf(ps, directory=directory)
    print('plot_median_speed_ecdf')
    plot_median_speed_ecdf(ps, directory=directory)
    print('plot_irq_speeds_ecdf')
    plot_irq_speeds_ecdf(ps, directory=directory)
    print('plot_turns_ecdf')
    plot_turns_ecdf(ps, directory=directory)
    print('plot_median_turn_ecdf')
    plot_median_turn_ecdf(ps, directory=directory)
    print('plot_rolls_ecdf')
    plot_rolls_ecdf(ps, directory=directory)
    print('plot_duration_ecdf')
    plot_duration_ecdf(ps, directory=directory)
    print('plot_length_ecdf')
    plot_length_ecdf(ps, directory=directory)
    print('plot_displacement_ecdf')
    plot_displacement_ecdf(ps, directory=directory)
    print('plot_jump_lengths_ecdf')
    plot_jump_lengths_ecdf(ps, directory=directory)
    print('plot_arrest_coefficients')
    plot_arrest_coefficients(ps, directory=directory)
    print('plot_meandering_ecdf')
    plot_meandering_ecdf(ps, directory=directory)
    print('plot_representative_tracks')
    plot_representative_tracks(ps, directory=directory)
    print('plot_turn_vs_speed')
    plot_turn_vs_speed(ps, directory=directory)

    print('plot_track_displacement_angle_xy')
    plot_track_displacement_angle_xy(ps, directory)
    # Graphs generated here can be illegible and very computationally expensive to generate for datasets with large
    # numbers of tracks. >40 and it's not really worth doing. 
    # plot_track_intantaneous_msd(ps, directory=directory, plot_heatmaps=False, plot_tracks_2D=False, 
    #                             manual_inspect_tracks=False, plot_instantaneous_msd_timeseries=False)
    print('plot_msd_errorbars')
    plot_msd_errorbars(ps, directory=directory)
    print('plot_msd_allT')
    plot_msd_allT(ps, directory=directory)
    print('plot_msd_vs_max_dt')
    plot_msd_vs_max_dt(ps, directory=directory)
    print('plot_dt_displacements')
    plot_dt_displacements(ps, directory=directory)
    print('plot_displacement_autocorrelation')
    plot_displacement_autocorrelation(ps, directory=directory)
  
  

    
    print('plot_instantaneous_fmi')
    plot_instantaneous_fmi(ps, directory=directory)
    print('plot_track_fmi')
    plot_track_fmi(ps, directory=directory)
    
    print('analyse_cell_entries_imaging_volume')
    for p in ps:
        p.analyse_cell_entries_imaging_volume()
    # plot_tracks(ps)  # Shows the plot, requires user interaction to close it. 

    # Primarily used for testing and detecting biases in metrics. 
    if bias_testing_graphs:
        print('plot_mi_time')
        plot_mi_time(ps, directory=directory)
        print('plot_turns_vs_timestep')
        plot_turns_vs_timestep(ps, directory=directory)
        print('plot_medSpd_observations')
        plot_medSpd_observations(ps, directory=directory)
        print('plot_medTurn_observations')
        plot_medTurn_observations(ps, directory=directory)
        print('plot_length_displacement_scat')
        plot_length_displacement_scat(ps, directory=directory, graph=True)
        print('plot_track_fmi_biases')
        plot_track_fmi_biases(ps, directory=directory)
        print('plot_medSpd_irqSpd_scat')
        plot_medSpd_irqSpd_scat(ps, directory)

    # Cumulative distribution plots are way more informative than histograms. 
    if histograms:
        plot_speed_histo(ps, directory=directory)
        plot_duration_histo(ps, directory=directory)
        plot_displacement_hist(ps, directory=directory)
        plot_median_turn_hist(ps, directory=directory)
        plot_length_histo(ps, directory=directory)
        plot_jump_lengths_hist(ps, directory=directory)
        plot_meandering_hist(ps, directory=directory)
        plot_median_speed_hist(ps, directory=directory)
        plot_irq_speeds_hist(ps, directory=directory)
        plot_turns_hist(ps, directory=directory)
        plot_rolls_hist(ps, directory=directory)
        plot_median_roll_hist(ps, directory=directory)


def output_track_start_end_points(profiles, directory):
    """ Writes a CSV of each track's start and end points, and the angle (arctan2) describing that displacement. """
    records = []  # List of dictionaries. 
    for prof in profiles:
        for track in prof.tracks:
            record = { 'track_id': track.ID,
                       'start_x': track.positions[0].x,
                       'start_y': track.positions[0].y,
                       'start_z': track.positions[0].z,
                       'end_x': track.positions[-1].x,
                       'end_y': track.positions[-1].y,
                       'end_z': track.positions[-1].z,
                       'dx': track.displacement_vector[0],
                       'dy': track.displacement_vector[1],
                       'dz': track.displacement_vector[2],
                       'arctan2(dx,dy)': track.get_arctan_xy()
                       }
            records.append(record)
    df = pandas.DataFrame(records)    
    df.to_csv(directory + '/track_start_end_points.csv', index=False)


def plot_track_displacement_angle_xy(profiles, directory, bins=36):
    """ Rose plot showing the orientation of track displacement in the xy plane as angles. """
    angles = []
    for prof in profiles:
        for track in prof.tracks:
            angles.append(track.get_arctan_xy())  # Returns radians
    
    # re: endpoint, need last bin location for np.histogram
    # bins + 1 here because this specifies boundary locations, of which there are one more than the number of bins 
    # (last bin has start AND end location)
    bin_locations = np.linspace(0.0, 2 * np.pi, bins + 1, endpoint=True)  
    bin_width = (2 * np.pi) / bins
    
    bin_hist, _ = np.histogram(angles, bin_locations)
    bin_hist = bin_hist / sum(bin_hist)  # Turn into a probability density
    
    bottom = 0  # Bin does not need to start from zero. Setting this to non-zero creates a "donut"

    ax = plt.subplot(111, polar=True)
    bin_start_theta = bin_locations[:-1]  # Exclude the last location. 
    ax.bar(bin_start_theta, bin_hist, width=bin_width, bottom=bottom, edgecolor='black')
    ax.grid(True)
    ax.set_yticklabels([])
    ax.set_axisbelow(True)  # Draw grid behind bars. 
    
    plt.savefig(directory + '/track_displacement_angle_xy.png', dpi=300)
    plt.savefig(directory + '/track_displacement_angle_xy.eps')


def plot_tracks(ps):
    """
    Plot all tracks on a 3D plot, and display. Is not (currently) saved automatically, because a 3D plot would need
    to be reoriented by the user to highlight what they wish. The start of tracks are coloured blue and finish red.
    """
    all_tracks = [t for p in ps for t in p.tracks]  # All tracks from all profiles
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    all_x = []
    all_y = []
    all_z = []
    for t in all_tracks:
        for p in t.positions:
            if p.tracked():
                all_x.append(p.x)
                all_y.append(p.y)
                all_z.append(p.z)
    max_x = max(all_x)
    max_y = max(all_y)
    max_z = max(all_z)
    min_x = min(all_x)
    min_y = min(all_y)
    min_z = min(all_z)
    # Draw the boundaries of the imaging volume.
    xs = [min_x, max_x, max_x, min_x, min_x]
    ys = [min_y, min_y, max_y, max_y, min_y]
    zs = [min_z, min_z, min_z, min_z, min_z]
    ax.plot(xs, ys, zs, c='k')
    xs = [min_x, max_x, max_x, min_x, min_x]
    ys = [min_y, min_y, max_y, max_y, min_y]
    zs = [max_z, max_z, max_z, max_z, max_z]
    ax.plot(xs, ys, zs, c='k')
    ax.plot([min_x, min_x], [min_y, min_y], [min_z, max_z], c='k')
    ax.plot([max_x, max_x], [min_y, min_y], [min_z, max_z], c='k')
    ax.plot([min_x, min_x], [max_y, max_y], [min_z, max_z], c='k')
    ax.plot([max_x, max_x], [max_y, max_y], [min_z, max_z], c='k')
    # Plot the tracks

    for track in all_tracks:
        xs = [p.x for p in track.positions if p.tracked()]
        ys = [p.y for p in track.positions if p.tracked()]
        zs = [p.z for p in track.positions if p.tracked()]
        n = len(xs)
        # use this to colour each trajectory on the path differently, to see direction of travel
        colors = unique_colour_per_series(n)
        for i in list(range(len(xs)-1)):
            ax.plot(xs[i:i+2], ys[i:i+2], zs[i:i+2], c=colors[i])

    teleport_starts = [ts for p in ps for ts in p.teleport_starts]
    teleport_ends = [ts for p in ps for ts in p.teleport_ends]
    for point in teleport_starts:
        ax.plot([point.x], [point.y], [point.z], c='b', marker='o')

    for point in teleport_ends:
        ax.plot([point.x], [point.y], [point.z], c='r', marker='o')
    plt.show()


def plot_speed_histo(profiles, directory):
    speeds = []
    for p in profiles:
        speeds.extend(p.collate_speeds())
    plot_utils.plot_histogram(data=speeds,
                              filename=os.path.join(directory, 'speeds-hist'),
                              xlabel='Speeds (um/min)', ylabel='Frequency', title='Cell Speeds')


def plot_speed_ecdf(profiles, directory, xmax=None):
    speeds = []
    for p in profiles:
        speeds.extend(p.collate_speeds())
    plot_utils.plotCDFs([speeds], xlabel='Speeds (um/min)', title='Collated Speeds',
                        filename=os.path.join(directory, 'speeds-ecdf'), xmax=xmax)


def write_distribution_medians(profiles, directory):
    speeds = []
    turns = []
    mi = []  # meandering index
    displacements = []
    for p in profiles:
        speeds.extend(p.collate_speeds())
        turns.extend(p.collate_turns())
        mi.extend(p.collate_meanders())
        displacements.extend(p.collate_displacements())
    with open(os.path.join(directory, 'median_distribution_values.txt'), 'w') as f:
        f.write('Median values of all speeds and turns across all tracks and times collated together\n')
        f.write('Median translational speed = {:.2f} microns/min\n'.format(np.median(speeds)))
        f.write('Median turn speed = {:.2f} degrees/min\n'.format(np.median(turns)))
        f.write('Median meandering index = {:.2f}\n'.format(np.median(mi)))
        f.write('Median displacement = {:.2f}\n'.format(np.median(displacements)))


def plot_individual_track_speed_hists(profiles, directory):
    direc = os.path.join(directory, 'individual_track_speeds')
    if not os.path.exists(direc):
        os.makedirs(direc)
    tID = 0
    for p in profiles:
        for t in p.tracks:
            graphLoc = os.path.join(direc, 'track{:d}'.format(tID))
            print(graphLoc)
            plot_utils.plot_histogram(data=t.speeds,
                                      filename=graphLoc,
                                      xlabel='Speeds (um/min)',
                                      ylabel='Frequency',
                                      title='Track Speeds; MI=' + str(t.meander))
            tID = tID+1


def plot_duration_histo(profiles, directory):
    duration = []
    for p in profiles:
        duration.extend([t.duration_min for t in p.tracks])
    plot_utils.plot_histogram(data=duration,    # convert durations in seconds to minutes,
                              filename=os.path.join(directory, 'durations-hist'),
                              xlabel='Duration (minutes)',
                              ylabel='Frequency',  title='Track Durations')


def plot_duration_ecdf(profiles, directory):
    duration = [t.duration_min for p in profiles for t in p.tracks]
    plot_utils.plotCDFs([duration], xlabel='Duration (minutes)', title='Track Durations',
                        filename=os.path.join(directory, 'durations-ecdf'))


def plot_length_histo(profiles, directory):
    # get track lengths. Note that broken tracks do not have a length assigned.
    length = [t.length for pro in profiles for t in pro.tracks if t.length]
    plot_utils.plot_histogram(data=length,
                              filename=os.path.join(directory, 'lengths-hist'),
                              xlabel='Length (micrometers)',
                              ylabel='Frequency', title='Track Lengths')


def plot_length_ecdf(profiles, directory):
    # get track lengths. Note that broken tracks do not have a length assigned.
    length = [t.length for pro in profiles for t in pro.tracks if t.length]
    plot_utils.plotCDFs([length], xlabel='Length (micrometers)', title='Track Lengths',
                        filename=os.path.join(directory, 'lengths-ecdf'))


def plot_jump_lengths_hist(profiles, directory):
    # all the individual jumps. Note that not all positions have displacements assigned (broken tracks).
    l = [pos.displacement for pro in profiles for t in pro.tracks for pos in t.positions if pos.displacement]
    plot_utils.plot_histogram(data=l,
                              filename=os.path.join(directory, 'jumps-hist'),
                              xlabel='Jumps (micrometers)',
                              ylabel='Frequency', title='Displacements in each time-step')


def plot_jump_lengths_ecdf(profiles, directory):
    # all the individual jumps. Note that not all positions have displacements assigned (broken tracks).
    l = [pos.displacement for pro in profiles for t in pro.tracks for pos in t.positions if pos.displacement]
    plot_utils.plotCDFs([l], xlabel='Jumps (micrometers)', title='Jump Lengths',
                        filename=os.path.join(directory, 'jumps-ecdf'))


def plot_arrest_coefficients(profiles, directory):
    arcoef = [t.arrest_coefficient for p in profiles for t in p.tracks if t.arrest_coefficient is not None]
    plot_utils.plotCDFs([arcoef], xlabel='Arrest coefficient', filename=os.path.join(directory, 'arrest_coefficient'))


def plot_displacement_hist(profiles, directory):
    displacement = [t.displacement for pro in profiles for t in pro.tracks]
    plot_utils.plot_histogram(data=displacement,
                              filename=os.path.join(directory, 'displacements-hist'),
                              xlabel='Displacement (micrometers)',
                              ylabel='Frequency', title='Track Displacements')


def plot_displacement_ecdf(profiles, directory):
    displacement = [t.displacement for pro in profiles for t in pro.tracks]
    plot_utils.plotCDFs([displacement], xlabel='Displacement (micrometers)', title='Track Displacements',
                        filename=os.path.join(directory, 'displacements-ecdf'))


def plot_track_displacement_timeseries(profiles, directory=None, graph_path=None, legend=False, 
                                       xlabel='Absolute time (min)', ylabel='Displacement (um)',
                                       xlim=None, ylim=None):
    """ Plots each track's displacements as a separate time series on the same graph. """
    plt.clf()
    # Unique colour for each track.
    num_tracks = len([t for p in profiles for t in p.tracks])
    colours = unique_colour_per_series(num_tracks)
    ti = 0  # Index of track, use to index colours. 
    for p in profiles:
        for t in p.tracks:
            track_displacements = [pos.displacement for pos in t.positions]
            track_timestamps_min = [pos.time_m for pos in t.positions]
            plt.plot(track_timestamps_min, track_displacements, color=colours[ti], label=t.ID)
            ti += 1
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    if legend:
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, fontsize='xx-small',
                   ncol=5, mode="expand", borderaxespad=0.)    
    if not graph_path:
        graph_path = os.path.join(directory, 'track_displacement_timeseries')
    plt.savefig(graph_path + '.png', dpi=600)


def plot_meandering_hist(profiles, directory):
    # all tracks' meandering indexes. Note that broken tracks don't have meandering indexes assigned, as true
    # length of track is not known.
    meander = [t.meander for pro in profiles for t in pro.tracks if t.meander]
    plot_utils.plot_histogram(data=meander,
                              filename=os.path.join(directory, 'meanderingIndexes-hist'),
                              xlabel='Meandering index',
                              ylabel='Frequency', xLim=(0.0, 1.0), title='Track Meandering Indexes')


def plot_meandering_ecdf(profiles, directory):
    # all tracks' meandering indexes. Note that broken tracks don't have meandering indexes assigned, as true
    # length of track is not known.
    meander = [t.meander for pro in profiles for t in pro.tracks if t.meander]
    plot_utils.plotCDFs([meander], xlabel='Meandering index', title='Track Meandering Indexes',
                        filename=os.path.join(directory, 'meanderingIndexes-ecdf'))


def plot_median_speed_hist(profiles, directory):
    medianSpeed = []
    for p in profiles:
        medianSpeed.extend([t.median_speed() for t in p.tracks])
    plot_utils.plot_histogram(data=medianSpeed,
                              filename=os.path.join(directory, 'speedMedian-hist'),
                              xlabel='Median cell speed (um/min)',
                              ylabel='Frequency', title='Track Median Speeds')


def plot_median_speed_ecdf(profiles, directory):
    medianSpeed = []
    for p in profiles:
        medianSpeed.extend([t.median_speed() for t in p.tracks if not math.isnan(t.median_speed())])
    plot_utils.plotCDFs([medianSpeed], xlabel='Median cell speed (um/min', title='Median Track Speeds',
                        filename=os.path.join(directory, 'speedMedian-ecdf'))


def plot_median_turn_hist(profiles, directory):
    medianTurn = [t.median_turn() for p in profiles for t in p.tracks if not math.isnan(t.median_turn())]
    plot_utils.plot_histogram(data=medianTurn,
                              filename=os.path.join(directory, 'turnMedian-hist'),
                              xlabel='Turn velocity (degrees/min)',
                              ylabel='Frequency', title='Track Median Turns')


def plot_median_turn_ecdf(profiles, directory):
    median_turn = [t.median_turn() for p in profiles for t in p.tracks if not math.isnan(t.median_turn())]
    plot_utils.plotCDFs([median_turn], xlabel='Turn velocity (degrees/min)', title='Track Median Turn Speeds',
                        filename=os.path.join(directory, 'turnMedian-ecdf'))


def plot_median_roll_hist(profiles, directory):
    medianRoll = [t.median_roll() for p in profiles for t in p.tracks if not math.isnan(t.median_roll())]
    plot_utils.plot_histogram(data=medianRoll,
                              filename=os.path.join(directory, 'rollMedian-hist'),
                              xlabel='Roll velocity (degrees/min)',
                              ylabel='Frequency', title='Track Median Rolls')


def plot_irq_speeds_hist(profiles, directory):
    irqSpeeds = []
    for p in profiles:
        irqSpeeds.extend([t.irq_speed() for t in p.tracks if not math.isnan(t.irq_speed())])
    plot_utils.plot_histogram(data=irqSpeeds,
                              filename=os.path.join(directory, 'speedIRQ-hist'),
                              xlabel='IRQ of cell speeds (um/min)',
                              ylabel='Frequency',
                              title='IRQ of Track Speeds')


def plot_irq_speeds_ecdf(profiles, directory):
    irqSpeeds = []
    for p in profiles:
        irqSpeeds.extend([t.irq_speed() for t in p.tracks if not math.isnan(t.irq_speed())])
    plot_utils.plotCDFs([irqSpeeds], xlabel='IRQ of cell speeds (um/min)', title='IRQ of Track Speeds',
                        filename=os.path.join(directory, 'speedIRQ-ecdf'))


def plot_turns_hist(profiles, directory):
    turns = []
    for p in profiles:
        turns.extend(p.collate_turns())
    plot_utils.plot_histogram(data=turns,
                              filename=os.path.join(directory, 'turns'),
                              xlabel='Turn velocity (degrees/min)',
                              ylabel='Frequency', title='Cell Turn Velocities')


def plot_turns_ecdf(profiles, directory, xmax=None):
    turns = []
    for p in profiles:
        turns.extend(p.collate_turns())
    plot_utils.plotCDFs([turns], xlabel='Turn velocity (degrees/min)', title='Cell Turn Velocities',
                        filename=os.path.join(directory, 'turns-ecdf'), xmax=xmax)


def plot_rolls_hist(profiles, directory):
    # not all positions have roll data, if the track was broken.
    rolls = []
    for p in profiles:
        rolls.extend(p.collate_rolls())
    plot_utils.plot_histogram(data=rolls,
                              filename=os.path.join(directory, 'rolls'),
                              xlabel='Roll velocity (degrees/min)',
                              ylabel='Frequency', title='Cell Roll Velocities')


def plot_rolls_ecdf(profiles, directory):
    # not all positions have roll data, if the track was broken.
    rolls = []
    for p in profiles:
        rolls.extend(p.collate_rolls())
    plot_utils.plotCDFs([rolls], xlabel='Roll velocity (degrees/min)', title='Cell Roll Velocities',
                        filename=os.path.join(directory, 'rolls-ecdf'))


def plot_medSpd_irqSpd_scat(profiles, directory):
    medSpd = [t.median_speed() for p in profiles for t in p.tracks]
    irqSpd = [t.irq_speed() for p in profiles for t in p.tracks]
    rho, pval = scipy.stats.spearmanr(medSpd, irqSpd)
    plt.clf()
    plt.scatter(medSpd, irqSpd)
    plt.xlabel('Median speed (um/min)')
    plt.ylabel('Speed IRQ (um/min)')
    plt.title('rho = {:5.2f}; pval = {:5.2f}'.format(rho, pval))
    plt.savefig(os.path.join(directory, 'medSpd_vs_irqSpd.png'), dpi=600)
    plt.savefig(os.path.join(directory, 'medSpd_vs_irqSpd.eps'), dpi=600)


def plot_medSpd_medTurn_scat(profiles, directory):
    medSpd = [t.median_speed() for p in profiles for t in p.tracks]
    medTurn = [t.median_turn() for p in profiles for t in p.tracks]
    rho, pval = scipy.stats.spearmanr(medSpd, medTurn)
    plt.clf()
    plt.scatter(medSpd, medTurn)
    plt.xlabel('Median speed (um/min)')
    plt.ylabel('Median turn (degrees)')
    plt.title('rho = {:5.2f}; pval = {:5.2f}'.format(rho, pval))
    plt.savefig(os.path.join(directory, 'medSpd_vs_medTurn.png'), dpi=600)
    plt.savefig(os.path.join(directory, 'medSpd_vs_medTurn.eps'), dpi=600)


def plot_medSpd_observations(profiles, directory):
    """
    How significant is the effect that fast moving cells leave the viewing area quickly? Does this bias the data
    towards slower moving cells?
    :return:
    """
    med_spd = []   # the track's median speed
    observations = []   # number of recorded positions for the track
    for p in profiles:
        for t in p.tracks:
            med_spd.append(t.median_speed())
            observations.append(len(t.positions))
    rho, pval = scipy.stats.spearmanr(med_spd, observations)
    plt.clf()
    plt.scatter(med_spd, observations)
    plt.xlabel('Median track translational speed (um/min)')
    plt.ylabel('Number of track observations')
    plt.title('rho = {:5.2f}; pval = {:5.4f}'.format(rho, pval))
    font = {'size': 18}
    plt.rc('font', **font)
    plt.gca().grid(True)   # turn on grid lines.
    # change the width of the plot boundary lines
    ax = plt.gca()
    [i.set_linewidth(2.0) for i in ax.spines.values()]
    plt.savefig(os.path.join(directory, 'medSpd_vs_observations.png'), dpi=300)
    plt.savefig(os.path.join(directory, 'medSpd_vs_observations.eps'), dpi=300)
    return rho, pval


def plot_medTurn_observations(profiles, directory):
    """
    How significant is the effect of directional cells leave the viewing area quickly? Does this bias the data
    towards less directional cells?
    :return:
    """
    med_turn = []   # the track's median speed
    observations = []   # number of recorded positions for the track
    for p in profiles:
        for t in p.tracks:
            med_turn.append(t.median_turn())
            observations.append(len(t.positions))
    rho, pval = scipy.stats.spearmanr(med_turn, observations)
    plt.clf()
    plt.scatter(med_turn, observations)
    plt.xlabel('Median track turn speed (degrees/min)')
    plt.ylabel('Number of track observations')
    plt.title('rho = {:5.2f}; pval = {:5.4f}'.format(rho, pval))
    font = {'size': 18}
    plt.rc('font', **font)
    plt.gca().grid(True)   # turn on grid lines.
    # change the width of the plot boundary lines
    ax = plt.gca()
    [i.set_linewidth(2.0) for i in ax.spines.values()]
    plt.savefig(os.path.join(directory, 'medTurn_vs_observations.png'), dpi=300)
    plt.savefig(os.path.join(directory, 'medTurn_vs_observations.eps'), dpi=300)
    return rho, pval


def plot_mi_length_scat(profiles, directory):
    mi = []
    l = []
    for p in profiles:
        for t in p.tracks:
            if t.meander and t.length:
                mi.append(t.meander)
                l.append(t.length)
    rho, pval = scipy.stats.spearmanr(mi, l)
    plt.clf()
    plt.scatter(l, mi)
    plt.xlabel('Length (um)')
    plt.ylabel('Meandering index')
    plt.title('rho = {:5.2f}; pval = {:5.2f}'.format(rho, pval))
    plt.savefig(os.path.join(directory, 'MI_vs_length.png'), dpi=600)
    plt.savefig(os.path.join(directory, 'MI_vs_length.eps'), dpi=600)


def plot_mi_displacement_scat(profiles, directory):
    mi = []
    d = []
    for p in profiles:
        for t in p.tracks:
            if t.meander and t.displacement:
                mi.append(t.meander)
                d.append(t.displacement)
    rho, pval = scipy.stats.spearmanr(mi, d)
    plt.clf()
    plt.scatter(d, mi)
    plt.xlabel('Displacement (um)')
    plt.ylabel('Meandering index')
    plt.title('rho = {:5.2f}; pval = {:5.2f}'.format(rho, pval))
    plt.savefig(os.path.join(directory, 'MI_vs_displacement.png'), dpi=600)
    plt.savefig(os.path.join(directory, 'MI_vs_displacement.eps'), dpi=600)


def plot_length_displacement_scat(profiles, filter_min=False, filter_max=False, graph=False, directory=None, name=None):
    l = []
    d = []
    for p in profiles:
        for t in p.tracks:
            if t.length and t.displacement:   # check that the track has this data (broken tracks may not).
                # if the filter is false, proceed. If filter is true (or has a non zero value), then proceed only
                # if the track length exceeds filter value.
                if not filter_min or (filter_min and t.length >= filter_min) and \
                   not filter_max or (filter_max and t.length <= filter_max):
                    l.append(t.length)
                    d.append(t.displacement)
    if l and d:     # ensure there is some data to actually graph.
        rho, pval = scipy.stats.spearmanr(l, d)
        if graph:
            plt.clf()
            plt.scatter(l, d)
            plt.xlabel('Length (um)')
            plt.ylabel('Displacement (um)')
            plt.title('rho = {:5.2f}; pval = {:5.2f}'.format(rho, pval))
            if name is None:
                name = 'length_vs_displacement'
            plt.savefig(os.path.join(directory, name + '.png'), dpi=600)
            plt.savefig(os.path.join(directory, name + '.eps'), dpi=600)
        return {'rho': rho, 'p': pval}
    else:
        return {'rho': float('nan'), 'p': float('nan')}


def plot_turn_vs_speed(profiles, directory=None, xlim=None, ylim=None):
    """
    Plots a density scatter plot of translational speeds against turn speeds. This is done for all recorded
    movements of all cells in the supplied profile(s).

    :param profiles: list of Profile objects
    :param directory: directory on file system to which data will be written
    :param xlim: tuple or list of two items, [min, max], the ranges to show for the xaxis
    :param ylim: tuple or list of two items, [min, max], the ranges to show for the yaxis
    :return:
    """
    data = []   # tuple (speed, turn)
    for p in profiles:
        for t in p.tracks:
            for j in range(1, len(t.positions)):
                # first item has no turn angle. it is the angle from prev vector to this one.
                pi = t.positions[j-1]
                pj = t.positions[j]
                if pj.speed and pi.speed and pj.turn:
                    spd_i = pi.speed
                    spd_j = pj.speed
                    spd = (spd_i + spd_j) / 2.0
                    turn = pj.turn
                    data.append((spd, turn))   # store as a tuple.
    spds = [t[0] for t in data]
    turns = [t[1] for t in data]
    rho, p = scipy.stats.spearmanr(spds, turns)
    plt.clf()
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    # going to make a density scatter plot. This calculates the density using KDE.
    xy = np.vstack([spds, turns])
    z = scipy.stats.gaussian_kde(xy)(xy)
    plt.scatter(spds, turns, c=z, s=100, edgecolor='')
    plt.xlabel('Translational speeds (um/min)')
    plt.ylabel('Turn speeds (degrees/min)')
    plt.title('rho = {:5.2f}; pval = {:.5f}'.format(rho, p))
    plt.gca().grid(True)   # turn on grid lines.
    font = {'size': 18}
    plt.rc('font', **font)
    plt.savefig(os.path.join(directory, 'speeds_vs_turns.png'), dpi=300)
    plt.savefig(os.path.join(directory, 'speeds_vs_turns.eps'), dpi=300)


def plot_turns_vs_timestep(profiles, directory=None):
    data = []  # tuple, (turn speed, timestep)
    for p in profiles:
        for t in p.tracks:
            for j in t.positions[1:]:   # first item has no turn information.
                if j.turn is not None and not np.isnan(j.turn):
                    data.append((j.turn, t.timestep_min))
    turns = [t[0] for t in data]
    timesteps = [t[1] for t in data]
    if np.min(timesteps) == np.max(timesteps):
        print('Cannot process turn speed vs dt correlations: all data have the same dt.')
        return
    rho, p = scipy.stats.spearmanr(turns, timesteps)
    print("turn speed and dt correlations = {r:.2f}, p value = {p:.5f}".format(r=rho, p=p))
    plt.clf()
    # going to make a density scatter plot. This calculates the density using KDE.
    xy = np.vstack([turns,timesteps])
    z = scipy.stats.gaussian_kde(xy)(xy)
    plt.scatter(turns, timesteps, c=z, s=100, edgecolor='')
    plt.xlabel('Turn speeds (degrees/min)')
    plt.ylabel('Time step (min)')
    plt.title('rho = {:5.2f}; pval = {:5.2f}'.format(rho, p))
    plt.gca().grid(True)   # turn on grid lines.
    font = {'size': 18}
    plt.rc('font', **font)
    plt.savefig(os.path.join(directory, 'turns_vs_timestep.png'), dpi=300)
    plt.savefig(os.path.join(directory, 'turns_vs_timestep.eps'), dpi=300)


def plot_msd_errorbars(profiles, labels=[], directory=None, profile_colours=None, max_dt_prop=0.2, 
                       xlim=None, ylim=None, grid=True):
    """
    If profiles is a list (of Profiles), each will be drawn as a separate series.
    Standard Error mean is indicated by filled areas on the plot. 
    The legend indicates slopes of linear regression of msd vs time interval in log-log space. 
    max_dt_prop is the maximum time interval over which MSD is calculated, specified as a proportion [0, 1] of the 
    maximally observed time span of any track in the given profile (ie, max msd time interval is profile-specific. )

    """
    if isinstance(profiles, profile.Profile):
        profiles = [profiles]  # Set up correct formatting
    
    if not profile_colours:
        # +1 at end to avoid an extreme yellow colour if there are only two. 
        profile_colours = plt.get_cmap('plasma')(np.linspace(0, 1.0, len(profiles) + 1))
    
    plt.clf()
    # Plot data for each profile
    for prof_index, p in enumerate(profiles):
        msds_all = dict()
        for t in p.tracks:
            deltaT_displacements_sq = t.get_deltaT_displacements_sq()
            for dt in deltaT_displacements_sq.keys():
                if dt not in msds_all.keys():
                    msds_all[dt] = []
                msds_all[dt].extend(deltaT_displacements_sq[dt])
        # Calculate mean values only for delta t displacements within the specified range.
        dts = sorted(msds_all.keys())
        cutoff = dts[-1] * max_dt_prop  # This cut off is specific to each dataset, not global. 
        dts = np.array([dt for dt in dts if dt <= cutoff])  # Filter
        dts_min = dts / 60.
        msds = np.array([np.mean(msds_all[dt]) for dt in dts])
        sems = np.array([scipy.stats.sem(msds_all[dt]) for dt in dts])  # Standard error mean. 
        
        # Linear regression on lot transformed data.
        dts_log = np.log(dts)
        msds_log = np.log(msds)
        slope, intercept, r, p, stderr = scipy.stats.linregress(x=dts_log, y=msds_log)

        lab = labels[prof_index] if labels else ''  
        lab += ' ({:.2f})'.format(slope)  # Add linear regression slope
    
        plt.fill_between(dts_min, msds-sems, msds+sems, color='lightgrey')
        plt.plot(dts_min, msds, color=profile_colours[prof_index], label=lab)

    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.ylabel("Mean squared displacement (um)")
    plt.xlabel("Time interval (minutes)")
    plt.gca().grid(grid)  # Turn on grid lines.
    if labels:
        plt.legend(loc='lower right')    
    plt.savefig(directory + '/msd_error_linear.eps', dpi=600)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(directory + '/msd_error_log.eps', dpi=600)


def plot_msd_allT(profiles, directory=None, separateProfiles=True, xlim=None, ylim=None):
    """
    :param profiles:
    :param directory:
    :param separateProfiles:
    :param xlim: tuple or list of two items, [min, max], the ranges to show for the xaxis
    :param ylim: tuple or list of two items, [min, max], the ranges to show for the yaxis
    :return:
    """
    msds_all = dict()
    # Collate all data together.
    for p in profiles:
        for t in p.tracks:
            deltaT_displacements_sq = t.get_deltaT_displacements_sq()
            for dt in deltaT_displacements_sq.keys():
                if dt not in msds_all.keys():
                    msds_all[dt] = []
                msds_all[dt].extend(deltaT_displacements_sq[dt])
    # Calculate mean values only for delta t displacements within the specified range.
    dts = sorted(msds_all.keys())
    cutoff = dts[-1] * 0.25
    dts = [dt for dt in dts if dt <= cutoff]
    msds = [np.mean(msds_all[dt]) for dt in dts]

    # Calculate slope.
    # Log-transform the data to find the slope.
    dts_log = [math.log(x) for x in dts]
    msds_log = [math.log(x) for x in msds]
    # Linear regression on lot transformed data.
    slope, intercept, r, p, stderr = scipy.stats.linregress(x=dts_log, y=msds_log)
    # Data to plot regression line on the plot. Plot will be log-transformed, so need to raise all values here to
    # e, otherwise it ends up double-log-transformed.
    linear = [math.exp((x * slope) + intercept) for x in dts_log]
    # plot.
    plt.clf()
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    if separateProfiles:
        for p in profiles:
            p_msds = {}
            for t in p.tracks:
                deltaT_displacements = t.get_deltaT_displacements_sq()
                for dt in deltaT_displacements_sq.keys():
                    if dt not in p_msds.keys():
                        p_msds[dt] = []
                    p_msds[dt].extend(deltaT_displacements_sq[dt])
            times = sorted(p_msds.keys())
            if len(times) > 0:
                cutoff = times[-1] * 0.25
                times = [t for t in times if t <= cutoff]
                vals = [np.mean(p_msds[dt]) for dt in times]
                plt.plot(times, vals, 'grey')
    else:
        # plots the empirical data of all profiles collated together. This only works if all profiles employed a
        # similar timestep value, otherwise you're best off plotting each one individually
        # (use separateProfiles=True)
        plt.plot(dts, msds)
    plt.plot(dts, linear, 'r', linewidth=2.0)
    plt.xlabel('Time (s)')
    plt.ylabel('MSD')
    plt.title('Slope = {s:.2f}'.format(s=slope))
    plt.xscale('log')
    plt.yscale('log')
    plt.gca().grid(True)   # turn on grid lines.
    ax = plt.gca()
    ax.set_xlim([dts[0]*0.9, ax.get_xlim()[1]])
    plt.savefig(os.path.join(directory, 'msd-allT.png'), dpi=600)
    plt.savefig(os.path.join(directory, 'msd-allT.eps'), dpi=600)


def plot_msd_vs_max_dt(profiles, directory=None, upper_max_dt=float('inf')):
    """ Plots how MSD slope differs with the maximum time interval used in regression. """
    max_dts, slopes = profile.Profile.calculate_msd_vs_max_dt(profiles, upper_max_dt)
    plt.clf()
    plt.plot(max_dts, slopes)
    plt.xlabel('Maximum time interval (s)')
    plt.ylabel('MSD slope')
    plt.gca().grid(True)  # turn on grid lines.
    plt.savefig(os.path.join(directory, 'msd_vs_max_dt.png'), transparent=True, dpi=600)


def plot_track_msd(profiles, directory, max_dt=None, max_dt_prop=0.25, xlim=None, ylim=None):
    """
    :param profiles:
    :param directory:
    :param xlim: tuple or list of two items, [min, max], the ranges to show for the xaxis
    :param ylim: tuple or list of two items, [min, max], the ranges to show for the yaxis
    :return:
    """
    msds_all = dict()
    # List of dictionaries, representing MSD stats (dictionary) for each track (list). 
    # Each dictionary (MSD) has as key a dt value which points to a list of MDS values for that given dt. 
    # list( dict [dt] -> list [msd values] )
    track_msds = list()  
    # Collate all data together.
    max_track_duration = np.max([t.duration_sec for p in profiles for t in p.tracks])
    plt.clf()
    for p in profiles:        
        for t in p.tracks:            
            track_msd = []  # MSD series for the given track. 
            deltaT_displacements = t.get_deltaT_displacements_sq()
            track_dts = sorted(deltaT_displacements_sq.keys())
            if max_dt is None:  # Left this code here in case we return to calculating this for each track individually.
                max_dt = max_track_duration * max_dt_prop
                                
            track_dts = [dt for dt in track_dts if dt < max_dt]
            for dt in track_dts:
                # Mean squared displacement value for this given dt. 
                msd_dt = np.mean(deltaT_displacements_sq[dt])
                track_msd.append(msd_dt)
            plt.plot(track_dts, track_msd)
                   
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.xlabel('Time (s)')
    plt.ylabel('MSD')
    # plt.title('Slope = {s:.2f}'.format(s=slope))
    plt.xscale('log')
    plt.yscale('log')
    plt.gca().grid(True)   # turn on grid lines.
    # ax = plt.gca()
    # ax.set_xlim([dts[0]*0.9, ax.get_xlim()[1]])
    plt.savefig(os.path.join(directory, 'msd-tracks.png'), dpi=600)
    plt.savefig(os.path.join(directory, 'msd-tracks.eps'), dpi=600)


def calculate_instantaneous_msd(track, half_window_width):
    """
    The window is centred on a jump (not a position). 
    Half window width refers to the number of jumps either side of AND INCLUDING the window centre. 
    If you have 14 positions, then there are 13 jumps. 
    
    With a half window width of 4, the first window centre covers to it's left positions 1-4 (starting from 1, e.g. Imaris)
    The first window covers to it's right 6-9 inclusive. 
    Window is centred on 4. 

    With a half window width of four, and 14 positions, 14 - 4 - 4 = 8 positions are lost. 
    That leaves 6 windows, and 6 window centres. 
    Yet there are only 13 jumps. This is why you have 4 dead jumps at the start, and 3 at the end. 

    'x' refers to the index used in slicing. To fully slice [0, 1, 2] you need list[0:3]
    '@' refers to windo centre. 
    14 positions:
    1   2   3   4   5   6   7   8   9  10  11  12  13  14  
    | a . b . c . @ . c . b . a |   x
        | a . b . c . @ . c . b . a |   x
            | a . b . c . @ . c . b . a |   x
                | a . b . c . @ . c . b . a |   x
                    | a . b . c . @ . c . b . a |   x
                        | a . b . c . @ . c . b . a |   x
                            | a . b . c . @ . c . b . a |   x

    7 window centres (carry the central @s down):
    .   .   .   . @ . @ . @ . @ . @ . @ . @ .   .   .   . 

    13 jumps:                                           |
    | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11| 12| 13|



    Trust me on this, I've literally spent hours of debugging and eventually resorting to pen and paper to figure this 
    out. This method, and the plotting in plot_track_intantaneous_msd below work as they should. 


    See Banigan 15, PLoS Computational Biology.
    Captures the distribution of displacements for a given time interval (dt).
    This is performed for subsets of positional data, derived from a sliding window applied over all the positional data.
    Both displacements along each component spatial axis {x, y, z} and net displacement are calculated.

    :return: dict[window centre] -> dict[dt] -> list[displacements found at that dt within the given window]
    """
    # Both of these data structure work thusly:
    # dict[window centre] -> dict[dt] -> list[displacements found at that dt within the given window]
    window_data_dt_disp = dict()  # Returnable from this function. 
    window_data_dt_disp_sq = dict()

    observations = len(track.positions)
    for window_centre in range(half_window_width-1, observations - half_window_width):
        window_start = window_centre - (half_window_width-1)
        # Used in slicing. Hence, not the index of the last item, but one beyond it. 
        window_end = window_centre + half_window_width  + 1
        subset_positions = track.positions[window_start : window_end]

        # All of the form dict[dt] -> list[displacements observed for that dt]
        deltaT_displacements = dict()
        deltaT_displacements_sq = dict()
        deltaT_displacements_X = dict()  # Currently unused, but could be in future. 
        deltaT_displacements_Y = dict()
        deltaT_displacements_Z = dict()

        # Iterate over all possible delta Ts, in terms of recorded data (not actual time).
        for d in range(1, len(subset_positions) - 1):
            delta_t = d * track.timestep_s

            # Process all possible temporal locations where interval delta_t occurred.
            for i in range(0, len(subset_positions) - d):
                # Can't be done for positions where cell was untracked.
                if subset_positions[i+d].tracked() and subset_positions[i].tracked():
                    dx = subset_positions[i+d].x - subset_positions[i].x
                    dy = subset_positions[i+d].y - subset_positions[i].y
                    dz = subset_positions[i+d].z - subset_positions[i].z
                    vector = (dx, dy, dz)
                    disp = geometry.vector_length(vector)
                    disp_sq = disp * disp

                    # Initialise if not already done so for this delta_t. 
                    # By doing this heare we only include a record for delta_t if there is associated data.
                    if delta_t not in deltaT_displacements.keys():
                        deltaT_displacements[delta_t] = []
                        deltaT_displacements_sq[delta_t] = []
                        deltaT_displacements_X[delta_t] = []  # Displacements in each component axis
                        deltaT_displacements_Y[delta_t] = []
                        deltaT_displacements_Z[delta_t] = []

                    deltaT_displacements[delta_t].append(disp)
                    deltaT_displacements_sq[delta_t].append(disp_sq)
                    # Component displacements along each axis.
                    deltaT_displacements_X[delta_t].append(dx)
                    deltaT_displacements_Y[delta_t].append(dy)
                    deltaT_displacements_Z[delta_t].append(dz)

        window_data_dt_disp[window_centre] = deltaT_displacements
        window_data_dt_disp_sq[window_centre] = deltaT_displacements_sq
    return window_data_dt_disp, window_data_dt_disp_sq


def plot_track_intantaneous_msd(profiles, half_window_width=5,  # All possible dts in the window used.
                                directory=None, graph_path=None,  # Supply either directory or graph_path
                                xlim=None, ylim=None, smooth=False,
                                msd_thresholds=[-1000, 0.8, 1.2],
                                heatmap_abs_time=True,  # If False, xlabel=time step.                                  
                                consequtive_top_threshold_length_upgrade=4,
                                plot_instantaneous_msd_timeseries=True,
                                plot_heatmaps=True,  # Jumps vs instantaneous timesteps for every track. 
                                manual_inspect_tracks=False,  # Plot 3D interactive trajectory for each track in turn.
                                # Not worth doing if there are more than 40 odd tracks. Can be time consuming. 
                                plot_tracks_2D=True  
                                ):
    """
    Plots a graph of the MSD slope for each track, at each point of time. 
    This is done by applying a sliding window over the track's positional data, calculating the MSD vs dt values 
    possible within that subset of positional data, and deriving the slope from linear regression on this data. 
    Multiple tracks are plotted on the same graph. 

    The data are plotted as a heatmap, with absolute experimental time on the x axis, and tracks stacked on the y axis.
    Two heatmaps are made, one of the MSD values themselves, and the other thresholded at value given in msd_thresholds.
    Finally, CSV files detailing the proportion of observations each track spends between each threshold are drawn.
    A similar CSV is created for all track observations combined into a single analysis. 

    Heatmaps can only be plotted when all tracks hare the same time step. 

    :param profiles:
    :param directory:
    :param xlim: tuple or list of two items, [min, max], the ranges to show for the xaxis
    :param ylim: tuple or list of two items, [min, max], the ranges to show for the yaxis
    :return:
    """
    if not (directory or graph_path):
        raise Exception("Please supply either a directory (default graph name used) or graph_path (graph name supplied)")

    num_tracks = len([t for p in profiles for t in p.tracks])
    max_time_min = np.max([pos.time_m for p in profiles for t in p.tracks for pos in t.positions])
    max_time_index = np.max([pos.timeID for p in profiles for t in p.tracks for pos in t.positions])
    timesteps_min = [track.timestep_min for p in profiles for track in p.tracks]  # Timestep of each track analysed. 
    timesteps_min = np.unique(timesteps_min)
    identical_timesteps = len(timesteps_min) == 1  # Cannot process heatmaps if different timesteps used. How to plot?
    
    timestep_min = np.min(timesteps_min)  # Smallest time step found
    num_timesteps = math.ceil(max_time_min / timestep_min)  # Largest time by smallest timestep, accommodates everything
    hm_data = np.empty([num_tracks, num_timesteps + 1])
    h, w = hm_data.shape
    hm_data[:] = np.nan  # Assign default value, omits colour for that square in the heatmap images. 
    max_threshold_value = len(msd_thresholds) -1  # Thresholding applies '0', '1', '2', ... to hm data. 

    if plot_instantaneous_msd_timeseries:
        plt.clf()  # For plotting a time series of instantaneous MSD, one series per track.         
        colours = unique_colour_per_series(num_tracks)  # Assign a unique colour to each track.
    ti = 0  # Track index, used to index colours. 
    for p in profiles:                
        for track in p.tracks:
            _, window_data_dt_disp_sq = calculate_instantaneous_msd(track, half_window_width)
            window_centres = sorted(window_data_dt_disp_sq.keys())
            # Convert to absolute time, using time step length and start time of the track (first position).
            window_centre_times = (np.array(window_centres) * track.timestep_s) + track.positions[0].time_s

            window_slopes = []
            for window_centre in window_centres:
                # Dictionary at the window centre is dict[dt] -> list[displacements for given dt]
                window_dts = sorted(window_data_dt_disp_sq[window_centre].keys())                
                window_msds = [np.mean(window_data_dt_disp_sq[window_centre][dt]) for dt in window_dts]
                # Calculate slope; linear regression on log transformed data.
                dts_log = np.log(window_dts)
                msds_log = np.log(window_msds)
                slope, intercept, r, p, stderr = scipy.stats.linregress(x=dts_log, y=msds_log)
                window_slopes.append(slope)
                if np.isnan(slope):  # This shouldn't happen. Print debug data and throw Exception. 
                    print('\n\n')
                    print(window_dts)
                    print(window_msds)
                    print(window_data_dt_disp_sq[window_centre])
                    raise Exception("NAN slope?!?")
            
            if smooth:
                # Smooth the window_slopes data. 
                window_slopes = lowess(endog=np.array(window_slopes), exog=np.array(window_centre_times), frac=(1./3.), 
                                       return_sorted=False)            
            
            window_location_indices = [i + int(track.positions[0].timeID) 
                                       - 1  # Imaris counts time indicies from 1, but np.array starts from 0
                                       for i in window_centres]
            hm_data[ti, window_location_indices] = window_slopes

            if len(window_slopes) >= 2 and plot_instantaneous_msd_timeseries:  # Don't try to plot nothing. 
                # Plot the MSD slope for this track as calculated over the present sliding window location. 
                plt.plot(window_centre_times / 60., # Convert to minutes.
                         window_slopes, label='{:.0f}'.format(round(track.ID)), color=colours[ti])  
            ti += 1

    if plot_instantaneous_msd_timeseries:
        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)
        plt.xlabel('Time (min)')
        plt.ylabel('MSD slope')
        # plt.gca().grid(True)   # Turn on grid lines.
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, fontsize='x-small',
                ncol=5, mode="expand", borderaxespad=0.)
        # ax = plt.gca()
        # ax.set_xlim([0, ax.get_xlim()[1]])
        if not graph_path:
            graph_path = os.path.join(directory, 'msd-tracks_sliding_window')
        plt.savefig(graph_path + '.png', dpi=600)
        plt.savefig(graph_path + '.eps', dpi=600) 
        plt.close() 

    # ------ Perform thresholding on hm_data. 
    track_ids = ['{:.0f}'.format(t.ID) for p in profiles for t in p.tracks]
    thresh_data = np.empty(hm_data.shape) 
    thresh_data[:] = np.nan  # Default values, doesn't draw colour on image. 

    # Capture all data for which instantaneous MSD values can be calculated.     
    for ti, thresh in enumerate(msd_thresholds):
        indices = np.where(thresh <= hm_data) 
        thresh_data[indices[0], indices[1]] = ti

    # Additional top tier of directionality, persisting directional motion, assigned to consequtive timesteps of length
    # greater than user-defined threshold.
    max_threshold_value += 1  # New max threshold value, sustained directionality. 
    for ti in list(range(thresh_data.shape[0])):        
        indices = np.where(thresh_data[ti, :] == max_threshold_value - 1)  # Extract instances of prev max threshold.
        time_step_indices = indices[0]       
        for k, g in groupby(enumerate(time_step_indices), lambda x: x[0]-x[1]):            
            subset_indices = list(map(itemgetter(1), g))
            if len(subset_indices) > consequtive_top_threshold_length_upgrade:
                thresh_data[ti, subset_indices] = max_threshold_value
        
    # Save thresholded MSD values
    thresh_df = pandas.DataFrame(data=thresh_data, index=track_ids, 
                                 # Start from 1. Data point registers between times 0 and 1. Imaris counts from 1 onwards.
                                 columns=['ti_{:d}-{:d}'.format(i+1, i+2)  
                                          for i in list(range(hm_data.shape[1])) ])
    thresh_df.to_csv(directory + '/msd-tracks_instantaneous_threshold_RAW' + str(msd_thresholds) + '.csv')

    # ------ Plot heatmaps, both raw MSD slopes and thresholded. 
    if identical_timesteps and plot_heatmaps:
        # Perform plot as a heatmap of MSD values    
        h, w = hm_data.shape
        # figaspect has GOT to be the most useful and obviously absent function related to matplotlib. 
        # Fit the figure size around the DATA, not squishing the data within some arbitrarily defined figure size. 
        #
        # Difficult to adjust the figure size to accommodate the heatmap plot, which varies with the specific dataset
        # being analysed. Current settings of h ** 2.0 still results in overlapping for 200 tracks or more. 
        # There's something non-liner going on here that I can't untangle. Best advice is, if you're going to 
        # plot more than 200, either split them into several profiles and plot separately, or mess around with these
        # values directly. 
        fig_w, fig_h = plt.figaspect( (h ** 2) / w ) 
        fig = plt.figure(figsize=(fig_w,fig_h))
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        cmap = plt.get_cmap('viridis')
        cmap.set_bad('w',1.)  # Set bad color to white. EPS cannot handle transparency, which is used for bad colours. 
        mesh = ax.pcolormesh(hm_data, cmap=cmap)
        plt.colorbar(mesh)
        # Invert the axis, put origin on top left (not bottom). First track at top. 
        ax.set_ylim(ax.get_ylim()[::-1])  
        if heatmap_abs_time:
            plt.xlabel('Absolute time (min)')
            locs, labels = plt.xticks()  # Get locations and labels
            time_stamps = np.array(locs * timestep_min)
            time_stamps = np.round(time_stamps, 1)  # Round to 1 decimal place
            plt.xticks(locs, time_stamps, rotation=45)
            plt.xlim([0, max_time_index])  # Constrain to sensible range, matplotlib can lose the plot here. 
        else:
            plt.xlabel('Time step')
        plt.yticks(np.arange(len(track_ids)) + 0.5, track_ids)
        plt.ylabel('Track')
        plt.savefig(directory + '/msd-tracks_instantaneous.png', dpi=600, bbox_inches='tight')
        plt.savefig(directory + '/msd-tracks_instantaneous.eps', bbox_inches='tight')
        plt.close()

        # Plot heatmap based on thresholded values
        fig = plt.figure(figsize=(fig_w,fig_h))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])        
        cmap = plt.get_cmap('plasma')
        cmap.set_bad('w',1.)  # Set bad color to white. EPS cannot handle transparency, which is used for bad colours. 
        mesh = ax.pcolormesh(thresh_data, cmap=cmap, vmax=max_threshold_value)
        # Invert the axis, put origin on top left (not bottom). First track at top. 
        ax.set_ylim(ax.get_ylim()[::-1]) 
        if heatmap_abs_time:
            plt.xlabel('Absolute time (min)')
            locs, labels = plt.xticks()  # Get locations and labels
            time_stamps = np.array(locs * timestep_min)
            time_stamps = np.round(time_stamps, 1)  # Round to 1 decimal place
            plt.xticks(locs, time_stamps, rotation=45)
            plt.xlim([0, max_time_index])  # Constrain to sensible range, matplotlib can lose the plot here. 
        else:
            plt.xlabel('Time step')
        plt.yticks(np.arange(len(track_ids)) + 0.5, track_ids)
        plt.ylabel('Track')
        plt.savefig(directory + '/msd-tracks_instantaneous_threshold' + str(msd_thresholds) + '.eps', 
                    bbox_inches='tight')
        plt.savefig(directory + '/msd-tracks_instantaneous_threshold' + str(msd_thresholds) + '.png', 
                    dpi=600, bbox_inches='tight')       
        plt.close()

    # ------
    # Make CSV of proportions and absolute numbers of observations each track spends at each instantaneous MSD threshold.
    rows = []
    # Set up column names, in part for maintaining desired order. 
    column_names = ['track ID', 'total observations']
    abs_names = []
    prop_names = []
    str_thresholds = [str(i) for i in msd_thresholds]
    str_thresholds.append('directed motion')
    for thresh in str_thresholds:
        abs_names.append('absolute observations > ' + str(thresh) + ' and < next column')
        prop_names.append('proportion observations > ' + str(thresh) + ' and < next column')
    column_names.extend(abs_names)
    column_names.extend(prop_names)

    for track_i in list(range(thresh_data.shape[0])):  # Cycle through each track. 
        track_data = thresh_data[track_i, :]  # Extract this specific track from 2D matrix. 
        total_observations = np.count_nonzero(~np.isnan(track_data))
        if total_observations > 0:  # Avoid div by zero. This track won't be included in the csv file. 
            row = {
                'track ID': track_ids[track_i],
                'total observations': total_observations,
            } 

            for ti, thresh in enumerate(str_thresholds):
                name = 'absolute observations > ' + str(thresh) + ' and < next column'
                abs_count = np.count_nonzero(track_data == ti) 
                row[name] = abs_count

                name = 'proportion observations > ' + str(thresh) + ' and < next column'
                prop = abs_count / total_observations
                row[name] = prop

            rows.append(row)
    df = pandas.DataFrame(rows, columns=column_names)
    df.to_csv(directory + '/msd-tracks_instantaneous_threshold_proportions.csv')

    # Same again, but for all tracks combined into a single analysis. Operates over all observations.
    # Use same column names. 
    total_observations = np.count_nonzero(~np.isnan(thresh_data))
    if total_observations > 0:
        row = {
            'track ID': 'All tracks combined',
            'total observations': total_observations
        }
        for ti, thresh in enumerate(str_thresholds):
            name = 'absolute observations > ' + str(thresh) + ' and < next column'
            abs_count = np.count_nonzero(thresh_data == ti)
            row[name] = abs_count

            name = 'proportion observations > ' + str(thresh) + ' and < next column'
            prop = abs_count / total_observations
            row[name] = prop
        df = pandas.DataFrame([row], columns=column_names)
        df.to_csv(directory + '/msd-all_data_instantaneous_threshold_proportions.csv')
    else:
        print("WARNING! No instantaneous MSD values could be recorded. Half window size too large?")

    if plot_tracks_2D or manual_inspect_tracks:
        # ------ Plot tracks in 2D and (optionally, for manual inspection) 3D
        class TrackData:
            """ Short convenience class to collect together these data. Calculated below for 3D and then used again for
            for plotting 2D. """
            def __init__(self, id):
                self.id = id
                self.xs = []
                self.ys = []
                self.zs = []
                self.cols = []  # Colours to plot each segment. 
                self.untracked_xs = []  # Used to plot large red circles showing where untracked positions are located. 
                self.untracked_ys = []
                self.untracked_zs = []
        
        threshold_colours = plt.get_cmap('plasma')(np.linspace(0, 1.0, max_threshold_value + 1))
        all_tracks_data = []
        for p in profiles:                
            for ti, track in enumerate(p.tracks):
                track_data = TrackData(track.ID)
                for pi, position in enumerate(track.positions):  # First and last positions MUST be tracked. 
                    if position.tracked():
                        track_data.xs.append(position.x)
                        track_data.ys.append(position.y)
                        track_data.zs.append(position.z)
                    else:  # If not tracked, carry last tracked item forward. 
                        track_data.xs.append(track_data.xs[-1])
                        track_data.ys.append(track_data.ys[-1])
                        track_data.zs.append(track_data.zs[-1])
                        track_data.untracked_xs.append(track_data.xs[-1])
                        track_data.untracked_ys.append(track_data.ys[-1])
                        track_data.untracked_zs.append(track_data.zs[-1])

                    # pi starts from zero, but the track may start anywhere in the time series
                    shifted_pi = int(pi + track.positions[0].timeID -1)  # -1 because imaris counts from 1, not 0  
                    threshold_value = thresh_data[ti, shifted_pi]
                    if not np.isnan(threshold_value):
                        c = threshold_colours[int(threshold_value)]
                        track_data.cols.append(c)
                    else:
                        track_data.cols.append('grey')
                all_tracks_data.append(track_data)  # Store for use below in plotting 2D trajectories.

                # Interactive 3D plots of track trajectory. One graph per track.     
                if manual_inspect_tracks:
                    fig_3d = plt.figure()
                    ax = fig_3d.gca(projection='3d')                    
                    for i in list(range(len(track_data.xs))):
                        col = track_data.cols[i] if i != 0 else 'red'  # Indicate the start of the track with a red segment. 
                        # l[1:1+1] gives one item. We need two to plot a line. 
                        ax.plot(track_data.xs[i:i+2], track_data.ys[i:i+2], track_data.zs[i:i+2], 
                                c=col)
                    # Indicate untracked positions with larger red dots. 
                    for x, y, z in zip(track_data.untracked_xs, track_data.untracked_ys, track_data.untracked_zs):
                        ax.scatter(x, y, z, c='red')
                    ax.set_title(track_data.id)
                    ax.set_xlabel('X position')
                    ax.set_ylabel('Y position')
                    ax.set_zlabel('Z position')
                    plt.show()
                    plt.close()

    if plot_tracks_2D:
        # Plot and save tracks in 2D, with segments coloured by MSD threshold value
        # Center the tracks on the origin
        plt.clf()
        linewidth = 2.0  
        fontsize = 18
        for track in all_tracks_data:
            start_x = track.xs[0]
            start_y = track.ys[0]
            start_z = track.zs[0]  # Not actually needed on a 2D plot, but included here for completeness. 
            centred_xs = [x - start_x for x in track.xs]
            centred_ys = [y - start_y for y in track.ys]
            centred_zs = [z - start_z for z in track.zs]
            for i in list(range(len(centred_xs))):
                # l[1:1+1] gives one item. We need two to plot a line, comprising two items. . 
                plt.plot(centred_xs[i:i+2], centred_ys[i:i+2], color=track.cols[i], linewidth=linewidth, label=track.id)
        # Mark the origin. 
        plt.plot([0], [0], marker="8", markersize=6, color="black")
        plt.xlabel('X displacement (um)')
        plt.xticks(rotation=45)
        plt.ylabel('Y displacement (um)')
        plt.gca().grid(True, linewidth=linewidth)   # Turn on grid lines.
        plt.gca().set_aspect('equal', adjustable='box')
        # font = {'size': fontsize}
        # plt.rc('font', **font)
        # Change the width of the plot boundary lines
        ax = plt.gca()
        [i.set_linewidth(linewidth) for i in ax.spines.values()]
        graph_path = os.path.join(directory, 'trajectories_instantaneous_msd_centred')
        plt.savefig(graph_path + '.png', dpi=600)
        plt.savefig(graph_path + '.eps', dpi=600)  
        plt.close()

        plt.clf()
        linewidth = 2.0  
        fontsize = 18
        for track in all_tracks_data:
            for i in list(range(len(track.xs))):
                # l[1:1+1] gives one item. We need two to plot a line, comprising two items. . 
                plt.plot(track.xs[i:i+2], track.ys[i:i+2], color=track.cols[i], linewidth=linewidth, label=track.id)
        plt.xlabel('X displacement (um)')
        plt.xticks(rotation=45)
        plt.ylabel('Y displacement (um)')
        plt.gca().grid(True, linewidth=linewidth)   # Turn on grid lines.
        plt.gca().set_aspect('equal', adjustable='box')
        # font = {'size': fontsize}
        # plt.rc('font', **font)
        # Change the width of the plot boundary lines
        ax = plt.gca()
        [i.set_linewidth(linewidth) for i in ax.spines.values()]
        graph_path = os.path.join(directory, 'trajectories_instantaneous_msd')
        plt.savefig(graph_path + '.png', dpi=600)
        plt.savefig(graph_path + '.eps', dpi=600) 
        plt.close() 


def plot_track_timeseries_msd_lowless(profiles, half_window_width,  # All possible dts in the window used. 
                                  directory=None, graph_path=None,  # Supply either directory or graph_path
                                  xlim=None, ylim=None):
    """
    Plots a graph of the MSD slope for each track, at each point of time. 
    This is done by applying a sliding window over the track's positional data, calculating the MSD vs dt values 
    possible within that subset of positional data, and deriving the slope from linear regression on this data. 
    Multiple tracks are plotted on the same graph. 

    :param profiles:
    :param directory:
    :param xlim: tuple or list of two items, [min, max], the ranges to show for the xaxis
    :param ylim: tuple or list of two items, [min, max], the ranges to show for the yaxis
    :return:
    """
    if not (directory or graph_path):
        raise Exception("Please supply either a directory (default graph name used) or graph_path (graph name supplied)")

    # Can be used to constrain the MSD analysis for all tracks to a common maximum dt. 
    # max_track_duration = np.max([t.duration_sec for p in profiles for t in p.tracks])
    plt.clf()
    for p in profiles:        
        for track in p.tracks:
            _, window_data_dt_disp_sq = calculate_instantaneous_msd(track, half_window_width)
            window_centres = sorted(window_data_dt_disp_sq.keys())
            # Convert to absolute time, using time step length and start time of the track (first position).
            window_centre_times = (np.array(window_centres) * track.timestep_s) + track.positions[0].time_s

            window_slopes = []
            for window_centre in window_centres:
                # Dictionary at the window centre is dict[dt] -> list[displacements for given dt]
                window_dts = sorted(window_data_dt_disp_sq[window_centre].keys())
                window_slopes = []  # Set below. 
                # if lowess_regression:
                ly = []  # Inputs supplied to lowess function. 
                lx = []
                for dt in window_dts:
                    msds = window_data_dt_disp_sq[window_centre][dt]
                    ly.extend(msds)
                    lx.extend([dt] * len(msds))  # Register a dt for each msd value

                # Several msd values will exist for each dt. Must supply pairings of all x and y values. 
                # With return_sorted=False a np.array is returned with predicted y values for each x value supplied. 
                # This will include several of the same y values, because the same x (dt) values are suppleid repatedly.
                smoothed_msd_slope = lowess(endog=np.array(ly), exog=np.array(lx), frac=(2./3.), 
                                            delta=0.001 * (np.max(lx) - np.min(lx)), return_sorted=False)
                window_dts, unique_indices = np.unique(ly, return_index=True)
                window_slopes = smoothed_msd_slope[unique_indices]
                print(window_slopes)
                
                # else:
                # window_msds = [np.mean(window_data_dt_disp_sq[window_centre][dt]) for dt in window_dts]
                # # Calculate slope; linear regression on log transformed data.
                # dts_log = np.log(window_dts)
                # msds_log = np.log(window_msds)
                # slope, intercept, r, p, stderr = scipy.stats.linregress(x=dts_log, y=msds_log)
                # window_slopes.append(slope)
            
            # Plot the MSD slope for this track as calculated over the present sliding window location. 
            plt.plot(window_centre_times / 60., window_slopes)  # Convert to minutes.

    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.xlabel('Time (min)')
    plt.ylabel('MSD slope')
    plt.gca().grid(True)   # Turn on grid lines.
    # ax = plt.gca()
    # ax.set_xlim([0, ax.get_xlim()[1]])
    if not graph_path:
        graph_path = os.path.join(directory, 'msd-tracks_sliding_window')
    plt.savefig(graph_path + '.png', dpi=600)
    plt.savefig(graph_path + '.eps', dpi=600)    


def plot_displacement_autocorrelation(profiles, directory='', max_dt=600, min_corr=1e-3):
    """ see Banigan 15, PLOS Computational Biology, for example of how this is graphed. """
    dac = profile.Profile.collate_displacement_autocorrelation(profiles)
    plt.clf()
    linewidth = 2.0
    keys = sorted(dac.keys())  # needed to extract data in order
    keys = [k for k in keys if k <= max_dt]   # only interested in data for time < 10 min
    # data sorted by key
    data = [dac[k] for k in keys]

    upperErr = []   # length of upper error bar from median (not the absolute value itself)
    lowerErr = []   # length of lower error bar from median (not absolute value itself)
    medianData = []   # the median data to be plotted.
    for d in data:
        median = np.median(d)
        err = np.percentile(a=d, q=[25., 75.])   # find interquartile range
        lowerErr.append(median - err[0])   # find distance of these points from median
        upperErr.append(err[1] - median)
        medianData.append(median)
    lowerErr = np.asarray(lowerErr)   # life is easier working with numpy arrays
    upperErr = np.asarray(upperErr)
    # can't plot negative values on log plot. This cuts off error values at a given threshold.
    ylower = np.maximum(min_corr, medianData - lowerErr)
    lowerErr = medianData - ylower

    keys = [k / 60.0 for k in keys]   # convert seconds to minutes.
    plt.errorbar(keys, medianData, yerr=[lowerErr, upperErr], fmt='o', elinewidth=linewidth)
    plt.xlim((0.0, max_dt / 60.))
    plt.ylim((min_corr, 1.0))
    plt.yscale('log')
    plt.xlabel('Time interval (min)')
    plt.ylabel('Displacement autocorrelation')
    plt.gca().grid(True, linewidth=linewidth)   # turn on grid lines.
    font = {'size': 18}
    plt.rc('font', **font)
    ax = plt.gca()
    # change the width of the plot boundary lines
    [i.set_linewidth(linewidth) for i in ax.spines.values()]
    plt.savefig(os.path.join(directory, 'dispAutoCorreation-allT.png'), dpi=600)
    plt.savefig(os.path.join(directory, 'dispAutoCorreation-allT.eps'), dpi=600)


def plot_dt_displacements(profiles, dts=None, min_dt=None, max_dt=None, spread=None, directory=''):
    """
    Plots cumulative distribution functions of displacements for given durations.
    In the ethos of Banigan 2015, PLoS Computational Biology.
    However, we use CDF rather than probability density functions.
    :param profiles:
    :param dts: if used, a list of keys (in seconds) to plot.
    :param min/max_dt: max and min filters to graph
    :param spread: integer, number of dts to graph, selected as an even spread of those available.
    :param directory:
    :return:
    """
    def acquire_format_data(deltaT_disp, dts):
        plt.clf()
        keys = sorted(deltaT_disp.keys())

        if dts is None:
            dt = keys[0]
            dts = [2 * dt, 4 * dt, 7 * dt, 12 * dt, 40 * dt]

        if dts:
            keys = [k for k in keys if k in dts]  # Filter, use only those keys selected by user.
        else:
            # Currently unreachable code. Looks messy for datasets with many time steps.
            # Filter time intervals to graph.
            if min_dt:
                keys = [k for k in keys if k >= min_dt]
            if max_dt:
                keys = [k for k in keys if k <= max_dt]
            if spread and spread < len(keys):
                # function will will select m items evenly spread over range of zero to v
                f = lambda m, v: [i * v // m + v // (2 * m) for i in range(m)]
                inds = f(spread, len(keys))  # which indices to graph
                # deep copy because we're manipulating positions.
                keys = [keys[i] for i in inds]

        distros = list()  # CDF plots require list of lists, one for each distribution to be plotted
        names = list()
        for key in keys:
            distros.append(deltaT_disp[key])
            names.append('{:.2f}'.format(key / 60.))  # include conversion to minutes
        return distros, names

    distros, names = acquire_format_data(profile.Profile.collate_deltaT_displacements(profiles), dts=dts)
    if len(distros) > 0:
        x_symlog_thresh = 0.5  # Avoid plotting zero on log graph
        plt.clf()
        plot_utils.plotCDFs(distros, names, xlog=True, xmin=x_symlog_thresh, xlabel='Displacement (um)',
                            filename=directory + '/dt_displacements')

        # Plot displacements in component axes. These can have both positive and negative data, so symlog scaling used.
        distros_x, names = acquire_format_data(profile.Profile.collate_deltaT_displacements_X(profiles), dts=dts)
        x_min = 1e0
        x_symlog_thresh = max(np.min(np.abs(np.concatenate(distros_x))), x_min)  # avoid plotting zero on log graph
        plt.clf()
        plot_utils.plotCDFs(distros_x, names, xlog='symlog', x_symlog_thresh=x_symlog_thresh, xlabel='Displacement (um)',
                            filename=directory + '/dt_displacements_X')

        distros_y, names = acquire_format_data(profile.Profile.collate_deltaT_displacements_Y(profiles), dts=dts)
        x_symlog_thresh = max(np.min(np.abs(np.concatenate(distros_y))), x_min)  # avoid plotting zero on log graph
        plt.clf()
        plot_utils.plotCDFs(distros_y, names, xlog='symlog', x_symlog_thresh=x_symlog_thresh, xlabel='Displacement (um)',
                            filename=directory + '/dt_displacements_Y')

        distros_z, names = acquire_format_data(profile.Profile.collate_deltaT_displacements_Z(profiles), dts=dts)
        x_symlog_thresh = max(np.min(np.abs(np.concatenate(distros_z))), x_min)  # avoid plotting zero on log graph
        plt.clf()
        plot_utils.plotCDFs(distros_z, names, xlog='symlog', x_symlog_thresh=x_symlog_thresh, xlabel='Displacement (um)',
                            filename=directory + '/dt_displacements_Z')
    else:
        print("WARNING! Unable to draw dt_displacement graphs. Likely, the dts selected for plotting did not correspond"
              "with the time step of the input data set.")

def plot_representative_tracks(Ps,
                               directory=None, graph_path=None,  # One of these must be selected.
                               n=40, xlim=None, ylim=None, sampling="displacement",
                               tick_interval=None, same_xy_lim=False, same_xy_aspect=True,
                               mark_origin=False, centre_on_origin=True, legend=False, fontsize=18.0, seed=0):
    """
    plots a representative subset of tracks from the given profiles. Track relative displacements are plotted,
    rather than absolute position (hence they are zeroed for t=0).
    Tracks are selected to cover the full range of displacements with equal interval.
    :param directory: filesystem location to save graphs to. Graphs are named automaticallyself.
    :param graph_path: instead of directory, specify the name of the graph. Leave off ".png" or ".eps", it is added automatically.
    :param n: number of tracks to draw
    :param xlim: tuple or list of two items, [min, max], the ranges to show for the xaxis
    :param ylim: tuple or list of two items, [min, max], the ranges to show for the yaxis
    :param sampling: "displacement" yields even sampling of track displacements in dataset. Anything else yields random
            sampling
    :param tick_interval: double, the interval on x and y axes. If set to None, will auto-determine.
    """
    import copy
    linewidth = 2.0    
    all_ts = []
    for prof in Ps:   # Collapse all tracks from all profiles into one list from which representatives are drawn.
        all_ts.extend(prof.tracks)

    random.seed(seed)
    if sampling == "displacement":  # Sets up sampling of tracks to uniformly sample range of displacements.
        all_ts.sort(key=lambda tra: tra.displacement)
    else:
        if len(all_ts) > n:  # Only shuffle if not drawing all the tracks. 
            random.shuffle(all_ts)  # Random sampling of tracks.

    if len(all_ts) > n:
        # Function will will select m items evenly spread over range of zero to v
        f = lambda m, v: [i * v // m + v // (2 * m) for i in range(m)]
        inds = f(n, len(all_ts))
        # Deep copy because we're manipulating positions.
        samples = [copy.deepcopy(all_ts[i]) for i in inds]
    else:
        samples = [copy.deepcopy(all_ts[i]) for i in range(len(all_ts))]
    # Optionally plot all tracks as starting from the origin. 
    if centre_on_origin:
        for t in samples:   # Calculate the relative position from time zero.
            start_x = t.positions[0].x
            start_y = t.positions[0].y
            start_z = t.positions[0].z
            for pos in t.positions:
                pos.x -= start_x
                pos.y -= start_y
                pos.z -= start_z
    plt.clf()
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    # Give every line a unique color
    colors = unique_colour_per_series(len(samples))    
    if sampling == 'displacement':
        # Ensure no correlation between track displacement and color.
        random.shuffle(colors)
    for i, t in enumerate(samples):
        plt.plot([pos.x for pos in t.positions], [pos.y for pos in t.positions], c=colors[i], linewidth=linewidth,
                 label=t.ID)
    if mark_origin:
        plt.plot([0], [0], marker="8", markersize=6, color="black")
    plt.xlabel('X displacement (um)')
    plt.xticks(rotation=45)
    plt.ylabel('Y displacement (um)')
    if tick_interval:
        # Ensure that all plots have the same major tick. This eases comparison between graphs from different datasets.
        start, end = plt.gca().get_xlim()
        plt.gca().xaxis.set_ticks(np.arange(start, end, tick_interval))
        start, end = plt.gca().get_ylim()
        plt.gca().yaxis.set_ticks(np.arange(start, end, tick_interval))
    elif same_xy_lim:
        start_x, end_x = plt.gca().get_xlim()
        start_y, end_y = plt.gca().get_ylim()
        plt.gca().set_xlim([min(start_x, start_y), max(end_x, end_y)])
        plt.gca().set_ylim([min(start_x, start_y), max(end_x, end_y)])
    plt.gca().grid(True, linewidth=linewidth)   # Turn on grid lines.
    # Ensure the same scale on both axes, so aren't compressed.
    if same_xy_aspect:
        plt.gca().set_aspect('equal', adjustable='box')
    if fontsize:
        font = {'size': fontsize}
        plt.rc('font', **font)
    # Change the width of the plot boundary lines
    ax = plt.gca()
    [i.set_linewidth(linewidth) for i in ax.spines.values()]
    if legend:
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, fontsize='xx-small',
                   ncol=5, mode="expand", borderaxespad=0.) 
    if not graph_path:
        graph_path = os.path.join(directory, 'representative_tracks')
    plt.savefig(graph_path + '.png', dpi=300)
    plt.savefig(graph_path + '.eps', dpi=300)


def plot_mi_time(Ps, directory=None):
    """ 
    Scatter plot of meandering index vs duration of track. There is frequently a bias; short tracks have higher MIs
    """
    linewidth = 2.0
    fontsize = 18
    data = []   # list of tuples, {time, mi}
    for prof in Ps:
        for tra in prof.tracks:
            dur = tra.positions[-1].time_m - tra.positions[0].time_m  # calculate duration of track from abs time
            data.append((dur, tra.meander))
    plt.clf()
    times = [d[0] for d in data]
    mis = [d[1] for d in data]
    plt.scatter(times, mis)
    rho, p = scipy.stats.spearmanr(times, mis)
    plt.title('rho = {rho:.2f}'.format(rho=rho))
    plt.gca().set_ylim([0.0, 1.0])
    xmin, xmax = plt.gca().get_xlim()
    plt.gca().set_xlim([0.0, xmax])
    plt.xlabel('Track duration (min)')
    plt.ylabel('Meandering index')
    font = {'size': fontsize}
    plt.rc('font', **font)
    plt.gca().grid(True, linewidth=linewidth)   # turn on grid lines.
    # change the width of the plot boundary lines
    ax = plt.gca()
    [i.set_linewidth(linewidth) for i in ax.spines.values()]
    plt.savefig(os.path.join(directory, 'mi_over_time.png'), dpi=300)
    plt.savefig(os.path.join(directory, 'mi_over_time.eps'), dpi=300)


def plot_displacement_time(Ps, directory=None):
    linewidth = 2.0
    fontsize = 18
    data = []   # list of tuples, {time, mi}
    for prof in Ps:
        for tra in prof.tracks:
            time = tra.positions[-1].time_s - tra.positions[0].time_s  # calculate duration of track from abs time
            time /= 60.0  # convert seconds into minutes
            data.append((time, tra.displacement))
    plt.clf()
    times = [d[0] for d in data]
    displacements = [d[1] for d in data]
    plt.scatter(times, displacements)
    rho, p = scipy.stats.spearmanr(times, displacements)
    plt.title('rho = {rho:.2f}'.format(rho=rho))
    xmin, xmax = plt.gca().get_xlim()
    plt.gca().set_xlim([0.0, xmax])
    plt.xlabel('Track duration (min)')
    plt.ylabel('Displacement (um)')
    font = {'size'   : fontsize}
    plt.rc('font', **font)
    plt.gca().grid(True, linewidth=linewidth)   # turn on grid lines.
    # change the width of the plot boundary lines
    ax = plt.gca()
    [i.set_linewidth(linewidth) for i in ax.spines.values()]
    plt.savefig(os.path.join(directory, 'displacement_over_time.png'), dpi=300)
    plt.savefig(os.path.join(directory, 'displacement_over_time.eps'), dpi=300)


def plot_track_fmi(Ps, directory=None):
    data = []
    for prof in Ps:
        for tra in prof.tracks:
            if tra.fmi is not None:
                data.append(tra.fmi)
    if len(data) > 0:  # do this only if there is something to plot
        plot_utils.plotCDFs([data], xlabel='Forward migration index', title=None,
                            filename=os.path.join(directory, 'track_fmi'), xmax=None, plot_uniform=True)


def plot_instantaneous_fmi(Ps, directory=None):
    """ Plot instantaneous FMI, pooled from all tracks and across all profiles. """
    # data =[]
    data = [i for p in Ps for i in p.collate_instantaneous_fmi()]
    if len(data) > 0:  # do this only if there is something to plot
        plot_utils.plotCDFs([data], xlabel='Instantaneous forward migration index', title=None,
                            filename=os.path.join(directory, 'instant_fmi'), xmax=None, plot_uniform=True)


def plot_track_fmi_biases(Ps, directory=None):
    fmis = []
    lengths = []
    durations = []
    mean_speeds = []
    mean_turns = []
    mi = []  # meandering index
    for prof in Ps:
        for tra in prof.tracks:
            if tra.fmi is not None:
                fmis.append(tra.fmi)
                lengths.append(tra.length)
                durations.append(tra.duration_min)
                # compile mean speeds from positions
                mean_speeds.append(np.mean([pos.speed for pos in tra.positions if pos.speed is not None]))
                mean_turns.append(np.mean([pos.turn for pos in tra.positions if pos.turn is not None]))
                mi.append(tra.meander)
    if len(fmis) > 0:
        plt.clf()
        df = pandas.DataFrame.from_dict({'fmis': fmis, 'lengths': lengths})
        r, p = scipy.stats.pearsonr(fmis, lengths)
        sns.regplot(x='fmis', y='lengths', data=df)
        plt.ylabel('Track length (um)')
        plt.xlabel('Track forward migration index')
        plt.title('Pearsons R = {:.2f}, p = {:.4f}'.format(r, p))
        plt.savefig(directory+'/track_fmi_vs_length.png', dpi=300, bbox_inches='tight')

        plt.clf()
        df = pandas.DataFrame.from_dict({'fmis': fmis, 'durations': durations})
        r, p = scipy.stats.pearsonr(fmis, durations)
        sns.regplot(x='fmis', y='durations', data=df)
        plt.xlabel('Track forward migration index')
        plt.ylabel('Track duration (min)')
        plt.title('Pearsons R = {:.2f}, p = {:.4f}'.format(r, p))
        plt.savefig(directory+'/track_fmi_vs_durations.png', dpi=300, bbox_inches='tight')

        plt.clf()
        df = pandas.DataFrame.from_dict({'fmis': fmis, 'mean_speeds': mean_speeds})
        r, p = scipy.stats.pearsonr(fmis, mean_speeds)
        sns.regplot(x='fmis', y='mean_speeds', data=df)
        plt.xlabel('Track forward migration index')
        plt.ylabel('Mean track speed (um/min)')
        plt.title('Pearsons R = {:.2f}, p = {:.4f}'.format(r, p))
        plt.savefig(directory+'/track_fmi_vs_mean_speeds.png', dpi=300, bbox_inches='tight')

        plt.clf()
        df = pandas.DataFrame.from_dict({'fmis': fmis, 'mean_turns': mean_turns})
        r, p = scipy.stats.pearsonr(fmis, mean_turns)
        sns.regplot(x='fmis', y='mean_turns', data=df)
        plt.xlabel('Track forward migration index')
        plt.ylabel('Mean track turn speed (degrees/min)')
        plt.title('Pearsons R = {:.2f}, p = {:.4f}'.format(r, p))
        plt.savefig(directory+'/track_fmi_vs_mean_turns.png', dpi=300, bbox_inches='tight')

        plt.clf()
        df = pandas.DataFrame.from_dict({'fmis': fmis, 'mi': mi})
        r, p = scipy.stats.pearsonr(fmis, mi)
        sns.regplot(x='fmis', y='mi', data=df)
        plt.xlabel('Track forward migration index')
        plt.ylabel('Meandering index')
        plt.title('Pearsons R = {:.2f}, p = {:.4f}'.format(r, p))
        plt.savefig(directory+'/track_fmi_vs_meander.png', dpi=300, bbox_inches='tight')


def count_tracks_spots(profiles):
    n_tracks = sum([len(p.tracks) for p in profiles])
    n_spots = sum([len(t.positions) for p in profiles for t in p.tracks])
    av = float(n_spots) / float(n_tracks)
    print('{t} tracks in total; {n} spots in total. This gives {av:.2f} spots per track (mean).'
          .format(t=n_tracks, n=n_spots, av=av))


def analyse_broken_path_biases(profiles, directory):
    broken = []
    unbroken = []
    for p in profiles:
        for t in p.tracks:
            if t.broken:
                broken.append(t)
            else:
                unbroken.append(t)
    if not broken:
        print('No broken tracks, cannot draw broken vs unbroken CDF')
        return
    b_m_speed = [t.median_speed() for t in broken]
    u_m_speed = [t.median_speed() for t in unbroken]
    length_ks, _ = scipy.stats.ks_2samp(b_m_speed, u_m_speed)
    plot_utils.plotCDFs([b_m_speed, u_m_speed], names=['broken', 'unbroken'], xlabel='speed (um/min)',
                        filename=os.path.join(directory, 'broken_medSpd'),
                        title='Broken vs unbroken median track speeds')
    b_displacement = [t.displacement for t in broken]
    u_displacement = [t.displacement for t in unbroken]
    displacement_ks, _ = scipy.stats.ks_2samp(b_displacement, u_displacement)
    output = ('KS-differences in lengths = ' + str(length_ks) + '\n' +
              'KS-differences in displacements = ' + str(displacement_ks))
    print(output)
    with open(os.path.join(directory, 'broken_track_stats'), 'w') as f:
        f.write(output)


def unique_colour_per_series(num_series):
    """
    Used to assign unique colours to e.g. tracks in a plot. 
    Defined once here so that it need be changed only here to alter all tracks in all plots simultaneously.
    """
    colours = plt.get_cmap('gist_rainbow')(np.linspace(0, 1.0, num_series))
    return colours


def nicenumber(x, round):
    """ 
    Found online at https://stackoverflow.com/questions/43075617/python-function-equivalent-to-rs-pretty
    2019/09/14.
    Written by stackoverflow user "Bart"
    """
    exp = np.floor(np.log10(x))
    f   = x / 10**exp

    if round:
        if f < 1.5:
            nf = 1.
        elif f < 3.:
            nf = 2.
        elif f < 7.:
            nf = 5.
        else:
            nf = 10.
    else:
        if f <= 1.:
            nf = 1.
        elif f <= 2.:
            nf = 2.
        elif f <= 5.:
            nf = 5.
        else:
            nf = 10.

    return nf * 10.**exp


def pretty(low, high, n):
    """ 
    Found online at https://stackoverflow.com/questions/43075617/python-function-equivalent-to-rs-pretty
    2019/09/14.
    Written by stackoverflow user "Bart"
    """
    range = nicenumber(high - low, False)
    d     = nicenumber(range / (n-1), True)
    miny  = np.floor(low  / d) * d
    maxy  = np.ceil (high / d) * d
    return np.arange(miny, maxy+0.5*d, d)


def pretty_series(series, n):
    return pretty(np.min(series), np.max(series), n)