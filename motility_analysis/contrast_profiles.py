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

    Mark N. Read, 2016

Contrasts two data sets containing motility data. The results are written to an XML file, and to CSV files. Statistical
comparisons are made using the Kolgomorov-Smirnov (ks) statistic, and the Vargha-Delaney "A" statistic. Both are
non-parametric. The A statistic can be interpreted as "the probability that a randomly drawn sample from distribution A
is larger than a randomly drawn sample from distribution B." Another way interpretation is the proportion of
observations in distribution A that are larger than those of distribution B.
"""
import matplotlib
# matplotlib.use('Agg')   # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
# plt.ioff()  # turn off interactive mode
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import os
import math
import numpy as np
# used in writing data as XML files
from xml.etree.ElementTree import Element, SubElement
from xml.etree import ElementTree
from xml.dom import minidom
from scipy.stats import ks_2samp
from . import graphing
from . import profile
from . import plot_utils


def plotMSDs(summaries, filename, settings, plotRaw=True, plotFit=False):
    """
    Plot two mean square displacements (over time) on the same graph.
    plotRaw and plotFit will draw the raw data (only makes sense if the time steps for all tracks in the raw data are
    the same) and the fitted straight line respectively.
    """
    plt.clf()
    if settings.msd_xlim:
        plt.xlim(settings.msd_xlim)
    if settings.msd_ylim:
        plt.ylim(settings.msd_ylim)

    times = []
    msds = []
    for s in summaries:        
        t = [r[0] for r in s.msdResults['msd']]
        times.append(t)
        m = [r[1] for r in s.msdResults['msd']]
        msds.append(m)

    plot_lines = []  # for generating legend
    plot_labels = []

    if plotRaw:
        for i, s in enumerate(summaries):
            p, = plt.plot(times[i], msds[i], color=s.colour, linewidth=settings.plot_linewidth)
            plot_lines.append(p)
            plot_labels.append('{n:s} {sl:.2f}'.format(n=s.label, sl=s.msdResults['slope']))

    if plotFit:
        for i, s in enumerate(summaries):
            p, = plt.plot(s.msdResults['linearTimes'], s.msdResults['linearPlot'], s.colour,
                          linewidth=settings.plot_linewidth)
            plot_lines.append(p)
            plot_labels.append('{n:s} {sl:.2f}'.format(n=s.label, sl=s.msdResults['slope']))

    plt.legend(plot_lines, plot_labels, loc='upper left')
    plt.xlabel('Time (s)')
    plt.ylabel('MSD')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Mean Squared Displacements')
    plt.gca().grid(True, linewidth=settings.plot_linewidth)   # turn on grid lines.
    font = {'size': settings.plot_textsize}
    plt.rc('font', **font)
    # change the width of the plot boundary lines
    ax = plt.gca()
    [i.set_linewidth(settings.plot_linewidth) for i in ax.spines.values()]
    plt.savefig(filename + '.png', dpi=600)
    plt.savefig(filename + '.svg', dpi=600)


def plot_dt_distplacements(summaries, filename, settings):
    """
    see Banigan 15, PLOS Computational Biology, for details of this metric.
    Note that we use cumulative distribution plots, rather than probability density plots.

    There's a lot of information on this plot.
    A distribution is drawn for each of the given summaries, at each of the specified (settings, if set, available
    otherwise) time intervals.
    Summaries are indicated with line styles, time intervals with colours.
    There are only (currently) 4 line styles specified, and they will repeat if more summaries are supplied.
    """
    plt.clf()
    dt_disps = [s.collate_dt_displacements() for s in summaries]  # list( dict[dt] = list( displacements ) )

    # get keys (dts) common to all supplied profiles
    dts = set(dt_disps[0].keys())
    for dt_disp in dt_disps[1:]:
        dts = dts.intersection(dt_disp.keys())
    dts = sorted(list(dts))
    if settings.dt_disp_dts:  # Retain those selected by the user, if selected. 
        dts = [dt for dt in dts if dt in settings.dt_disp_dts]

    disps = []  # list(list(list(displacements)))  # first (outer) list is summary. Second is dt.
    for i, dt_disp in enumerate(dt_disps):
        disps.append(list())
        for dt in dts:
            if dt in dt_disp.keys():
                disps[i].append(dt_disp[dt])
            else:
                disps[i].append(list())  # empty list if this time interval not present

    line_styles = ['-', ':', '--', '-.']
    # select evenly spaced colours from the given spectrum
    cm = pylab.get_cmap('gist_rainbow')
    colours = list(cm(1.*i/len(dts)) for i in range(len(dts)))
    plot_lines = []
    plot_labels = []
    plt.clf()
    for i, s in enumerate(summaries):
        ls = line_styles[i % len(line_styles)]
        for j, dt in enumerate(dts):
            cdf = plot_utils.ecdf(disps[i][j])
            # insert one record at start of these lists to ensure curves extend to y=0.
            xs = [cdf[0][0]]
            xs.extend([record[0] for record in cdf])
            ys = [0]
            ys.extend([record[1] for record in cdf])
            p, = plt.plot(xs, ys, color=colours[j], linestyle=ls, linewidth=settings.plot_linewidth)
            # Only for the first dt, so label will show line style against label.
            # Too many combinations to show colours for dts also. These increase, so can be ascertained post hoc
            if j == 0:
                plot_lines.append(p)
                plot_labels.append(s.label)
    plt.legend(plot_lines, plot_labels, loc='lower right')
    plt.xscale('log')
    x_lim = (settings.dt_disp_min_displacement, plt.xlim()[1])
    plt.xlim(x_lim)
    plt.ylabel('Cumulative distribution')
    plt.xlabel('Log displacement')
    plt.gca().grid(True, linewidth=settings.plot_linewidth)  # turn on grid lines.
    font = {'size': settings.plot_textsize}
    plt.rc('font', **font)
    ax = plt.gca()
    # change the width of the plot boundary lines
    [i.set_linewidth(settings.plot_linewidth) for i in ax.spines.values()]

    plt.savefig(filename + '.png', dpi=600)
    plt.savefig(filename + '.svg', dpi=600)


def plot_disp_autocorr_diffs(summaries, filename, settings, 
                             dataseries=lambda s: s.collate_displacement_autocorrelation()[0]):
    """ see Banigan 15, PLOS Computational Biology, for example of how this is graphed. """

    def calc_errs(data):
        """ these calculations need to be performed several times. Formatting data for plotting. """
        upper_err = []  # length of upper error bar from median (not the absolute value itself)
        lower_err = []  # length of lower error bar from median (not absolute value itself)
        median_data = []  # the median data to be plotted.
        for d in data:
            median = np.median(d)
            err = np.percentile(a=d, q=settings.error_bar_quartiles)  # range of distribution
            lower_err.append(median - err[0])  # find distance of these points from median
            upper_err.append(err[1] - median)
            median_data.append(median)
        lower_err = np.asarray(lower_err)  # life is easier working with numpy arrays
        upper_err = np.asarray(upper_err)
        # can't plot negative values on log plot. This cuts off error values at a given threshold.
        ylower = np.maximum(settings.dac_mincorr, median_data - lower_err)
        lower_err = median_data - ylower
        return median_data, lower_err, upper_err

    dacs = [dataseries(s) for s in summaries]

    # get keys (dts) common to all supplied profiles
    dts = set(dacs[0].keys())
    for dac in dacs[1:]:
        dts = dts.intersection(dac.keys())
    dts = sorted(list(dts))
    dts = [dt for dt in dts if settings.dac_min_dt <= dt <= settings.dac_max_dt]
    datas = []
    for dac in dacs:
        data = [dac[dt] for dt in dts]
        datas.append(data)

    median = list()  # needed for box plots of distributions
    lower_err = list()
    upper_err = list()
    for data in datas:
        m, l, u = calc_errs(data)
        median.append(m)
        lower_err.append(l)
        upper_err.append(u)

    dts_min = [dt / 60. for dt in dts]  # convert from seconds to minutes
    # Shift each box plot slightly for clarify, prevent them being drawn atop one another
    # Spread defines the range over which boxes are shifted around their actual x location
    # This is defined as 15% of the range between each time interavl, and is then spread evenly between each box
    interval = dts_min[1] - dts_min[0] if len(dts_min) > 1 else 1.0
    total_spread = 0.25 * interval
    ind_spread = total_spread / (len(summaries) - 1)  # 1 gap for 2 items
    x_values = []  # list of lists. The x-values at which to plot boxes, for each data set
    for i in range(len(summaries)):
        xs = [dt - (total_spread / 2.) + (i * ind_spread) for dt in dts_min]
        x_values.append(xs)

    plt.clf()
    for i, s in enumerate(summaries):
        plt.errorbar(x_values[i], median[i], yerr=[lower_err[i], upper_err[i]], color=s.colour,
                     fmt='o', label=s.label, elinewidth=settings.plot_linewidth)

    plt.legend()
    min_y = max(plt.ylim()[0], settings.dac_mincorr)
    plt.ylim((min_y, 1.0))
    plt.yscale('log')
    plt.xlabel('Time interval (min)')
    plt.ylabel('Displacement autocorrelation')
    plt.gca().grid(True, linewidth=settings.plot_linewidth)   # turn on grid lines.
    font = {'size': settings.plot_textsize}
    plt.rc('font', **font)
    ax = plt.gca()
    # change the width of the plot boundary lines
    [i.set_linewidth(settings.plot_linewidth) for i in ax.spines.values()]

    plt.savefig(filename + '.png', dpi=600)
    plt.savefig(filename + '.svg', dpi=600)


def plot_msd_vs_max_dt_diffs(datasets, filename, names=None, colours=None):
    """

    :param datasets: list of tuples. Each tuple is ([max dts], [slopes])
    :param names: list of strings, one per item in datasets
    :param colours: list of colours, one per item in datasets
    :param filename:
    """
    plt.clf()
    plot_lines = []
    for ds in datasets:
        max_dts = ds[0]
        slopes = ds[1]
        p_l, = plt.plot(max_dts, slopes)
        plot_lines.append(p_l)
    if names is not None:
        plt.legend(plot_lines, names, loc='lower left')
    plt.xlabel('Maximum time interval (s)')
    plt.ylabel('MSD slope')
    plt.ylim([1.0, 1.5])
    plt.gca().grid(True)  # turn on grid lines.
    plt.savefig(filename, bbox_inches='tight', transparent=True, dpi=600)


def plot_track_displacement_angle_xy_diffs(summaries, filename, colours, bins=36):
    """ Rose plot showing the orientation of track displacement in the xy plane as angles. """    
    # re: endpoint, need last bin location for np.histogram
    # bins + 1 here because this specifies boundary locations, of which there are one more than the number of bins 
    # (last bin has start AND end location)
    bin_locations = np.linspace(0.0, 2 * np.pi, bins + 1, endpoint=True)  
    bin_width = (2 * np.pi) / bins
    
    angles = summaries[0].collate_displacement_angles_xy()
    bin_hist, _ = np.histogram(angles, bin_locations)
    bin_hist = bin_hist / sum(bin_hist)  # Turn into a probability density
    
    bottom = 0  # Bin does not need to start from zero. Setting this to non-zero creates a "donut"

    ax = plt.subplot(111, polar=True)
    bin_start_theta = bin_locations[:-1]  # Exclude the last location. 
    ax.bar(bin_start_theta, bin_hist, width=bin_width, bottom=bottom, edgecolor='black', facecolor=colours[0], 
           label=summaries[0].label)

    for summary, col in zip(summaries[1:], colours[1:]):
        angles = summary.collate_displacement_angles_xy()
        bin_hist, _ = np.histogram(angles, bin_locations)
        bin_hist = bin_hist / sum(bin_hist)  # Turn into a probability density. 
        bin_hist = list(bin_hist)
        # Step function needs first item appended onto end otherwise the radial plot does not go completely round. 
        bin_hist.append(bin_hist[0])

        ax.step(bin_locations, bin_hist, where='mid', color=col, linewidth=2, label=summary.label)

    ax.grid(True)
    ax.set_yticklabels([])
    ax.set_axisbelow(True)  # Draw grid behind bars. 
    ax.legend()

    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.savefig(filename, bbox_inches='tight')


class ProfileSummary:
    """
    Collates statistical profiles from one or more profiles.
    """
    def __init__(self, profiles, settings, label=None, colour=None):
        """
        Collect together stats for a given set of profiles.
        !!!NOTE that MSD calculations are performed twice. First, here, to extract MSD-dt values for all available
        dt values (no cutoff). When this ProfilesSummary is contrasted with another, a common cutoff value is chosen
        to suit both, and MSD slopes are calculated again (based on filtered data).
        """
        self.label = label
        self.colour = colour  # used for plotting
        self.profiles = profiles
        self.settings = settings
        # Deferred calculations on request
        self._dur = None  
        self._len = None  
        self._dis = None  
        self._mea = None  
        self._spd = None  
        self._mspd = None  
        self._iqrSpd = None
        self._turn = None
        self._mTurn = None
        self._iqrTurn = None
        self._roll = None
        self._mRoll = None
        self._iqrRoll  = None
        self._fmi = None
        self._instant_fmi = None
        self._displacement_angles_xy = None

        # MSD is calculated twice. First over all possible dt values (no cutoff). An appropriate cutoff is then
        # selected from two (or more?) profiles together, and each profile's MSD is then re-calculated using this
        # cutoff.
        self._msdData = None
        self._lastTimeSample = None
        self.msdResults = None  # This is set after comparison with other profile.

        # Displacement autocorrelation
        self._disp_autocorr = None        
        self._normalised_dac = None

        self._dt_displacements = None

    def collate_durations(self):
        if self._dur is None:
            self._dur = [i for p in self.profiles for i in p.duration]
        return self._dur
    
    def collate_lengths(self):
        if self._len is None:
            self._len = [i for p in self.profiles for i in p.collate_lengths()]
        return self._len
    
    def collate_displacements(self):
        if self._dis is None:
            self._dis = [i for p in self.profiles for i in p.displacement]
        return self._dis
    
    def collate_meandering_indices(self):
        if self._mea is None:
            self._mea = [i for p in self.profiles for i in p.collate_meanders()]
        return self._mea

    def collate_speeds(self):
        if self._spd is None:
            self._spd = [i for p in self.profiles for i in p.collate_speeds()]
        return self._spd

    def collate_median_speeds(self):
        if self._mspd is None:
            self._mspd = [i for p in self.profiles for i in p.medSpd]
        return self._mspd
    
    def collate_iqr_speeds(self):
        if self._iqrSpd is None:
            self._iqrSpd = [i for p in self.profiles for i in p.irqSpd]
        return self._iqrSpd

    def collate_turns(self):
        if self._turn is None:
            self._turn = [i for p in self.profiles for i in p.collate_turns()]
        return self._turn

    def collate_median_turns(self):
        if self._mTurn is None:
            self._mTurn = [i for p in self.profiles for i in p.medTurn]
        return self._mTurn

    def collate_iqr_turns(self):
        if self._iqrTurn is None:
            self._iqrTurn = [i for p in self.profiles for i in p.irqTurn]
        return self._iqrTurn

    def collate_rolls(self):
        if self._roll is None:
            self._roll = [i for p in self.profiles for i in p.collate_rolls()]
        return self._roll

    def collate_median_rolls(self):
        if self._mRoll is None:
            self._mRoll = [i for p in self.profiles for i in p.medRoll]
        return self._mRoll

    def collate_iqr_rolls(self):
        if self._iqrRoll is None:
            self._iqrRoll = [i for p in self.profiles for i in p.irqRoll]
        return self._iqrRoll

    def collate_fmis(self):
        if self._fmi is None:
            self._fmi = [t.fmi for p in self.profiles for t in p.tracks if t.fmi is not None]
        return self._fmi

    def collate_instant_fmis(self):
        if self._instant_fmi is None:
            self._instant_fmi = [i for p in self.profiles for i in p.collate_instantaneous_fmi()]
        return self._instant_fmi

    def collate_displacement_angles_xy(self):
        if self._displacement_angles_xy is None:
            self._displacement_angles_xy = [t.get_arctan_xy() for p in self.profiles for t in p.tracks]
        return self._displacement_angles_xy

    def collate_msd(self):
        if self._msdData is None or self._lastTimeSample is None:
            msd_results = profile.Profile.calculate_msd(self.profiles, max_dt=float('inf'), method='allT')
            self._msdData = msd_results['msd']  # Contains tuple (<time><msd value>)
            self._lastTimeSample = self._msdData[-1][0]
        return self._msdData, self._lastTimeSample

    def collate_displacement_autocorrelation(self):
        if self._disp_autocorr is None:
            self._disp_autocorr = profile.Profile.collate_displacement_autocorrelation(self.profiles)
            # Normalised displacement autocorrelation (DAC).
            # To mitigate effect of "shudder" (random noise around positional data, owing to either instrument or changes
            # in cell shape), we adjust all DAC data such that the smallest interval (dt)'s DAC value is centred on 1.
            # The minimum dt specified in the settings should be large enough that the longer term DAC exponential decay
            # (seen as linear trend on log(DAC) graphs) has commenced.
            self._normalised_dac = dict()  # same format as disp_autcorr: dict[dt] = list(dac values)
            dts = sorted(self._disp_autocorr.keys())
            # Filter values to use
            dts = [dt for dt in dts if self.settings.dac_min_dt <= dt <= self.settings.dac_max_dt]  
            # Add quantity to each DAC value such that the median of the smallest interval (dt) = 1.
            smallest_dt = np.min(dts)
            min_dt_median_dac = np.median(self._disp_autocorr[smallest_dt])
            adjustment = 1. - min_dt_median_dac
            for dt in dts:
                adj = [v + adjustment for v in self._disp_autocorr[dt]]  # normalise to 1
                self._normalised_dac[dt] = adj
        return self._disp_autocorr, self._normalised_dac

    def collate_dt_displacements(self):
        if self._dt_displacements is None:
            self._dt_displacements = profile.Profile.collate_deltaT_displacements(self.profiles)
        return self._dt_displacements


class ProfilesContrastStats:
    def __init__(self, ps1, ps2, settings):
        """
        msd_time_cutoff value of None means unspecified, and will be calculated to be 25% of longest track found in
        either profile. Another value can be used instead.
        """
        # Kolmogorov-Smirnov test statistics.
        if settings.duration:
            self.dur_KSd, self.dur_KSp = ks_2samp(ps1.collate_durations(), ps2.collate_durations())
        else:
            self.dur_KSd, self.dur_KSp = float('nan'), float('nan')

        if settings.length:
            self.len_KSd, self.len_KSp = ks_2samp(ps1.collate_lengths(), ps2.collate_lengths())
        else:
            self.len_KSd, self.len_KSp = float('nan'), float('nan')

        if settings.displacement:
            self.dis_KSd, self.dis_KSp = ks_2samp(ps1.collate_displacements(), ps2.collate_displacements())
        else:
            self.dis_KSd, self.dis_KSp = float('nan'), float('nan')
        
        if settings.meander:
            self.mea_KSd, self.mea_KSp = ks_2samp(ps1.collate_meandering_indices(), ps2.collate_meandering_indices())
        else:
            self.mea_KSd, self.mea_KSp = float('nan'), float('nan')
        
        if settings.speeds:
            self.spd_KSd, self.spd_KSp = ks_2samp(ps1.collate_speeds(), ps2.collate_speeds())
            self.mspd_KSd, self.mspd_KSp = ks_2samp(ps1.collate_median_speeds(), ps2.collate_median_speeds())
            self.irqSpd_KSd, self.irqSpd_KSp = ks_2samp(ps1.collate_iqr_speeds(), ps2.collate_iqr_speeds())
        else:
            self.spd_KSd, self.spd_KSp = float('nan'), float('nan')
            self.mspd_KSd, self.mspd_KSp = float('nan'), float('nan')
            self.irqSpd_KSd, self.irqSpd_KSp = float('nan'), float('nan')
        
        if settings.turns:
            self.turn_KSd, self.turn_KSp = ks_2samp(ps1.collate_turns(), ps2.collate_turns())
            self.mTurn_KSd, self.mTurn_KSp = ks_2samp(ps1.collate_median_turns(), ps2.collate_median_turns())
            self.irqTurn_KSd, self.irqTurn_KSp  = ks_2samp(ps1.collate_iqr_turns(), ps2.collate_iqr_turns())
        else:
            self.turn_KSd, self.turn_KSp = float('nan'), float('nan') 
            self.mTurn_KSd, self.mTurn_KSp = float('nan'), float('nan') 
            self.irqTurn_KSd, self.irqTurn_KSp = float('nan'), float('nan') 
        
        if settings.rolls:
            self.roll_KSd, self.roll_KSp = ks_2samp(ps1.collate_rolls(), ps2.collate_rolls())
            self.mRoll_KSd, self.mRoll_KSp = ks_2samp(ps1.collate_median_rolls(), ps2.collate_median_rolls())
            self.irqRoll_KSd, self.irqRoll_KSp  = ks_2samp(ps1.collate_iqr_rolls(), ps2.collate_iqr_rolls())
        else:
            self.roll_KSd, self.roll_KSp = float('nan'), float('nan')
            self.mRoll_KSd, self.mRoll_KSp = float('nan'), float('nan')
            self.irqRoll_KSd, self.irqRoll_KSp = float('nan'), float('nan')

        self.fmi_KSd, self.fmi_KSp = float('nan'), float('nan')
        self.inst_fmi_KSd, self.inst_fmi_KSp = float('nan'), float('nan')
        if settings.fmi:
            # Both have to be >0 for comparison to make sense
            if np.min([len(ps1.collate_fmis()), len(ps2.collate_fmis())]) > 0: 
                self.fmi_KSd, self.fmi_KSp = ks_2samp(ps1.collate_fmis(), ps2.collate_fmis())
                self.inst_fmi_KSd, self.inst_fmi_KSp = ks_2samp(ps1.collate_instant_fmis(), ps2.collate_instant_fmis())
        
        if settings.displacement_angle_xy:
            self.displacement_angles_xy_KSd, self.displacement_angles_xy_KSp = ks_2samp(ps1.collate_displacement_angles_xy(), 
                                                                                        ps2.collate_displacement_angles_xy())
        else: 
            self.displacement_angles_xy_KSd, self.displacement_angles_xy_KSp = float('nan'), float('nan')
        
        # Other form of statistics
        # Correlation between length and displacement is dependent on the range of lengths input - higher correlation
        # for shorter lengths. Simulation track lengths tend to be longer than vivo, because simulation does not have
        # broken tracks. Hence, calculate correlations only for the range of track lengths generated by simulation.
        # The following finds the minimum length in each of the sets of profiles.
        if settings.length_displacement_corr:
            ps1min = 0.0
            ps1max = 20.0
            ps2min = 0.0
            ps2max = 20.0
            if ps1.collate_lengths():
                ps1min = min(ps1.collate_lengths())
                ps1max = max(ps1.collate_lengths())
            if ps2.collate_lengths():
                ps2min = min(ps2.collate_lengths())
                ps2max = max(ps2.collate_lengths())
            filter_min = max(ps1min, ps2min)
            filter_max = min(ps1max, ps2max)
            self.ps1_l_disp_corr = graphing.plot_length_displacement_scat(ps1.profiles, filter_min=filter_min,
                                                                        filter_max=filter_max, graph=False)
            self.ps2_l_disp_corr = graphing.plot_length_displacement_scat(ps2.profiles, filter_min=filter_min,
                                                                        filter_max=filter_max, graph=False)
            self.length_displacement_corr_diff = abs(self.ps1_l_disp_corr['rho'] - self.ps2_l_disp_corr['rho'])
        else:
            self.ps1_l_disp_corr = float('nan')
            self.ps2_l_disp_corr = float('nan')
            self.length_displacement_corr_diff = float('nan')

        if settings.msd:
            msd_data1, _ = ps1.collate_msd()
            msd_data2, _ = ps2.collate_msd()
            ps1.msdResults = profile.Profile.calculate_msd(msd=msd_data1, max_dt=settings.msd_max_dt)
            ps2.msdResults = profile.Profile.calculate_msd(msd=msd_data2, max_dt=settings.msd_max_dt)
            self.msdDiff = abs(ps1.msdResults['slope'] - ps2.msdResults['slope'])
        else:
            self.msdDiff = float('nan')

        # Displacement autocorrelation stats. For a given profile, each deltaT associates with a distribution of
        # correlation values. We perform a KS test at each, and then take the sum of these to derive a single metric.
        if settings.disp_autocorr:
            self.disp_autocorr_ks_cumulative = 0.
            ps1_dac, ps1_normalised_dac = ps1.collate_displacement_autocorrelation()
            ps2_dac, ps2_normalised_dac = ps2.collate_displacement_autocorrelation()
            dts = [dt for dt in ps1_dac.keys() if dt in ps2_dac.keys()]  # deltaT values common to both profiles
            dts = [dt for dt in dts if settings.dac_min_dt <= dt <= settings.dac_max_dt]  # filter values to use
            for dt in dts:
                acs1 = ps1_dac[dt]
                acs2 = ps2_dac[dt]
                d, p = ks_2samp(acs1, acs2)
                self.disp_autocorr_ks_cumulative = self.disp_autocorr_ks_cumulative + d

            # Contrast normalised DAC values.
            self.normalised_dac_ks_cumulative = 0.
            for dt in dts:
                d, p = ks_2samp(ps1_normalised_dac[dt], ps2_normalised_dac[dt])
                self.normalised_dac_ks_cumulative = self.normalised_dac_ks_cumulative + d
        else:
            self.disp_autocorr_ks_cumulative = float('nan')
            self.normalised_dac_ks_cumulative = float('nan')

        # Displacement distributions for given time intervals.
        # The sum of KS 'd' values contrasting distributions of displacements for given time intervals.
        if settings.dt_disp:
            self.dt_displacement_ks_cumulative = 0.
            # Find dts common to both data sets
            dts = [dt for dt in ps1.collate_dt_displacements().keys() if dt in ps2.collate_dt_displacements().keys()]
            if settings.dt_disp_dts:  # and retain only those selected by user
                dts = [dt for dt in dts if dt in settings.dt_disp_dts]
            for dt in dts:
                d, p = ks_2samp(ps1.collate_dt_displacements()[dt], ps2.collate_dt_displacements()[dt])
                self.dt_displacement_ks_cumulative = self.dt_displacement_ks_cumulative + d
        else:
            self.dt_displacement_ks_cumulative = float('nan')


    def write_to_file(self, file_name, label1, label2):
        f = open(file_name, 'w')
        f.write('# ' + label1 + ' vs ' + label2 + '\n')
        f.write('{:15s}, {:6s}\n'.format('metric', 'Dval'))
        f.write('{:15s}, {:6f}\n'.format('duration', self.dur_KSd))
        f.write('{:15s}, {:6f}\n'.format('length', self.len_KSd))
        f.write('{:15s}, {:6f}\n'.format('displacement', self.dis_KSd))
        f.write('{:15s}, {:6f}\n'.format('meander', self.mea_KSd))
        f.write('{:15s}, {:6f}\n'.format('speed', self.spd_KSd))
        f.write('{:15s}, {:6f}\n'.format('median spd', self.mspd_KSd))
        f.write('{:15s}, {:6f}\n'.format('IRQ spd', self.irqSpd_KSd))
        f.write('{:15s}, {:6f}\n'.format('turns', self.turn_KSd))
        f.write('{:15s}, {:6f}\n'.format('median turns', self.mTurn_KSd))
        f.write('{:15s}, {:6f}\n'.format('irq turns', self.irqTurn_KSd))
        f.write('{:15s}, {:6f}\n'.format('rolls', self.roll_KSd))
        f.write('{:15s}, {:6f}\n'.format('median rolls', self.mRoll_KSd))
        f.write('{:15s}, {:6f}\n'.format('irq rolls', self.irqRoll_KSd))
        f.write('{:15s}, {:6f}\n'.format('track FMI', self.fmi_KSd))
        f.write('{:15s}, {:6f}\n'.format('instant FMI', self.inst_fmi_KSd))
        f.write('{:15s}, {:6f}\n'.format('DAC cum ks', self.disp_autocorr_ks_cumulative))
        f.write('{:15s}, {:6f}\n'.format('Normal DAC cum ks', self.normalised_dac_ks_cumulative))
        f.write('{:15s}, {:6f}\n'.format('dt disp cum ks', self.dt_displacement_ks_cumulative))
        f.write('{:15s}, {:6f}\n'.format('disp angle XY', self.displacement_angles_xy_KSd))


    def write_response_stats(self, file_name):
        """
        writes an XML file to the file system containing the stats contrasting these two profiles.
        """
        def prettify(elem):
            """
            Return a pretty-printed XML string for the Element.

            Taken from: http://pymotw.com/2/xml/etree/ElementTree/create.html on 19/08/2014
            """
            rough_string = ElementTree.tostring(elem, 'utf-8')
            re_parsed = minidom.parseString(rough_string)
            return re_parsed.toprettyxml(indent="  ")

        def setup_response(name, parent, ksd, ksp=None):
            """
            Creates tags representing response data. Includes data for KS tests (p and d) and Vargha-Delaney A test.
            """
            child = SubElement(parent, name)
            child_ks_d = SubElement(child, 'ks-d')
            child_ks_d.text = str(ksd)
            if ksp is not None:
                child_ks_p = SubElement(child, 'ks-p')
                child_ks_p.text = str(ksp)

        responses = Element('responses')
        setup_response(name='duration',     parent=responses, ksd=self.dur_KSd,     ksp=self.dur_KSp)
        setup_response(name='displacement', parent=responses, ksd=self.dis_KSd,     ksp=self.dis_KSp)
        setup_response(name='length',       parent=responses, ksd=self.len_KSd,     ksp=self.len_KSp)
        setup_response(name='speed',        parent=responses, ksd=self.spd_KSd,     ksp=self.spd_KSp)
        setup_response(name='speedMedian',  parent=responses, ksd=self.mspd_KSd,    ksp=self.mspd_KSp)
        setup_response(name='speedIRQ',     parent=responses, ksd=self.irqSpd_KSd,  ksp=self.irqSpd_KSp)
        setup_response(name='meander',      parent=responses, ksd=self.mea_KSd,     ksp=self.mea_KSp)
        setup_response(name='turns',        parent=responses, ksd=self.turn_KSd,    ksp=self.turn_KSp)
        setup_response(name='turnsMedian',  parent=responses, ksd=self.mTurn_KSd,   ksp=self.mTurn_KSp)
        setup_response(name='turnsIRQ',     parent=responses, ksd=self.irqTurn_KSd, ksp=self.irqTurn_KSp)
        setup_response(name='rolls',        parent=responses, ksd=self.roll_KSd,    ksp=self.roll_KSp)
        setup_response(name='rollsMedian',  parent=responses, ksd=self.mRoll_KSd,   ksp=self.mRoll_KSp)
        setup_response(name='rollsIRQ',     parent=responses, ksd=self.irqRoll_KSd, ksp=self.irqRoll_KSp)
        setup_response(name='trackFMI',  parent=responses, ksd=self.fmi_KSd, ksp=self.fmi_KSp)
        setup_response(name='instantFMI', parent=responses, ksd=self.inst_fmi_KSd, ksp=self.inst_fmi_KSp)
        setup_response(name='displacement_angle_xy', parent=responses, ksd=self.displacement_angles_xy_KSd, 
                                                                       ksp=self.displacement_angles_xy_KSp)

        # MSD data.
        msdNode = SubElement(responses, 'msdDiff')
        msdNode.text = str(self.msdDiff)
        # Displacement autocorrelations
        dacNode = SubElement(responses, 'disp_autocorr_cumulative_KS')
        dacNode.text = str(self.disp_autocorr_ks_cumulative)
        # Normalised displacement autocorrelations
        norm_dac_node = SubElement(responses, 'normalised_DAC_cumulative_KS')
        norm_dac_node.text = str(self.normalised_dac_ks_cumulative)
        # dt displacements
        dt_disp_node = SubElement(responses, 'dt_disp_cumulative_KS')
        dt_disp_node.text = str(self.dt_displacement_ks_cumulative)
        #
        child = SubElement(responses, 'length_displacement_corr')
        c1 = SubElement(child, 'ps1')
        c2 = SubElement(child, 'ps2')
        c3 = SubElement(child, 'diff')
        c1.text = str(self.ps1_l_disp_corr)
        c2.text = str(self.ps2_l_disp_corr)
        c3.text = str(self.length_displacement_corr_diff)

        with open(file_name, 'w') as f:
            f.write(prettify(responses))


class ContrastSettings():
    """
    There are many ways profiles can be contrasted.
    Rather than running all analyses, this object stores context/problem-specific options (which should not be
    hard-coded into this module).
    """
    def __init__(self,
                 dt=None,  # Autodetect. 
                 disp_autocorr=True,  # Displacement autocorrelation
                 disp_autocorr_dts=None,  # What intervals to use
                 dac_min_dt=0.,  # Seconds
                 dac_max_dt=300,  # Seconds
                 dac_mincorr=1e-3,
                 # Note, can't perform this analysis (i.e. it will likely crash) if constituent profiles differ in 
                 # their time steps. 
                 dt_disp=True,
                 dt_disp_dts=None,  # List of items, in seconds
                 dt_disp_min_displacement=1.,
                 msd=True,
                 msd_max_dt=float('inf'),
                 msd_xlim=None,  # Tuple of two items, if used
                 msd_ylim=None,
                 speeds=True,
                 turns=True,
                 rolls=True,
                 displacement=True,  # From track start to end. Not recommended, as confounded by variable duration
                 displacement_angle_xy=True,
                 duration=True,
                 meander=True,
                 length=True,
                 length_displacement_corr=False,
                 fmi=True,  # Forward migration index                 
                 error_bar_quartiles=(25, 75),
                 plot_linewidth=2.,
                 plot_textsize=18
                 ):
        self.disp_autocorr = disp_autocorr
        self.disp_autocorr_dts = disp_autocorr_dts
        self.dac_min_dt = dac_min_dt
        self.dac_max_dt = dac_max_dt
        self.dac_mincorr = dac_mincorr
        self.dt_disp = dt_disp        
        if dt_disp_dts:            
            self.dt_disp_dts = dt_disp_dts
        elif dt:
            self.dt_disp_dts = [2 * dt, 4 * dt, 7 * dt, 12 * dt, 40 * dt]
        else:
        	self.dt_disp_dts = dt_disp_dts

        self.dt_disp_min_displacement = dt_disp_min_displacement
        self.msd = msd
        self.msd_xlim = msd_xlim
        self.msd_ylim = msd_ylim
        self.msd_max_dt = msd_max_dt
        self.speeds = speeds
        self.turns = turns
        self.rolls = rolls
        self.displacement = displacement
        self.displacement_angle_xy = displacement_angle_xy
        self.duration = duration
        self.meander = meander
        self.length = length
        self.length_displacement_corr = length_displacement_corr
        self.fmi = fmi
        self.error_bar_quartiles = error_bar_quartiles
        self.plot_linewidth = plot_linewidth
        self.plot_textsize = plot_textsize


class NoTracksException(Exception):
    pass


def _setup_profile(prof, label, colour, settings):
    if isinstance(prof, ProfileSummary):
        return prof

    if isinstance(prof, profile.Profile):
        prof = [prof]

    if [p.tracks for p in prof]:  # ensure there are tracks in the profile
        summary = ProfileSummary(prof, settings=settings, label=label, colour=colour)
    else:
        raise NoTracksException('WARNING! one of the supplied profiles contains no tracks.')
    return summary


def contrast(profile1, profile2, label1, label2, out_dir=None, draw_graphs=False, profile3=None,
             label3=None, p1_colour='b', p2_colour='g', p3_colour='r',
             p1_p2_contrast_name=None,
             p1_p3_contrast_name=None, settings=ContrastSettings()):
    """
    Contrast the two supplied profiles. An optional 3rd distribution can be provided, and it is contrasted to the first.


    :param profile1: first of the profiles to contrast. Can be a list of profiles, or ProfileSummary object
    :param profile2: second profile to contrast. Can be a list of profiles, or ProfileSummary object
    :param label1: name to give first profile for graphing.
    :param label2: name to give second profile for graphing.
    :param out_dir: directory to which output data should be written.
    :param draw_graphs: time consuming to draw graphs, set to True if you wish to have them.
    """
    print('contrasting profiles.')
    if p1_p2_contrast_name is None:
        p1_p2_contrast_name = 'contrast_{:s}_{:s}_statistics.xml'.format(label1, label2)
    if p1_p3_contrast_name is None and profile3 is not None:
        p1_p3_contrast_name = 'contrast_{:s}_{:s}_statistics.xml'.format(label1, label3)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    summary1 = _setup_profile(profile1, label1, colour=p1_colour, settings=settings)
    summary2 = _setup_profile(profile2, label2, colour=p2_colour, settings=settings)
    summaries = [summary1, summary2]
    if profile3 is not None:
        summary3 = _setup_profile(profile3, label3, colour=p3_colour, settings=settings)
        summaries.append(summary3)
    else:
        summary3 = None
        contrast1_3 = None

    # Contrast various summaries
    # Convention to calculate MSD only over short times. Here, 25% of the lowest max dt.
    if settings.msd:
        last_time_samples = []
        for s in summaries:
            _, lts = s.collate_msd()
            last_time_samples.append(lts)
        settings.msd_max_dt = 0.25 * np.min(last_time_samples)
        # settings.msd_max_dt = 0.25 * np.min([s._lastTimeSample for s in summaries])

    contrast1_2 = ProfilesContrastStats(summary1, summary2, settings)
    contrast1_2.write_to_file(file_name=os.path.join(out_dir, 'KS-stats.csv'),
                              label1=summary1.label, label2=summary2.label)
    contrast1_2.write_response_stats(os.path.join(out_dir, p1_p2_contrast_name))
    if summary3 is not None:
        contrast1_3 = ProfilesContrastStats(summary1, summary3, settings)
        contrast1_3.write_response_stats(os.path.join(out_dir, p1_p3_contrast_name))

    # can only calculate the following stats if there are tracks in each profile.
    graphable = not math.isnan(contrast1_2.dur_KSd)  # can't draw graphs if there's no data
    if draw_graphs and graphable:
        labels = [s.label for s in summaries]
        colours = [s.colour for s in summaries]

        if settings is None or settings.dt_disp:
            plot_dt_distplacements(summaries, filename=out_dir + '/dt_displacements', settings=settings)

        if settings is None or settings.duration:
            plot_utils.plotCDFs([s.collate_durations() for s in summaries], names=labels,
                                xlabel='Track duration (min)',
                                filename=os.path.join(out_dir, 'CDF-duration'),
                                title='Track Duration', colours=colours)
        if settings is None or settings.length:
            plot_utils.plotCDFs([s.collate_lengths() for s in summaries], names=labels,
                                xlabel='Track length (um)',
                                filename=os.path.join(out_dir, 'CDF-length'),
                                title='Track Length', colours=colours)
        if settings is None or settings.displacement:
            plot_utils.plotCDFs([s.collate_displacements() for s in summaries], names=labels,
                                xlabel='Track displacement (um)',
                                filename=os.path.join(out_dir, 'CDF-displacement'),
                                title='Track Displacement', colours=colours)
        if settings is None or settings.displacement_angle_xy:
            plot_track_displacement_angle_xy_diffs(summaries=summaries, colours=colours,
                                filename=os.path.join(out_dir, 'displacement_angles_xy'))
            plot_utils.plotCDFs([s.collate_displacement_angles_xy() for s in summaries], names=labels,
                                xlabel='Displacement angles XY plane (radians)',
                                filename=os.path.join(out_dir, 'CDF-displacement_angles_xy'),
                                colours=colours)                                
        if settings is None or settings.meander:
            plot_utils.plotCDFs([s.collate_meandering_indices() for s in summaries], names=labels,
                                xlabel='Meandering index',
                                filename=os.path.join(out_dir, 'CDF-meander'),
                                title='Track Meandering Indexes', colours=colours,
                                xmin=0.0, xmax=1.0)
        if settings is None or settings.speeds:
            plot_utils.plotCDFs([s.collate_speeds() for s in summaries], names=labels,
                                xlabel='Translation speed (um/min)',
                                filename=os.path.join(out_dir, 'CDF-speeds'),
                                title='Translation Speeds', colours=colours)
            plot_utils.plotCDFs([s.collate_median_speeds() for s in summaries], names=labels,
                                xlabel='Median translation speed (um/min)',
                                filename=os.path.join(out_dir, 'CDF-spdMed'),
                                title='Median Track Translation Speeds',
                                colours=colours)
            plot_utils.plotCDFs([s.collate_iqr_speeds() for s in summaries], names=labels,
                                xlabel='Translation speed IRQ (um/min)',
                                filename=os.path.join(out_dir, 'CDF-spdIRQ'),
                                title='Track Translation Speed IRQs', colours=colours)
        if settings is None or settings.turns:
            plot_utils.plotCDFs([s.collate_turns() for s in summaries], names=labels,
                                xlabel='Turn speed (degrees/min)',
                                filename=os.path.join(out_dir, 'CDF-turns'),
                                title='Turn Speeds', colours=colours)
            plot_utils.plotCDFs([s.collate_median_turns() for s in summaries], names=labels,
                                xlabel='Median turn speed (degrees/min)',
                                filename=os.path.join(out_dir, 'CDF-turnsMedian'),
                                title='Median Track Turn Speeds', colours=colours)
            plot_utils.plotCDFs([s.collate_iqr_turns() for s in summaries], names=labels,
                                xlabel='Turn speed IRQ (degrees/min)',
                                filename=os.path.join(out_dir, 'CDF-turnsIRQ'),
                                title='Track Turn Speed IRQs', colours=colours)
        if settings is None or settings.rolls:
            plot_utils.plotCDFs([s.collate_rolls() for s in summaries], names=labels,
                                xlabel='Roll speed (degrees/min)',
                                filename=os.path.join(out_dir, 'CDF-rolls'),
                                title='Roll Speeds', colours=colours)
            plot_utils.plotCDFs([s.collate_median_rolls() for s in summaries], names=labels,
                                xlabel='Median roll speed (degrees/min)',
                                filename=os.path.join(out_dir, 'CDF-rollsMedian'),
                                title='Median Track Roll Speeds',
                                colours=colours)
            plot_utils.plotCDFs([s.collate_iqr_rolls() for s in summaries], names=labels,
                                xlabel='Roll speed IRQ (degrees/min)',
                                filename=os.path.join(out_dir, 'CDF-rollsIRQ'),
                                title='Track Roll Speed IRQs', colours=colours)
        if settings is None or settings.fmi:
            if len(summary1.collate_fmis()) > 0 or len(summary2.collate_fmis()) > 0:
                plot_utils.plotCDFs([s.collate_fmis() for s in summaries], names=labels,
                                    xlabel='Track forward migration index',
                                    filename=os.path.join(out_dir, 'CDF-trackFMI'),
                                    colours=colours)
                plot_utils.plotCDFs([s.collate_instant_fmis() for s in summaries], names=labels,
                                    xlabel='Instant forward migration index',
                                    filename=os.path.join(out_dir, 'CDF-instantFMI'),
                                    colours=colours)
        if settings is None or settings.msd:
            plotMSDs(summaries, filename=out_dir + '/msdDiff', settings=settings, plotRaw=False, plotFit=True)
            plot_msd_vs_max_dt_diffs(
                [profile.Profile.calculate_msd_vs_max_dt(s.profiles, upper_max_dt=settings.msd_max_dt)
                 for s in summaries],
                filename=os.path.join(out_dir, 'msd_vs_max_dt_diff.png'), names=labels)

        if settings is None or settings.disp_autocorr:
            plot_disp_autocorr_diffs(summaries=summaries, filename=out_dir + '/disp_autocorr', settings=settings)
            plot_disp_autocorr_diffs(summaries=summaries, filename=out_dir + '/disp_autocorr_normalised',
                                     settings=settings, 
                                     dataseries=lambda s: s.collate_displacement_autocorrelation()[1])

    print('finished contrasting profiles.')
    return summary1, summary2, summary3, contrast1_2, contrast1_3

