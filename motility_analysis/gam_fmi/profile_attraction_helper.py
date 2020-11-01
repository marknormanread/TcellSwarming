"""
Execute from parent directory, as:
$> python3 -m imaris.profile_positive_attraction_tiled.py

(this is needed to import motility_analysis, which sits in the root project directory)

Script characterises the positive attraction forces attracting cells to the bolus. This is done by drawing 2D density
scatter graphs of x-position, track age, direction of approach to bolus, translational speed and turn speed with
one another.

Mark N. Read, 2017
"""
from contextlib import redirect_stdout # Used to redirect gam.summary() to the file system, rather than standard out.
import io  # Used to redirect gam.summary() to the file system, rather than standard out.
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm  # color maps
import math
import motility_analysis.profile
import numpy as np
import os
import pandas
from pygam import LinearGAM, s, te
import scipy.stats
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
import sklearn.linear_model
import seaborn as sns
import sys
import motility_analysis.graphing as graphing


matplotlib.rcParams['figure.autolayout'] = True



def chemotactic_response_handler(direct, profiles, arrest_coefficient_filter):
    """ 
    Python-based heatmaps of FMI, speeds, directionality vs distance and time. 
    Outputs measures in a csv file for downstream analysis in R. 
    Profiles can be either a motility_analyis.profile.Profile object, or a list of such objects. 
    """
    # Analyse the chemotactic response of these data sets. Can also draw graphs using python infrastructure. 
    outdir = '{:s}/arrest_coeff_filter_{:.2f}'.format(direct, arrest_coefficient_filter)
    os.makedirs(outdir, exist_ok=True)
    chemotactic_response = process(profiles, outdir)
    chemotactic_response.to_csv(outdir + '/chemotactic_response.csv', index=False)


def motility_graphing_handler(direct, profiles, arrest_coefficient_filter):
    outdir = '{:s}/arrest_coeff_filter_{:.2f}'.format(direct, arrest_coefficient_filter)
    graph_dir = outdir + '/motility_summary'
    os.makedirs(graph_dir, exist_ok=True)
    graphing.plot_profile_graphs(profiles, directory=graph_dir)


def representative_tracks_handler(direct, profiles):
    xlim = None  # e.g. [-125, 125]
    ylim = None
    os.makedirs(direct, exist_ok=True)
    for i in list(range(10)):  # Representative tracks are randomly selected; sample repeatedly with different seeds. 
        graphing.plot_representative_tracks(Ps=profiles, graph_path=direct + '/representative_tracks_' + str(i),
                                            n=40, xlim=None, ylim=None,
                                            sampling="random", tick_interval=None, same_xy_lim=False, mark_origin=True,
                                            seed=i)


def build_profile_handler(data_fp, arrest_coefficient_filter, graphs=False, trim_displacement=False, fmi_direction=(-1, 0, 0)):
    prof = motility_analysis.profile.build_profile(
        directory=data_fp, graphs=False,
        trim_displacement=trim_displacement,
        trim_arrest_coefficient=arrest_coefficient_filter, tracks=None, fmi_direction=fmi_direction,
        analyse_teleports=False, interpolate=True, style="whitegrid")
    if graphs:
        motility_graphing_handler(data_fp, prof, arrest_coefficient_filter)
    return prof



def angle_between(v1, v2):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'.
    Copied on 2017-09-21 from:
    https://newtonexcelbach.wordpress.com/2014/03/01/the-angle-between-two-vectors-python-version/
    """
    cosang = np.dot(v1, v2)
    sinang = np.linalg.norm(np.cross(v1, v2))
    rads = np.arctan2(sinang, cosang)
    degs = math.degrees(rads)
    return degs


def fishers_z(r):
    """
    Transforms data from [0, 1] to [-inf, inf] which behaves better with many statistical techniques.
    Operates on a list
    """
    z = [0.5 * np.log((1. + ri) / (1. - ri)) for ri in r]
    return np.array(z)


def fishers_z_inverse(z):
    """
    Operates on a list
    """
    r = [(np.exp(2. * zi) - 1) / (np.exp(2. * zi) + 1) for zi in z]
    return np.array(r)


class BinnedHeatmap():
    def __init__(self, xs, times, response, resolution=5):
        empirical_x_max = np.max(xs)
        empirical_x_min = np.min(xs)
        empirical_time_max = np.max(times)
        # calculate bin boundaries
        self.x_bins = np.linspace(empirical_x_min, empirical_x_max+1, resolution+1)
        self.y_bins = np.linspace(0., empirical_time_max+1, resolution+1)

        # creates a list of list of lists. Bear with me... first list represents bins for x, and contains lists
        # representing bins for y. The last hierarchy of lists is what we add data to, and they can be of variable
        # length. I  would have preferred to do this using numpy.ndarray, but couldn't find a way. ndarray data types
        # need to be of fixed length, and this isn't.
        grid = [[ [] for _ in list(range(resolution)) ] for _ in list(range(resolution))]
        # sort the response values into the bins
        # note that bin ids start from 1. Need to subtract 1 if bin data is stored in array/list
        xis = np.digitize(xs, self.x_bins)
        yis = np.digitize(times, self.y_bins)
        for xi, yi, r, xx, yy in zip(xis, yis, response, xs, times):
            grid[xi-1][yi-1].append(r)  # -1 because

        # calculate the mean value in each bin, store as 2D ndarray
        self.grid_medians = np.ndarray(shape=(resolution, resolution), dtype=float)
        for xi in list(range(len(grid))):
            for yi in list(range(len(grid[0]))):
                self.grid_medians[xi, yi] = np.median(grid[xi][yi])

    def predict(self, predict_points):
        z = []
        # note that bin ids start from 1. Need to subtract 1 if bin data is stored in array/list
        xis = np.digitize(predict_points[:, 0], self.x_bins)
        yis = np.digitize(predict_points[:, 1], self.y_bins)
        for xi, yi, x, y in zip(xis, yis, predict_points[:, 0], predict_points[:, 1]):
            # print('xi = {}, yi = {}, x = {}, y = {}'.format(xi, yi, x, y))
            z.append(self.grid_medians[xi-1, yi-1])
        return np.array(z)


def three_d_plots(response, xs, times, method, title, graph_fp, z_transform=False, contour_plot=True, n_splines=5,
                  cmap=cm.rainbow, centred_val=None):
    """
    Try to separate effects of distance from bolus and age of experiment on the given response. This involves
    fitting a 2d plane in 3d space (x=distance, y=time, z=response).

    Three methods are implemented.
    """
    empirical_x_max = np.max(xs)
    empirical_x_min = np.min(xs)
    empirical_time_max = np.max(times)
    if z_transform:
        response = fishers_z(response)  # operates over entire list
    X = np.column_stack((xs, times))  # combine data for input to regressors
    xx, yy = np.meshgrid(np.arange(empirical_x_min, empirical_x_max, 5.),
                         np.arange(0., empirical_time_max, 1.))
    # np_c code here return a list (representing samples) containing coordinates (representing features).
    predict_points = np.c_[xx.ravel(), yy.ravel()]
    r_score = ''
    # Linear regression
    if method is 'lr':
        lr = sklearn.linear_model.LinearRegression(fit_intercept=True, normalize=True)
        lr.fit(X, response)
        Z = lr.predict(predict_points)
        r_score = '{:.2f}'.format(lr.score(X, response))

    if method is 'lr2':  # quadratic linear regression
        data = np.c_[xs, times, response]
        # best-fit quadratic curve
        A = np.c_[np.ones(data.shape[0]), data[:, :2], np.prod(data[:, :2], axis=1), data[:, :2] ** 2]
        C, _, _, _ = scipy.linalg.lstsq(A, data[:, 2])
        XX = xx.flatten()
        YY = yy.flatten()
        # evaluate it on a grid
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX * YY, XX ** 2, YY ** 2], C).reshape(xx.shape)

    # Support vector regression
    elif method is 'svr':
        scaler = StandardScaler()
        X_ = scaler.fit_transform(X)  # scale features to zero mean unit variance
        svr = SVR(C=1., epsilon=0.2, kernel='poly', degree=5)
        svr.fit(X_, response)
        r_score = '{:.2f}'.format(svr.score(X, response))
        print("score = " + str(svr.score(X, response)))  # how good a fit is this?

        predict_points_ = scaler.transform(predict_points)  # scale to zero mean unit variance
        Z = svr.predict(predict_points_)

    # Generalized Additive Model
    elif method is 'gam':
        # combine data for input to regressors
        # set up features for interaction effets: x, y, x+y, x-y.
        X = np.column_stack((xs, times))
        # gam = LinearGAM(n_splines=n_splines).gridsearch(X, response)
        gam = LinearGAM(te(0, 1)).fit(X, response)        
        gam.summary()  # Prints to std out. 

        with io.StringIO() as buf, redirect_stdout(buf):
            gam.summary()
            output = buf.getvalue()
            with open(graph_fp + '.txt', 'w') as f:  # Write to file system.
                f.write(output)

        # set up features for interaction effets: x, y, x+y, x-y.
        predict_points = np.c_[predict_points[:,0], predict_points[:,1]]
        Z = gam.predict(predict_points)

    elif method is 'heatmap':
        hm = BinnedHeatmap(xs, times, response, resolution=5)
        Z = hm.predict(predict_points)
        print(predict_points)

    if z_transform:  # transform back
        Z = fishers_z_inverse(Z)

    vmin = None  # Ranges for plotting centred values with given colormap.
    vmax = None
    if centred_val is not None:
        # Either minimum or maximum Z value is furthest away from centred_val (not sure how this behaves if you're daft
        # and choose a centred value outside of the range [min(Z), max(Z)]). Set vmax and vmin to be equally far 
        # from centred_val, using biggest difference. 
        min_diff_z = math.fabs(centred_val - np.min(Z))
        max_diff_z = math.fabs(np.max(Z) - centred_val)
        bleed = min_diff_z if min_diff_z > max_diff_z else max_diff_z
        vmin = centred_val - bleed
        vmax = centred_val + bleed

    zz = Z.reshape(xx.shape)  # back to meshgrid coordinate format
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    p = plt.pcolor(xx, yy, zz, cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(p)
    if contour_plot:
        p_cont = plt.contour(xx, yy, zz, colors='k')
        plt.clabel(p_cont, fontsize=9, inline=1)  # add labels to the contour lines
    ax.set_xlabel('Distance from bolus (um)')
    ax.set_ylabel('Experiment time (min)')
    ax.grid()
    fig.savefig(graph_fp + '.png', dpi=300, bbox_inches='tight', transparent=True)
    fig.savefig(graph_fp + '.eps', bbox_inches='tight', transparent=True)


def multi_factorial_color_scatter(response, xs, times, title, graph_fp):
    """ Plot response as a coloured dot at the distance from bolus and experimental time. """
    plt.clf()
    cmap = cm.rainbow
    plt.scatter(xs, times, c=response, cmap=cmap)
    plt.title(title)
    plt.savefig(graph_fp, dpi=300, bbox_inches='tight', transparent=True)


def process(all_prof, outdir, speed_filter=2., x_filter=0.):
    if isinstance(all_prof, motility_analysis.profile.Profile):
        all_prof = [all_prof]
    # Extract parameters for plotting from the supplied profile(s).
    data = []  # tuple (x-location, speed, turn, displacement angle-to-bolus, time(min), forward migration index)
    for p in all_prof:
        for t in p.tracks:
            for pos_id, position in enumerate(t.positions[1:]):  # ignore first position, as no speed assigned to it
                # if speed unassigned, ignore; avoid those cells that actually stick to bolus, filter by distance
                if position.speed and position.speed >= speed_filter and position.turn \
                        and position.displacementVector and position.x > x_filter:
                    angle_to_bolus = angle_between([-1, 0, 0], position.displacementVector)
                    data.append(  # append tuple
                        (position.x,
                         position.speed,
                         position.turn,
                         angle_to_bolus,
                         position.time_m,  # absolute time w/r/t start of experiment
                         pos_id,  # relative time (id) w/r/t when track observation started
                         position.instant_fmi,
                         t.fmi
                         )
                        )
    xlocs = [t[0] for t in data]
    speeds = [t[1] for t in data]
    turns = [t[2] for t in data]
    angles = [t[3] for t in data]
    abs_times = [t[4] for t in data]  # absolute time since start of experiment
    track_obs_num = [t[5] for t in data]  # time relative to when track was first recorded, in observation number
    instant_fmis = [t[6] for t in data]
    track_fmis = [t[7] for t in data]

    # compile track forward migration index values. Record the FMI at the track starting position
    fmi_t = []  # time, distance from bolus and FMI are all recorded.
    fmi_x = []
    fmi = []
    for p in all_prof:
        for t in p.tracks:
            track_fmi = t.fmi
            fmi.append(track_fmi)
            fmi_t.append(t.positions[-1].time_m)
            fmi_x.append(t.positions[-1].x)  # end position

    # Compile data for returning (and possible writing to FS).
    chemotactic_response = pandas.DataFrame(data={'dist_bolus': xlocs,
                                                  'absolute_time': abs_times,
                                                  # 'track_id',
                                                  # 'fmi': ,
                                                  'instant_fmi': instant_fmis,
                                                  'track_fmi': track_fmis,
                                                  'angle_bolus': angles,
                                                  'speed': speeds,
                                                  'turn_speed': turns})
    
    # Corresponds to Jack's f metric - velocity w.r.t. the radius from the bolus; negative values = attraction to bolus
    chemotactic_response['f'] = chemotactic_response['speed'] * chemotactic_response['instant_fmi']
    f = chemotactic_response['f'].values.tolist()

    # Draw the graphs
    if True:
        multi_factorial_color_scatter(response=fmi, xs=fmi_x, times=fmi_t, title='Forward migration index',
                                      graph_fp=outdir + '/spatiotemporal_fmi')

        three_d_plots(response=fmi, xs=fmi_x, times=fmi_t, method='heatmap', title='Forward migration index',
                      graph_fp=outdir + '/hm_distance_vs_exp_time_vs_fmi',
                      z_transform=False, contour_plot=False, cmap=cm.seismic, centred_val=0.0)
        three_d_plots(response=fmi, xs=fmi_x, times=fmi_t, method='gam', title='Forward migration index',
                      graph_fp=outdir + '/gam_distance_vs_exp_time_vs_fmi', z_transform=False,
                      cmap=cm.seismic, centred_val=0.0)


        three_d_plots(response=instant_fmis, xs=xlocs, times=abs_times, method='gam', title='Instantaneous FMI',
                      graph_fp=outdir + '/gam_distance_vs_exp_time_vs_instant_fmi', z_transform=False, n_splines=4,
                      cmap=cm.seismic, centred_val=0.0)
        three_d_plots(response=instant_fmis, xs=xlocs, times=abs_times, method='heatmap', title='Instantaneous FMI',
                      graph_fp=outdir + '/hm_distance_vs_exp_time_vs_instant_fmi', z_transform=False, contour_plot=False,
                      cmap=cm.seismic, centred_val=0.0)

        three_d_plots(response=f, xs=xlocs, times=abs_times, method='heatmap', title='F',
                      graph_fp=outdir + '/hm_distance_vs_exp_time_vs_f',
                      z_transform=False, contour_plot=False, cmap=cm.seismic, centred_val=0.0)
        three_d_plots(response=f, xs=xlocs, times=abs_times, method='gam', title='F',
                      graph_fp=outdir + '/gam_distance_vs_exp_time_vs_f', z_transform=False,
                      cmap=cm.seismic, centred_val=0.0)

        three_d_plots(response=angles, xs=xlocs, times=abs_times, method='heatmap', title='Angle of approach to bolus (degrees)',
                      graph_fp=outdir + '/hm_distance_vs_exp_time_vs_taxis', z_transform=False, contour_plot=False,
                      cmap=cm.seismic, centred_val=90.0)
        three_d_plots(response=angles, xs=xlocs, times=abs_times, method='gam', title='Angle of approach to bolus (degrees)',
                      graph_fp=outdir + '/gam_distance_vs_exp_time_vs_taxis',
                      cmap=cm.seismic, centred_val=90.0)

        three_d_plots(response=speeds, xs=xlocs, times=abs_times, method='heatmap', title='Translational speeds (um/min)',
                      graph_fp=outdir + '/hm_distance_vs_exp_time_vs_speeds', z_transform=False, contour_plot=False,
                      cmap=cm.plasma)
        three_d_plots(response=speeds, xs=xlocs, times=abs_times, method='gam', title='Translational speeds (um/min)',
                      graph_fp=outdir + '/gam_distance_vs_exp_time_vs_speeds',
                      cmap=cm.plasma)

        three_d_plots(response=turns, xs=xlocs, times=abs_times, method='heatmap', title='Turn speeds (degrees/min)',
                      graph_fp=outdir + '/hm_distance_vs_exp_time_vs_turns', z_transform=False, contour_plot=False,
                      cmap=cm.inferno)
        three_d_plots(response=turns, xs=xlocs, times=abs_times, method='gam', title='Turn speeds (degrees/min)',
                      graph_fp=outdir + '/gam_distance_vs_exp_time_vs_turns',
                      cmap=cm.inferno)


    return chemotactic_response
