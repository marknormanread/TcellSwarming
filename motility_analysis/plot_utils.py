"""
Collects together convenient plotting routines.
"""
import matplotlib.pyplot as plt
import math
import numpy


def plot_histogram(data, filename, xlabel, ylabel, xLim=None, title=None):
    """ Generic method for plotting a histogram. """
    if not data:
        # list is empty!
        print('WARNING: cannot create plot ' + filename + ' as no data was supplied to plot.')
        return
    data = [d for d in data if not math.isnan(d)]   # remove 'nan' from data.
    plt.clf()
    plt.hist(data)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if xLim:
        plt.xlim(xLim)
    if title:
        plt.title(title)
    font = {'size': 18}
    plt.rc('font', **font)
    plt.savefig(filename + '.png', dpi=300)
    plt.savefig(filename + '.eps', dpi=300)


def ecdf(observations):
    """
    Empirical Cumulative Distribution Function.
    :param observations: list of observations, in any order.
    :return: a list of tuples: [ (<value>,<proportion of observations less than or equal to value>) ]
    Values will be given in ascending order.
    """
    s = sorted(observations)  # non-destructive, preserves original order of observations.
    n = len(s)  # number of observations.
    cdf = []
    for i, sample in enumerate(s):
        proportion = float(i + 1) / n
        tup = (sample, proportion)
        cdf.append(tup)
    return cdf


def plotCDFs(distributions, names=None, colours=None, xlabel=None, title=None, filename=None, fontsize=18,
             linewidth=2.0, xlog=False, xmin=None, xmax=None, line_styles=None, x_symlog_thresh=None,
             plot_uniform=False):
    """
    Calculates cumulative distribution functions for each distribution in the supplied list.
    :param distributions: list of lists. Each lower level list is a distribution to plot as empirical cdf
    :param names: name of each distribution, used for legend. List of strings.
    :param colours: list of objects interpretable as colours. If used, one per distribution. E.g. ['r','k','b','g']
    :param xlabel: xlabel of graph
    :param title: title of graph
    :param filename: if left as None, will show graph visually rather than save to filesystem. Leave off the extension
    (.png and similar), as these are added. Files saved as both eps and png.
    :param fontsize:
    :param linewidth:
    :param xlog: False, linear scale. True, log scale. 'symlog', symetrical log scale that has a section of linear
    around zero.
    :param x_symlog_thresh: Thresholds around zero for linear plot if 'symlog' selected.
    :param plot_uniform: if True, plot the distribution for a uniform distribution within the max and min values of the
    other distributions shown.
    """
    plt.clf()
    plots = []
    if line_styles is None:
        line_styles = ['-'] * len(distributions)
    legend_names = []
    if names is None:
        names = [None] * len(distributions)
    for i, (distro, name, ls) in enumerate(zip(distributions, names, line_styles)):
        if len(distro) > 0:
            cdf = ecdf(distro)
            # insert one record at start of these lists to ensure curves extend to y=0.
            D1X = [cdf[0][0]]
            D1X.extend([record[0] for record in cdf])
            D1Y = [0]
            D1Y.extend([record[1] for record in cdf])

            if colours:
                p1, = plt.plot(D1X, D1Y, linewidth=linewidth, color=colours[i])
            else:
                p1, = plt.plot(D1X, D1Y, linewidth=linewidth)
            plots.append(p1)
            legend_names.append(name)
    if plot_uniform:
        all_data = numpy.concatenate([numpy.array(d) for d in distributions])
        x_min = numpy.min(all_data)
        x_max = numpy.max(all_data)
        plt.plot((x_min, x_max), (0., 1.), linewidth=linewidth)
    if len(legend_names) > 0 and None not in legend_names:
        plt.legend(plots, legend_names, loc='lower right')
    plt.ylabel('Proportion')
    if xlabel:
        plt.xlabel(xlabel)
    if title:
        plt.title(title)
    plt.gca().grid(True, linewidth=linewidth)   # turn on grid lines.
    font = {'size': fontsize}
    plt.rc('font', **font)
    # change the width of the plot boundary lines
    ax = plt.gca()
    [i.set_linewidth(linewidth) for i in ax.spines.values()]
    if xmax is not None or xmin is not None:
        curXmin, curXmax = plt.gca().get_xlim()
        if xmax is not None:
            curXmax = xmax
        if xmin is not None:
            curXmin = xmin
        plt.gca().set_xlim([curXmin, curXmax])
    if xlog is True:
        plt.xscale('log')
    if xlog is 'symlog':
        if x_symlog_thresh:
            plt.xscale('symlog', linthreshx=x_symlog_thresh)
        else:
            plt.xscale('symlog')
    if filename:
        plt.savefig(filename + '.png', dpi=600)
        plt.savefig(filename + '.eps', dpi=600)
    else:
        plt.show()
