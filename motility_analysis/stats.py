# Module contains useful statistical tools and functions.
#
#
import matplotlib
matplotlib.use('Agg')   # Force matplotlib to not use any Xwindows backend. Must be called before any pyplot import.
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


def AtestSlow(dist1, dist2):
    """
    Implementation of Vargha-Delaney's A test, a non-parametric effect magnitude test. The distributions may be of
    different sizes in this implementation.

    This is not an efficient implementation for distributions with large numbers of samples, but its operation is
    intuitive.
    """
    equal = 0.0
    greater = 0.0
    for x in dist1:
        for y in dist2:
            if x == y:
                equal += 1
            elif x > y:
                greater += 1

    # multiplication of the number of samples in each distribution
    nm = len(dist1) * len(dist2)
    return (greater / nm) + ((0.5 * equal) / nm)


def Atest(dist1, dist2):
    """
    Fast implementation of Vargha-Delaney's A test, a non-parametric effect magnitude. Value between 0 and 1 returned,
    which represents the probability of a randomly selected sample from population A being bigger than a randomly
    selected sample from population B.

    In this implementation, the sizes of the two input distributions may differ. This is a much faster, but less
    intuitive, implementation than AtestSlow above. Here the two distributions are sorted into monotonically increasing
    values. As such, the location and number of samples in one distribution that are less than a member of the other
    distribution is known, and hence they do not need to be frequently revisited.
    """
    equal = 0.0
    greater = 0.0
    A = sorted(dist1)
    B = sorted(dist2)
    # where to start searching through B. All B[0..start-1] < a. Either B[start] > a or B[start]==a (and may the the
    # first in a string of B == a.
    start = 0
    for a in A:
        # don't start from the beginning. Values monotonically increase in both lists.
        for i in range(start, len(B)):
            if B[i] < a:		# this value smaller than a, so can start from the next value next time.
                start = i+1
            if B[i] == a:
                equal += 1
            if B[i] > a:		# no point looking further, B[i:end]>a.
                break
        greater += start
    # multiplication of the number of samples in each distribution
    nm = len(dist1) * len(dist2)
    return (greater / nm) + ((0.5 * equal) / nm)


def ecdfPlot(dist):
    """
    Calculates the X and Y coorindates for an empirical cumulative
    distribution plot. It is surprisingly difficult to find a
    good implementation (easy to install and use), so I've written
    my own.

    'dist' contains the empirically observed values. They do not
    need to be sorted.

    Returns a tuple, the X coordinates and the 'cdf' which contains
    monotonically values between 0 to 1. Both lists are of the same
    length, and can be used to plot the graph.
    """
    X = sorted(dist)	# X coorinates
    cdf = []			# will put the corresponding CDF values in here.
    n = len(dist)
    for xi in X:		# scan through each x-coordinate
        smallEq = 0.0	# count how many values in the distribution are smaller
        for j in dist:
            if xi >= j : smallEq += 1
        prop = smallEq / n # calculate proportion
        cdf.append(prop)
    return (X, cdf)


def ecdf(observations):
    """
    Empirical Cumulative Distribution Function.
    arguments:
      :observations : list of observations, in any order.

    Returns a list of tuples: [ (<value>,<proportion of observations less than or equal to value>) ]
    Values will be given in ascending order.
    """
    s = sorted(observations)	# non-destructive, preserves original order of observations.
    n = len(s)					# number of observations.
    cdf = []
    for i, sample in enumerate(s):
        proportion = float(i+1)/n
        tup = (sample, proportion)
        cdf.append(tup)
    return cdf


def prop_ecdf(ecdf, value, start=0):
    """
    Returns the proportion of observations in a sample at or below the specified value.
    arguments:
      :ecdf : array of tuple values, [(<value>,<proportion>)], sorted in ascending order of value, where prop gives the
              proportion of all observations in the sample at or below value. A structure of this sort is returned by
              this module's ecdf function.
      :value : value to be interrogated.
      :start : index in ecdf to start looking for a match. Used for making KS significantly more efficient.
    """
    VALUE = 0
    PROP = 1
    if not ecdf:
        raise Exception('No observations supplied.')
    proportion = 0.0
    # index is where to start looking next time. This is the index where the current match was found.
    for index in range(start, len(ecdf)):
        record = ecdf[index]
        if value >= record[VALUE]:
            proportion = record[PROP]
        else:
            break
    # by definition, index was the last record OVER the value, hence index - 1
    return proportion, max(0,index-1)	# ensure index is not less than zero.


def ks_2samp(dist1, dist2):
    """
    Calculates the Kolgomorov-Smirnov 2 sample test between the two supplied distributions. Distributions do not have
    to be in order.

    This implementation is not as fast as that of python's Scipy. I've tried to write this to be intuitively appealing.
    """
    VALUE = 0		# indexes of values and proportions of distribution lying below that value, in the tuple returned
    PROP = 1		# by ecdf.
    DISTRO = 2		# third item in tuple is boolean, stores if record came from ecdf1 or ecdf2.
    ecdf1 = ecdf(dist1)
    ecdf2 = ecdf(dist2)
    # tuples contain following: (<value>,<T/F>,<prop>). Where T/F is true when the value came from distro 1,
    # and false if it came from distro 2. Hence, <prop> corresponds to distro 1 if <T/F> is true, and vice versa.
    values = []
    for record in ecdf1:
        values.append( (record[VALUE], record[PROP], True) )
    for record in ecdf2:
        values.append( (record[VALUE], record[PROP], False) )
    # sort the array of tuples by the first value in each tuple.
    values.sort(key=lambda tup: tup[0])
    d = 0.0
    # where to start searching through ecdf records. Do not need to start from the beginning, since values are
    # monotonically increasing. Simply need to store how far through the records the previous search went, as this is
    # where to start for the next search.
    start1 = 0
    start2 = 0
    for record in values:
        if record[DISTRO]:
            p, start1 = prop_ecdf(ecdf2, record[VALUE], start1)
            i = abs(record[PROP] - p)
        else:
            p, start2 = prop_ecdf(ecdf1, record[VALUE], start2)
            i = abs(p - record[PROP])
        if i > d:
            d = i		# store new D value, if highest seen so far.
    return d


#
# #
# def main():
# 	"""
# 	Used for testing purposes.
# 	"""
# 	import random
# 	import time
#
# 	A = [random.random() for _ in range(100000)]
# 	B = [random.random() for _ in range(100000)]
#
# # 	A = [1,2,4,5,6,7,8,9,11]
# # 	B = [5,6]
# 	start = time.time()
# 	result = Atest(A,B)
# 	dur = str(time.time() - start)
# 	print 'new A = ' + str(result) + '; took ' + dur
#
#
# 	start = time.time()
# 	result = AtestSlow(A, B)
# 	dur = str(time.time() - start)
# 	print 'old A = ' + str(result) + '; took ' + dur
#
#
# # python main method hook.
# if __name__ == '__main__':
# 	main()
