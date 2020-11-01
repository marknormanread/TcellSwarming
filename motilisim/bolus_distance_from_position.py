# Reads a simulation-generated Position file, and a simulation parameters file (to find bolus location), and creates
# a new file describing each spot's distance from the bolus.
# This can be fed into the swarming metric.
#
# Also creates a file "bolus_infiltration.csv" describing the absolute number and proportion of total cells that have
# infiltrated the tumouroid.
#
# Assumes the bolus is cylindrical through Z, hence a circle in the XY plane. As such, Z positional values 
# are ignored, as the closest distance of a cell to the bolus is at the cell's Z-axis value anyway. 
#
# Mark N. Read, 2018

import numpy as np
import os
import pandas
import shutil
import sys
import xml.etree.ElementTree as ET

# Directory in which to operate.
direc = sys.argv[sys.argv.index('-i') + 1]

outdir = direc

pos_f = direc + '/_Position.csv'  # Must contain columns named "Position_X", "Position_Y" and "Time".
param_f = direc + '/parameters.xml'  # Used to read bolus and environmental characteristics (radius).

# If data set has higher temporal resolution than desired, subsampling is possible.
# This is done at the level of unique recorded time-points.
# Because of float imprecision, absolute times or time-intervals aren't used (e.g., 95000.0012 seconds).
# Rather, every n'th recorded time point is specified.
subsample = None
if '-subsample' in sys.argv:
    subsample = int(sys.argv[sys.argv.index('-subsample') + 1])
    ss_dir = direc+'/subsample_' + str(subsample)
    os.makedirs(ss_dir, exist_ok=True)
    shutil.copy(param_f, ss_dir + '/parameters.xml')
    outdir = ss_dir


# Exception raised if these are not set to non-None, below.
bolus_x = None
bolus_y = None
bolus_r = None

# Retrieve bolus characteristics from parameters file.
tree = ET.parse(param_f)
root = tree.getroot()

node = root.find('./Environment/BoundedCylinder/bolus/radius')
if node is not None and node.text is not None:
    bolus_r = int(node.text)
    print('bolus_r = {:d}'.format(bolus_r))

node = root.find('./Environment/BoundedCylinder/radius')
if node is not None and node.text is not None:
    environment_r = float(node.text)
    bolus_x = environment_r
    bolus_y = environment_r
    print('bolus_x = {:.2f}'.format(bolus_x))

if bolus_x is None or bolus_r is None:
    raise Exception("Bolus parameters not detected. Exiting.")


# Read in positional data
print('Reading in positional file (can take a bit of time for larger files.)')
pos = pandas.read_csv(pos_f, index_col=False)
times = pos['Time'].tolist()
uniq_times = sorted(list(set(times)))
# Subsample every n'th time observed in the data set
if subsample is not None:
    # Performed by sub-sampling the range of indices for the uniq_times list. 
    ss_times = [uniq_times[i] for i in list(range(len(uniq_times))) if i % subsample == 0]
else:
    ss_times = uniq_times

ss_pos = pos.loc[pos['Time'].isin(ss_times)]  # Filter positional data for sub-sampled times.


print('Calculating distances from bolus')
# Maths in pandas is super fast. 
orig_headers = list(ss_pos.columns.values)  # Used below, to save subsampled table.
ss_pos.loc[:, 'dx'] = ss_pos['Position_X'] - bolus_x 
ss_pos.loc[:, 'dy'] = ss_pos['Position_Y'] - bolus_y
ss_pos.loc[:, 'dx2'] = ss_pos['dx'] * ss_pos['dx']  # Squared distances from bolus
ss_pos.loc[:, 'dy2'] = ss_pos['dy'] * ss_pos['dy']
ss_pos.loc[:, 'sum'] = ss_pos['dx2'] + ss_pos['dy2']
ss_pos.loc[:, 'distance_centre'] = np.sqrt(ss_pos['sum'])
ss_pos.loc[:, 'distance_bolus'] = ss_pos['distance_centre'] - bolus_r

# Save output. 
print('Saving distance to bolus as csv file')
ss_pos[['distance_bolus', 'Time']].to_csv(outdir + "/distances_to_bolus.csv", index=False)
if subsample is not None:
    # Save only the original header outputs
    ss_pos[orig_headers].to_csv(outdir+'/_Position.csv', index=False)
print('Done')


print('Calculating bolus infiltration')
ss_pos.loc[:, 'bolus_infiltration'] = ss_pos['distance_bolus'] < 0  # Negative values are inside bolus
total_at_time = []  # One item per time point
abs_infiltrates_at_time = []  # One item per time point
prop_infiltrates_at_time = []  # One item per time point
for time in ss_times:
    mask = ss_pos['Time'] == time
    ss_pos_time = ss_pos.loc[mask, :]  # Filter out current time.
    n_time = sum(mask) # How many positions (cells) recorded at the current time?
    n_bolus_infiltrates = sum(ss_pos_time['bolus_infiltration'] == 1)
    prop_bolus = n_bolus_infiltrates / n_time if n_time > 0 else 0.

    total_at_time.append(n_time)
    abs_infiltrates_at_time.append(n_bolus_infiltrates)
    prop_infiltrates_at_time.append(prop_bolus)


print('Saving bolus infiltration as csv file')
# Save as CSV
df = pandas.DataFrame.from_dict({'Time': ss_times,
                                 'Total cells': total_at_time,
                                 'Total bolus infiltrates': abs_infiltrates_at_time,
                                 'Proportion bolus infiltrates': prop_infiltrates_at_time})

df.to_csv(direc + '/bolus_infiltration.csv', index=False)
print('Done')

