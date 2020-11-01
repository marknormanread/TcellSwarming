# Reenacts a simulation, saving each time frame as a graph (matplotlib) that can later be stitched into a movie. 
# This is much faster than using the simulation's own GUI, and does not require an X window. 
#
# Also plotted are the number of captured cells per time point, and the mean X, Y and Z positions of the cell
# population at each time point. 
#
# Mark N. Read, 2018
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import os
import pandas
import seaborn as sns
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

# Essential to supply this
directory = sys.argv[1]
stills_directory = directory + '/stills'

mkmovie_fp = '/Users/markread/dropbox_usyd/projects/biro/neutroswarm/motilisim/makemovie_supply_dir.sh'

shutil.rmtree(stills_directory, ignore_errors=True)
os.makedirs(stills_directory, exist_ok=True)

pos_df = pandas.read_csv(directory + '/_Position.csv')
pos_df.dropna(axis='index', how='any', inplace=True)

# Pull out unique time points
time_points = pos_df['Time'].unique()  # Minutes
time_points.sort()
time_points_h = time_points / 60.  # Hours, better for plotting

# Parameters for plotting
rw_colour = 'limegreen'
chemotactic_colour = 'blue'
bolus_infiltration_colour = 'orangered'
bolus_colour = 'lightgrey'
cell_size = 1.

# Retrieve bolus characteristics from parameters file.
param_fp = directory + '/parameters.xml'  # Used to read bolus and environmental characteristics (radius).
tree = ET.parse(param_fp)
root = tree.getroot()
node = root.find('./Environment/BoundedCylinder/bolus/radius')
if node is not None and node.text is not None:
    bolus_r = int(node.text)

node = root.find('./Environment/BoundedCylinder/radius')
if node is not None and node.text is not None:
    environment_r = float(node.text)
    bolus_x = environment_r
    bolus_y = environment_r

if bolus_x is None or bolus_r is None:
    raise Exception("Bolus parameters not detected. Exiting.")
pos_df.loc[:, 'dx'] = pos_df['Position_X'] - bolus_x 
pos_df.loc[:, 'dy'] = pos_df['Position_Y'] - bolus_y
pos_df.loc[:, 'dx2'] = pos_df['dx'] * pos_df['dx']  # Squared distances from bolus
pos_df.loc[:, 'dy2'] = pos_df['dy'] * pos_df['dy']
pos_df.loc[:, 'sum'] = pos_df['dx2'] + pos_df['dy2']
pos_df.loc[:, 'distance_centre'] = np.sqrt(pos_df['sum'])
pos_df.loc[:, 'distance_bolus'] = pos_df['distance_centre'] - bolus_r


# Cycle through each time point
for frame_id, time in enumerate(time_points):
    print('Processing time ' + str(time))
    mask = pos_df['Time'] == time
    # print(mask)
    frame_positions = pos_df.loc[mask, :]
    
    # Colour cells by state. 
    frame_positions.loc[:, 'Colour'] = rw_colour
    mask_chemotactic = frame_positions['Chemotactic'] == True    
    frame_positions.loc[mask_chemotactic, 'Colour'] = chemotactic_colour
    mask_inside_bolus = frame_positions['distance_bolus'] <= 0
    frame_positions.loc[mask_inside_bolus, 'Colour'] = bolus_infiltration_colour
    
    fig = plt.figure(figsize=[5,5], facecolor='white', clear=True)
    ax = fig.gca()

    plt.scatter(frame_positions['Position_X'], frame_positions['Position_Y'], c=frame_positions['Colour'], s=cell_size)

    # Colour the bolus 
    bolus = Circle(xy=[bolus_x, bolus_y], radius=bolus_r, fill=True, facecolor=bolus_colour, zorder=0)
    ax.add_artist(bolus)

    # Print the time on the bottom left corner.
    h, m = divmod(time, 60)
    plt.text(0.05 * environment_r, 0.05 * environment_r, '{h:02.0f}:{m:02.0f}'.format(h=abs(h), m=m), fontsize=16)

    # Turn off tick marks
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False) # labels along the bottom edge are off
    
    # Force dimensions plotted
    plt.xlim([0, environment_r * 2])
    plt.ylim([0, environment_r * 2])
    plt.savefig(stills_directory + '/img_{id:03d}'.format(id=frame_id), dpi=200)
    plt.close(fig)

print(directory)
subprocess.run([mkmovie_fp, directory])


captured_cells = []  # How many cells in each time-point
# X, Y and Z locations of cells at each time point. 
captured_cells_xs = []
captured_cells_ys = []
captured_cells_zs = []
mean_xs = []
mean_ys = []
mean_zs = []
# Cycle through each time point
for frame_id, time in enumerate(time_points):
    mask = pos_df['Time'] == time
    frame_positions = pos_df.loc[mask, :]
    cells_at_timepoint = len(frame_positions)
    captured_cells.append(cells_at_timepoint)
    mean_x = np.mean(frame_positions['Position_X'])
    mean_y = np.mean(frame_positions['Position_Y'])
    mean_z = np.mean(frame_positions['Position_Z'])
    mean_xs.append(mean_x)
    mean_ys.append(mean_y)
    mean_zs.append(mean_z)


sns.set_style('darkgrid')
plt.clf()
sns.lineplot(x=time_points_h, y=captured_cells)
plt.xlabel("Time (h)")
plt.ylabel("Cells")
plt.savefig(directory + '/total_captured_cells.png', bbox_inches='tight', dpi=300)

plt.clf()
sns.lineplot(x=time_points_h, y=mean_xs)
plt.xlabel("Time (h)")
plt.ylabel("Mean cell X location (um)")
plt.savefig(directory + '/cell_mean_xs.png', bbox_inches='tight', dpi=300)

plt.clf()
sns.lineplot(x=time_points_h, y=mean_ys)
plt.xlabel("Time (h)")
plt.ylabel("Mean cell Y location (um)")
plt.savefig(directory + '/cell_mean_ys.png', bbox_inches='tight', dpi=300)

plt.clf()
sns.lineplot(x=time_points_h, y=mean_zs)
plt.xlabel("Time (h)")
plt.ylabel("Mean cell Z location (um)")
plt.savefig(directory + '/cell_mean_zs.png', bbox_inches='tight', dpi=300)