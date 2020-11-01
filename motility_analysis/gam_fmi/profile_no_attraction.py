"""
Execute from parent directory, as:
$> python3 -m imaris.profile_no_attraction.py

(this is needed to import motility_analysis, which sits in the root project directory)

Mark N. Read, 2017
"""
import motility_analysis.profile as profile
import motility_analysis.graphing as graphing

arrest_coefficient_filter = 0.8  # retain only tracks spending <80% of their time going slow

d1 = 'imaris/no_attraction/CTLsAlone_Feb2016/Imimage_Jorgeg_20141004_Colgel_Vehicle_Isotype_Position1'
d2 = 'imaris/no_attraction/CTLsAlone_Feb2016/Imimage_Jorgeg_20141004_Colgel_Vehicle_Isotype_Position2'
d3 = 'imaris/no_attraction/CTLsAlone_Feb2016/Imimage_Jorgeg_20150809_CTLs_Colgel_Vehicle_Isotype_Position1'
d4 = 'imaris/no_attraction/CTLsAlone_Feb2016/Imimage_Jorgeg_20150809_CTLs_Colgel_Vehicle_Isotype_Position2'
d5 = 'imaris/no_attraction/CTLsAlone_Feb2016/Imimage_Jorgeg_20150823_CTLs_Colgel_Vehicle_Isotype_Position1'
d6 = 'imaris/no_attraction/CTLsAlone_Feb2016/Imimage_Jorgeg_20150823_CTLs_Colgel_Vehicle_Isotype_Position2'

graphs = False
fmi_direction = (-1, 0, 0)
p1 = profile.build_profile(directory=d1, graphs=graphs, trim_duration=False,
                           trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
                           analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
p2 = profile.build_profile(directory=d2, graphs=graphs, trim_duration=False,
                           trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
                           analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
p3 = profile.build_profile(directory=d3, graphs=graphs, trim_duration=False,
                           trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
                           analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
p4 = profile.build_profile(directory=d4, graphs=graphs, trim_duration=False,
                           trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
                           analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
p5 = profile.build_profile(directory=d5, graphs=graphs, trim_duration=False,
                           trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
                           analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)
p6 = profile.build_profile(directory=d6, graphs=graphs, trim_duration=False,
                           trim_arrest_coefficient=arrest_coefficient_filter, tracks=None,
                           analyse_teleports=False, interpolate=True, timestep_s=None, fmi_direction=fmi_direction)

all_prof = [p1, p2, p3, p4, p5, p6]

graphing.plot_msd_vs_max_dt(all_prof, directory='imaris/no_attraction/CTLsAlone_Feb2016/summary')
graphing.plot_profile_graphs(all_prof, directory='imaris/no_attraction/CTLsAlone_Feb2016/summary')

