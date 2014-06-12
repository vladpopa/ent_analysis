# before running: need to modify bomex.py and model.cfg according to the case
# (bomex.py is used for plume_statistics.py and should contain the same model configuration as model.cfg,
# which cloudtracker uses)

# this runs cloudtracker and produces profiles for each id'ed plume
# (it runs cloudtracker in the current directory, which is currently configured to track plumes)
# note: the plume tracker still needs to be tested (try a dry case)
execfile('run_analysis.py')

# this outputs a .pkl file that contains a dictionary of sub-cloud plume statistics 
# (lifetime, cross-sectional area, and time of condensation)
execfile('plume_statistics.py')

# plot histogram of cross-sectional area of the sub-cloud plumes
execfile('plot_area.py')
