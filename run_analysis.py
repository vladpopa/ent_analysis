import glob, os, sys
sys.path.append(os.getcwd() + '/lib/')
sys.path.append(os.getcwd() + '/cloudtracker/')

import multiprocessing as mp
from multiprocessing import Pool
PROC = 12

import model_param as mc
import cloudtracker.main

from conversion import *
from time_profiles import make_profiles
from id_profiles import core_profiles

def run_conversion(filelist):
	# bin3d2nc conversion
	pool = mp.Pool(PROC)
	pool.map(convert.main, filelist)
	
	# generate_tracking
	# Wrap the module for multi-processing
	pool = mp.Pool(PROC)
	pool.map(conv_wrapper, enumerate(filelist))
	
def conv_wrapper((time_step, filename)):
	generate_tracking.main(time_step, filename)
	
def run_cloudtracker():
	# Change the working directory for cloudtracker
	os.chdir('./cloudtracker')
	
	model_config = mc.model_config
	
	model_config['input_directory'] = model_config['data_directory'] + 'tracking/'
	cloudtracker.main.main(model_config)
	
	# Return to entrainment analysis directory
	os.chdir('../')	

def run_profiler(filelist):
	if(time_profiles):
		### time_profiles
		# Return to entrainment analysis directory
		os.chdir('./time_profiles')	
		
		pool = mp.Pool(PROC)
		pool.map(time_profiles_wrapper, enumerate(filelist))
		
		# Return to entrainment analysis directory
		os.chdir('../')	
	
	if(id_profiles):
		### id_profiles
		os.chdir('./id_profiles')
		core_profiles.main('core')
	
def time_profiles_wrapper((time_step, filename)):
	make_profiles.main(time_step, filename)

def main():
	filelist = glob.glob('%s/variables/*.nc' % mc.data_directory)
	filelist.sort()

	if(conversion):
		### File conversion (.bin3d -> netCDF)
		run_conversion(filelist)
	
	if(cloudtracker):
		### Cloudtracker
		run_cloudtracker()
	
	### Additional Profiles
	run_profiler(filelist)
	
if __name__ == '__main__':
	main()
	
	print 'Entrainment analysis completed'
	
