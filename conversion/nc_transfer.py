#!/usr/bin/python                                                                                             
import os, glob, shutil
import ent_analysis.lib.model_param as mc

def transfer(filename):
	path, name = os.path.split(filename)
	f = name.split('_')
	print "Transferring... " + name
	
	# netCDF file destination base folder
	dst = mc.data_directory
	if('CORE' in f):
		shutil.move(filename, '%s/core_entrain/' % (dst))
	elif ('CLOUD' in f):
		shutil.move(filename, '%s/condensed_entrain/' % (dst))
	else:
		shutil.move(filename, '%s/variables/' % (dst))

	return

if __name__ == "__main__":
	transfer()