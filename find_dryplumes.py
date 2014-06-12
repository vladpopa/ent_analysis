import cPickle
import glob
import bomex as mc
import os
import numpy as np

plume_input_files = os.path.join(mc.analysis_dir, 'cloudtracker/pkl', \
            'plume_*.pkl')
plume_file_list = glob.glob(plume_input_files)


ids = np.array([])

print 'finding max id'
for t, filename in enumerate(plume_file_list):
   plumes = cPickle.load(open(filename, 'rb'))
   ids = np.append(ids, plumes.keys())
   
max_id = int(np.max(ids))

#max_id = 2841
#keep track of plumes that have condensed points
condensed = np.ones(max_id+1) < 0
lifetimes = np.zeros(max_id+1)

print 'finding dry plumes'
for t, filename in enumerate(plume_file_list):
    print 't = {}'.format(t)
    plumes = cPickle.load(open(filename, 'rb'))
    for id, plume in plumes.iteritems():
        lifetimes[id] = lifetimes[id] + 1 
        #if the plume has condensed points flag it
        if len(plume['condensed']) and not condensed[id]:
            #print(id, 'condensed')
            condensed[id] = True

dry_plume_ids = np.where(~condensed)[0]
condensed_plume_ids = np.where(condensed)[0]


print 'dry plume ids are: {}'.format(dry_plume_ids) 
print 'number of dry plumes: {}'.format(len(dry_plume_ids))
print 'number of plumes that form clouds: {}'.format(len(condensed_plume_ids))
print 'total number of plumes: {}'.format(max_id)

print 'number of clusters that classify as noise: {}'.format(len(lifetimes[lifetimes < 2]))

print 'average lifetime of dry plumes: {}'.format(np.mean(lifetimes[~condensed]))
print 'average lifetime of plumes that saturate: {}'.format(np.mean(lifetimes[condensed]))


   
    

