#!/usr/bin/env python

from pylab import *
import cPickle

data = cPickle.load(open('/tera/vpopa/bomex/analysis/cloudtracker/pkl/graph_events.pkl', 'rb'))

print "Number of Clouds: ", len(data)

start_with_split=0
start_with_break=0
start_new=0
end_with_break=0
end_with_merge=0
start_and_end_broken=0
start_and_end_split_merge=0
core_clouds = 0
condensed_clouds = 0
clean_end = 0

complete_clouds = []


for i in data:
    complete = False
    times = data[i].keys()
    if data[i]['has_core']:
        core_clouds = core_clouds + 1
    times.remove('has_core')
    if data[i]['has_condensed']:
        condensed_clouds = condensed_clouds + 1
    times.remove('has_condensed')
    times.sort()
    
    if 'split_start' in data[i][times[0]]:
        start_with_split = start_with_split + 1
        if 'merge_end' in data[i][times[-1]]:
            start_and_end_split_merge = start_and_end_split_merge + 1
    elif 'break_start' in data[i][times[0]]:
        start_with_break = start_with_break + 1
    else:
        start_new = start_new + 1
        complete = True

    if 'break_end' in data[i][times[-1]]:
        end_with_break = end_with_break + 1
        complete = False
    elif 'merge_end' in data[i][times[-1]]:
        end_with_merge = end_with_merge + 1
        complete = False
    else:
        clean_end = clean_end + 1
        
    if ('break_start' in data[i][times[0]]) and ('break_end' in data[i][times[-1]]):
        start_and_end_broken = start_and_end_broken + 1
        
    if complete: complete_clouds.append(i)
    

print "Split start: ", start_with_split
print "Break start: ", start_with_break
print "New start: ", start_new

print "Break end: ", end_with_break
print "Merge end: ", end_with_merge
print "Clean end: ", clean_end

print "Start/End broken: ", start_and_end_broken
print "Start/End split/merge: ", start_and_end_split_merge
print "Condensed clouds: ", condensed_clouds
print "Core clouds: ", core_clouds

cPickle.dump(complete_clouds, open('pkl/complete_clouds.pkl','wb'))

