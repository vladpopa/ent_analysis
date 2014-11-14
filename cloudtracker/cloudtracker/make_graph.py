#!/usr/bin/env python
import cPickle
import networkx

full_output=False

def full_output(cloud_times, cloud_graphs, merges, splits, MC):
    cPickle.dump(cloud_times, open('pkl/cloud_times.pkl','wb'))

    n = 0
    clouds = {}
    for subgraph in cloud_graphs:
        # Clouds include plume only clusters
        events = {'has_condensed': False, 'has_core': False}
        for node in subgraph:
            node_events = []
            t = int(node[:8])
            if t == 0: node_events.append('break_start')
            if t == (MC['nt']-1): node_events.append('break_end')

            if node in merges:
                node_events.append('merge_end')
            if node in splits:
                node_events.append('split')
                
            info = subgraph.node[node]
            if info['merge']:
                node_events.append('merge')
            if info['split']:
                node_events.append('split_start')
                
            events['has_condensed'] = (info['condensed'] > 0) or events['has_condensed']
            events['has_core'] = (info['core'] > 0) or events['has_core']

            if t in events:
                events[t].extend(node_events)
            else:
                events[t] = node_events[:]

        # Store cloud events and increment cloud id
        clouds[n] = events
        n = n + 1

    cPickle.dump(clouds, open('pkl/graph_events.pkl', 'wb'))

def make_graph(MC):
    graph = networkx.Graph()

    merges = {}
    splits = {}

    for t in range(MC['nt']):
        clusters = cPickle.load(open('pkl/clusters_%08g.pkl' % t, 'rb'))

        for id in clusters:
            # Make dictionaries of every split and every merge event that occurs
            # in a cluster's lifecycle
            m_conns = clusters[id]['merge_connections']
            s_conns = clusters[id]['split_connections']
            core = len(clusters[id]['core'])
            condensed = len(clusters[id]['condensed'])
            plume = len(clusters[id]['plume'])
            attr_dict = {'merge': m_conns,
                         'split': s_conns,
                         'core': core,
                         'condensed': condensed,
                         'plume': plume}
                         
            for item in m_conns:
                node1 = '%08g|%08g' % (t, id)
                node2 = '%08g|%08g' % (t-1, item)
                merges[node2] = node1
            for item in s_conns:
                node1 = '%08g|%08g' % (t, id)
                node2 = '%08g|%08g' % (t, item)
                splits[node2] = node1            
        
            # Construct a graph of the cluster connections
            graph.add_node('%08g|%08g' % (t, id), attr_dict = attr_dict)
            if clusters[id]['past_connections']:
                for item in clusters[id]['past_connections']:
                    graph.add_edge('%08g|%08g' % (t-1, item),
                                   '%08g|%08g' % (t, id))

    # Iterate over every cloud in the graph
    for subgraph in networkx.connected_component_subgraphs(graph):
        # Find the duration over which the cloud_graph has cloudy points.
        times = set()
        for node in subgraph:
            if subgraph.node[node]['condensed'] > 0:
                times.add(int(node[:8]))

        # If cloud exists for less than 5 minutes, check if it has split events
        # If it has split events, remove them and reconnect the cloud
        if (len(times) < 5*60/MC['dt']) & (len(times) > 0):
            for node in subgraph:
                if subgraph.node[node]['split']:
                    item = subgraph.node[node]['split'].pop()
                    t = int(node[:8])
                    graph.add_edge(node, '%08g|%08g' % (t, item))
 
    for subgraph in networkx.connected_component_subgraphs(graph):
        # Find the duration over which the cloud_graph has cloudy points.
        times = set()
        for node in subgraph:
            if subgraph.node[node]['condensed'] > 0:
                times.add(int(node[:8]))

        # If a cloud exists less than 5 minutes, check for merge events
        if (len(times) < 5*60/MC['dt']) & (len(times) > 0):
            for node in subgraph:
                if subgraph.node[node]['merge']:
                    item = subgraph.node[node]['merge'].pop()
                    t = int(node[:8])
                    graph.add_edge(node, '%08g|%08g' % (t-1, item))

    cloud_times = []    
    cloud_graphs = []
    cloud_noise = []
    for subgraph in networkx.connected_component_subgraphs(graph):
        # If a graph exists less than 2 minutes, classify it as noise
        # Otherwise, put it in cloud_graphs
        plume_time = set()
        condensed_time = set()
        core_time = set()
        # Cloud total volumes over lifetime
        plume_volume = 0
        condensed_volume = 0
        core_volume = 0
        for node in subgraph.nodes():
            condensed_vol = subgraph.node[node]['condensed'] 
            if condensed_vol > 0:
                condensed_volume = condensed_volume + condensed_vol
                time = int(node[:8])
                condensed_time.add(time)

            core_vol = subgraph.node[node]['core'] 
            if core_vol > 0:
                core_volume = core_volume + core_vol
                time = int(node[:8])
                core_time.add(time)

            plume_vol = subgraph.node[node]['plume'] 
            if plume_vol > 0:
                plume_volume = plume_volume + plume_vol
                time = int(node[:8])
                plume_time.add(time)

        # If a graph exists less than 2 minutes, classify it as noise.     
        if len(plume_time) < 2*60/MC['dt']:
            cloud_noise.append(subgraph)
        else:
            plume_time = list(plume_time)
            plume_time.sort()
            condensed_time = list(condensed_time)
            condensed_time.sort()
            core_time = list(core_time)
            core_time.sort()
            cloud_graphs.append((plume_volume, subgraph, plume_time, condensed_time, core_time))
    
    # Cloud with largest condensed volume sorted largest to smallest; volume not 
    # used any further
    cloud_graphs.sort()
    cloud_graphs.reverse()
    cloud_times = [item[2:] for item in cloud_graphs]  
    cloud_graphs = [item[1] for item in cloud_graphs]
    if full_output: full_output(cloud_times, cloud_graphs, merges, splits, MC)

    return cloud_graphs, cloud_noise