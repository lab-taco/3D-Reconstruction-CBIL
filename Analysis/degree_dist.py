from matplotlib import pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import numpy as np
import scipy.stats as stats
import math
import os
import sys
from vpython import *
from function2 import * #load data
MF_colors=np.array([[192,192,192], [135,121,78], [109,135,100]])/256
GC_colors = ['red', 'orange', 'yellow']    

#num_synapses_range
def input_processing(data_name):
    load = data_load(data_name)
    synapse_data=load[:-1] #load[-1] : GC number
    num_GC=load[-1][1]

    ed_list=[]
    #mf_color_map=[]
    for mf in synapse_data:
        #mf_color_map.append([mf[-1].x, mf[-1].y, mf[-1].z])
        for synapse in mf[:-1]:
            edge=['M%s'%synapse[0], 'G%s'%synapse[1]]
            ed_list.append(edge)

    node_GC=[]
    node_MF=[]
    for ind_gc in range(num_GC):
        node_GC.append('G%s'%ind_gc)
    for ind_mf in range(len(synapse_data)):
        node_MF.append('M%s'%ind_mf)
    print('# GCs', len(node_GC))
    print('# MFs', len(node_MF))
    
    return node_GC, node_MF, ed_list

def normal_dist():
    mu = 0
    variance = 1
    sigma = math.sqrt(variance)
    x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
    plt.plot(x, stats.norm.pdf(x, mu, sigma))

def bins_labels(bins, **kwargs):
    bin_w = (max(bins) - min(bins)) / (len(bins) - 1)
    plt.xticks(np.arange(min(bins)+bin_w/2, max(bins), bin_w), bins, **kwargs)
    plt.xlim(bins[0], bins[-1])

def broadness_difference(node_GC, node_MF, ed_list):
    unshuffled=neuralnet4(node_GC, node_MF, ed_list, \
            MF_subgroup='total', shuffle=False, rt_ratio=True, print_stat=False)
    std_unshfld=np.std(unshuffled)

    print('std_unshfld', std_unshfld)

    std_diff_list=[]
    for i in range(20):
        #print(i, '-th shuffling...')
        shuffled=neuralnet4(node_GC, node_MF, ed_list, \
            MF_subgroup='total', shuffle=True, rt_ratio=True, print_stat=False)
        std_diff = np.std(shuffled) - std_unshfld
        std_diff_list.append(std_diff)
        print(i, '-th shuffling...: std_shuffled', np.std(shuffled), 'std_diff', std_diff)
    
    return np.mean(std_diff_list)
    
def range_for_plot(distribution, range):
    range_max=int(max(abs(distribution)))+1
    dist_range= np.arange(-range_max, range_max)
    return dist_range

def CDF_comparison_of_shuffling(node_GC, node_MF, ed_list):
    unshuffled=neuralnet4(node_GC, node_MF, ed_list, \
            MF_subgroup='total', shuffle=False, rt_ratio=True, print_stat=False)
    shuffled=neuralnet4(node_GC, node_MF, ed_list, \
                MF_subgroup='total', shuffle=True, rt_ratio=True, print_stat=False)

    x1, val1 = plot_degree_dist(unshuffled, plot=False, CDF_return2=True)
    x2, val2 = plot_degree_dist(shuffled  , plot=False, CDF_return2=True)
    
    plt.title('plot together')
    plt.plot(x1, val1)
    plt.plot(x2, val2)
    plt.show()


def plot_degree_dist(distribution, normal_dist=False, plot=True, \
                    moments_return=False, CDF_return2=False):
    #x=[1,2,3,4,5]
    #y=[2,4,6,8,10]
    #print('variance:', np.var(distribution))
    
    if normal_dist==True:
        mu = 0
        #variance = 8.36514
        variance = 1
        sigma = math.sqrt(variance)
        x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
        plt.plot(x, stats.norm.pdf(x, mu, sigma))

    print_dist_stats(distribution)
    
    from statsmodels.distributions.empirical_distribution import ECDF
    ecdf = ECDF(distribution)
    if plot==True:        
        plt.title('ECDF')
        plt.plot(ecdf.x, ecdf.y)
        plt.show()
    #sys.exit()

    range_max=int(max(abs(distribution)))+1
    log_range_= np.arange(-range_max, range_max)
    #log_range_= np.arange(int(min(distribution)-2), int(max(distribution)+3))
    if CDF_return2: return ecdf.x, ecdf.y
    if plot==True:
        plt.title('Log scale Degree of frequency distibution')
        #print('log min:', min(log_counts), 'log max:', max(log_counts), \
        #        'mean:', np.mean(log_counts), 'var:', np.var(log_counts), 'std:', np.std(log_counts))
        #print('log min:', np.round(min(log_counts),3), 'log max', np.round(max(log_counts),3), \
        
        dist_of_probability=True        
        plt.xlim(-range_max, range_max)
        plt.hist(distribution, bins=log_range_, density=dist_of_probability, align='mid', cumulative=False)
        plt.show()
        

    if moments_return==True:
        std=np.round(np.std(distribution),4)
        mean=np.round(np.mean(distribution),4)
        return mean, std
    

''' 
#sample drawing
def neuralnet():
    from networkx.drawing.nx_agraph import graphviz_layout

    G = nx.DiGraph()
    ed = [[0, 4, -1],
     [0, 5, -1],
     [1, 4, -1],
     [1, 5, -1],
     [2, 4, -1],
     [2, 5, 10],
     [4, 3, -1],
     [5, 3, 100]]

    G.add_weighted_edges_from(ed)
    pos = graphviz_layout(G, prog='dot', args="-Grankdir=LR")
    nx.draw(G,with_labels=True,pos=pos, font_weight='bold')
    plt.show()
'''

#refer networkx_Test.py instead of this

def migration_timing(cells, colors):
    color_map = []    
    stat=False
    if stat==True:
        e, m, l =0, 0, 0
    for index, cell in enumerate(cells):
        if index<int(len(cells)/3): #early MF
            color=colors[0]
            if stat==True: e+=1
            #print('Early', index)
        elif index<int(len(cells)*2/3): #mid MF
            color=colors[1]
            if stat==True: m+=1
            #print('Mid', index)
        else:                               #late MF
            color=colors[2]
            if stat==True: l+=1
            #print('Late', index)
        color_map.append(color)        
    if stat==True: print('Colormap, Early: %d,  Mid: %d,  Late: %d '%(e, m, l))

    
    return color_map

def node_pruning(indices, color_map, toatal_cells):
    #remove_list=[]
    remove_indices=[]
    for index, cell in enumerate(toatal_cells):
        if not index in indices:
            #remove_list.append(cell)
            remove_indices.append(index)
    
    target_cells=np.ndarray.tolist(np.delete(toatal_cells, remove_indices))
    color_map=np.ndarray.tolist(np.delete(color_map, remove_indices, axis=0))

    return target_cells, color_map

def edge_pruning(remaining_nodes, edges):
    remove_list_edge=[]
    for edge in edges:
        if not edge[0] in remaining_nodes and edge not in remove_list_edge:
            remove_list_edge.append(edge)
    eds=[edg for edg in edges if not edg in remove_list_edge]
    #print('remove_list_edge', len(remove_list_edge))
    return eds

def print_dist_stats(distribution):
    print('{:<5}'.format('min:'), '{:>10}'.format(np.round(min(distribution),3)),\
          '{:<5}'.format('max:'), '{:>10}'.format(np.round(max(distribution),3)),\
          '{:<6}'.format('mean:'), '{:>10}'.format(np.round(np.mean(distribution),4)),\
          '{:<5}'.format('var:'), '{:>15}'.format(np.round(np.var(distribution),4)),\
          '{:<5}'.format('std:'), '{:>10}'.format(np.round(np.std(distribution),4))) 


def frequency_distribution(node_MF, node_GC, edges, gc_color_map):
    total_syn=[]
    for mf in node_MF:
        syn_for_mf=[]
        for edge in edges:
            if edge[0]==mf:
                syn_for_mf.append(edge)
        total_syn.append(syn_for_mf)

    ratios=[]
    for syns_per_mf in total_syn:
        early=0
        late=0
        for syn in syns_per_mf:
            ind_GC=node_GC.index(syn[1])
            if gc_color_map[ind_GC]==GC_colors[0]: early+=1
            elif gc_color_map[ind_GC]==GC_colors[2]: late+=1
        #ratio= (early+0.1)/(late+0.1)
        #noise1 = abs(np.random.normal(0, 1e-01))
        #noise2 = abs(np.random.normal(0, 1e-01))
        noise1 = abs(np.random.normal(0, 1))
        noise2 = abs(np.random.normal(0, 1))
        ratio= (early+noise1)/(late+noise2)
        #ratio = np.log(ratio)
        ratios.append(ratio)
    log_ratios=np.log(ratios)
    #log_ratios=ratios
    return log_ratios

def frequency_distribution_of_difference(node_MF, node_GC, edges, gc_color_map):
    total_syn=[]
    for mf in node_MF:
        syn_for_mf=[]
        for edge in edges:
            if edge[0]==mf:
                syn_for_mf.append(edge)
        total_syn.append(syn_for_mf)

    differences=[]
    for syns_per_mf in total_syn:
        early=0
        late=0
        for syn in syns_per_mf:
            ind_GC=node_GC.index(syn[1])
            if gc_color_map[ind_GC]==GC_colors[0]: early+=1
            elif gc_color_map[ind_GC]==GC_colors[2]: late+=1
        #ratio= (early+0.1)/(late+0.1)
        difference= early-late
        differences.append(difference)
    print('differences shape', np.shape(differences))
    return differences

def frecuency_MF(node_MF, node_GC, edges, gc_color_map):
    total_syn=[]
    for mf in node_MF:
        syn_for_mf=[]
        for edge in edges:
            if edge[0]==mf:
                syn_for_mf.append(edge)
        total_syn.append(syn_for_mf)

    counts=[]
    total_early=0
    total_mid=0
    total_late=0
    for syns_per_mf in total_syn:
        #early=0
        #mid=0
        #late=0
        for syn in syns_per_mf:
            ind_GC=node_GC.index(syn[1])
            if gc_color_map[ind_GC]==GC_colors[0]: total_early+=1
            elif gc_color_map[ind_GC]==GC_colors[1]: total_mid+=1
            elif gc_color_map[ind_GC]==GC_colors[2]: total_late+=1
        #print(early, late)
        #total_early+=early
        #total_mid+=mid
        #total_late+=late
        #count = [early, late]
        #counts.append(count)
    #avg_early=np.average(counts[0])
    #avg_late=np.average(counts[1])    
    #return total_early, total_late, avg_early, avg_late
    return total_early, total_mid, total_late

def frecuencyGC(node_MF, node_GC, edges, gc_color_map):
    Counts=[]
    for gc in node_GC:
        early=0
        late=0
        for edge in edges:
            if edge[1]==gc:
                ind=node_GC.index(gc)
                if gc_color_map(ind)==GC_colors[0]: early+=1
                elif gc_color_map[ind]==GC_colors[2]: late+=1  
        Counts.append([early, late])

def target_MF_indexing(node_MF, MF_subgroup):
    if MF_subgroup==None:
        print('Non MF subgroup selected, selecting 5 samples randomly')
        indices = np.random.choice(range(len(node_MF)), 5, replace=False)
    elif MF_subgroup=='early':
        indices = [index for index in range(len(node_MF))if index < int(len(node_MF)/3)]
    elif MF_subgroup=='mid':
        indices = [index for index in range(len(node_MF)) if index < int(len(node_MF)*2/3)\
            and index >=int(len(node_MF)/3)]
    elif MF_subgroup=='late':
        indices = [index for index in range(len(node_MF)) if index >= int(len(node_MF)*2/3)]
    elif MF_subgroup=='total':
        indices = [index for index in range(len(node_MF))]
    #print('# indices MFs',MF_subgroup,':', len(indices))
    return indices

def shuffling(node_GC, node_MF, edges, gc_color_map, eval_shf_func=False): #shuffle the synapses between target GC groups and MF groups    
    if eval_shf_func:
        print('edge before', [eg for eg in edges if gc_color_map[node_GC.index(eg[1])]==GC_colors[0]])
    total_early_gc_edges2=[]
    for mf in node_MF: # For selected MFs (sub) group
        early_gc_edges=[]
        for edge in [eg for eg in edges if eg[0]==mf]: # For whole edges with the mf            
            if gc_color_map[node_GC.index(edge[1])]==GC_colors[0]: #if synapsed with early GCs                
                early_gc_edges.append(edge)                               
        total_early_gc_edges2.append(early_gc_edges) #early GC edges list per MF
    if eval_shf_func:
        print('Original edges', np.shape(edges))
        for i in total_early_gc_edges2:
            print(i)
    edges=np.array(edges)        
    for edge_list in total_early_gc_edges2:
        if len(edge_list)>0:
            for edge in edge_list:
                edges=np.delete(edges, np.where(np.all(edges==edge, axis=1)), axis=0)
    if eval_shf_func: print('Deleted edges', np.shape(edges))
        
    np.random.shuffle(total_early_gc_edges2)

    if eval_shf_func:
        print('shuffled')
        for i in total_early_gc_edges2:
            print(i)     
    
    for ind, edge_list in enumerate(total_early_gc_edges2):
        if len(edge_list)>0:
            for edge in edge_list:                
                edge[0]='M%s'%ind
                #print(np.shape(edges))
                #print(np.shape([edge]))
                #edges=np.append(edges, [edge], axis=0)
    if eval_shf_func:
        print('sorted')
        for i in total_early_gc_edges2:
            print(i)     
    
    rebuilt_edges=[]    
    for ind, edge_list in enumerate(total_early_gc_edges2):
        if len(edge_list)>0:
            for edge in edge_list:                
                rebuilt_edges.append(edge)
    if eval_shf_func:
        print('reshaped')
        print(rebuilt_edges)
    edges=np.concatenate((rebuilt_edges, edges), axis=0)    
    edges=np.ndarray.tolist(edges)
    
    if eval_shf_func: 
        print('Shuffled edges', np.shape(edges))    
        print('early edge final', [eg for eg in edges if gc_color_map[node_GC.index(eg[1])]==GC_colors[0]])
    #sys.exit()
    return edges

def print_connectivity_stats(total_early, total_mid, total_late, MF_subgroup, node_MF):
    avg_early = np.around(total_early/len(node_MF),3)
    avg_mid = np.around(total_mid/len(node_MF),3)
    avg_late = np.around(total_late/len(node_MF), 3)    
    print('For ', MF_subgroup, 'MF subgroup,', '# node_MF', len(node_MF))
    print('# dend. from Early GC:', total_early, 'Average #Early GC per MF:', avg_early)
    print('# dend. from Mid GC:', total_mid, 'Average #Mid GC per MF:', avg_mid)
    print('# dend. from Late GC:', total_late, 'Average #Late GC per MF:', avg_late) 

def draw_network(gc_color_map, mf_color_map, edges, node_GC, node_MF):
    color_map=gc_color_map+mf_color_map    
    nodes= node_GC + node_MF

    G = nx.DiGraph() # Define the Graph and add the nodes/edges
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    if len(edges)>100:
        raise Exception("Too many edges to plot")
    nx.draw(G, pos=nx.bipartite_layout(G, node_MF),node_size=140,\
         scale= 2, node_color=color_map, with_labels=True)
    plt.show()

def neuralnet4(node_GC, node_MF, edges, MF_subgroup=None, \
            plot_frequency_dist=False, draw_net=False, shuffle=False, \
                rt_ratio=False, print_stat=True):    
    mf_color_map = migration_timing(node_MF, MF_colors)    
    gc_color_map = migration_timing(node_GC, GC_colors)
    
    indices = target_MF_indexing(node_MF, MF_subgroup)
    node_MF, mf_color_map = node_pruning(indices, mf_color_map, node_MF)
    edges                 = edge_pruning(node_MF, edges)
    if shuffle: edges = shuffling(node_GC, node_MF, edges, gc_color_map)
    
    if draw_net: draw_network(gc_color_map, mf_color_map, edges, node_GC, node_MF)
    
    ratios=frequency_distribution(node_MF, node_GC, edges, gc_color_map)    
    if plot_frequency_dist: plot_degree_dist(ratios)        

    total_early, total_mid, total_late = frecuency_MF(node_MF, node_GC, edges, gc_color_map)
    if print_stat: print_connectivity_stats(total_early, total_mid, total_late, MF_subgroup, node_MF)
    #print('variance:', np.var(ratios))    
    if rt_ratio: return ratios


def random_networks(Num_GCs, Num_MFs):    
    
    node_GC=[]
    node_MF=[]
    for gc in range(Num_GCs):        
        node_GC.append('G%s'%gc)
    for mf in range(Num_MFs):        
        node_MF.append('M%s'%mf)

    ed_list=[]
    for i in range(Num_GCs):
        synapse_partner=np.random.choice(Num_MFs, size=4, replace=False)
        syns=sorted(synapse_partner)
        for j in syns:            
            ed_list.append(['M%s'%j, 'G%s'%i])

    print('# GCs', len(node_GC))
    print('# MFs', len(node_MF))
    print('# edges(synapses)', len(ed_list))
    
    return node_GC, node_MF, ed_list




