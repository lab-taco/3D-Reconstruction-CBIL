import numpy as np
import sys

def connectivity_statistics(Cells, Cell_type='GC'):
    num_syn_partners= [len(c.synapse_partners) for c in Cells]
    print('per', Cell_type,\
         '  mean:', np.mean(num_syn_partners), 'std:', np.round(np.std(num_syn_partners),2))
    


def extract_edges(MFs): # edges element structure: [Ind MF, Ind GC], Type: <Class, 'List'>
    edges_for_all_MFs=[]
    for ind_mf, mf in enumerate(MFs):
        ind_synapsed=mf.synapse_partners
        edges_mf=[] #element structure: [index_mf, index_partner_gcs]
        for ind_gc in ind_synapsed:
            edges_mf.append([ind_mf, ind_gc])
            #edges_mf.append([index, GCs.index(syn_partner), [GCs[GCs.index(syn_partner)].color] ]) #color?       
        edges_mf.append(mf.color)        
        edges_for_all_MFs.append(edges_mf)
    return edges_for_all_MFs


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