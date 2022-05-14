from .Sample_Cell_Gen import *
import numpy as np
import sys

def connectivity_statistics(Cells, Cell_type='GC'):
    num_syn_partners= [len(c.synapse_partners) for c in Cells]
    print('per', Cell_type,\
         '  mean:', np.mean(num_syn_partners), 'std:', np.round(np.std(num_syn_partners),2))
    

"""
def extract_edges(MFs): # edges element structure: [Ind MF, Ind GC], Type: <Class, 'List'>
    edges_for_all_MFs=[]
    for ind_mf, mf in enumerate(MFs):
        indices_synapsed=mf.synapse_partners
        edges_mf=[] #element structure: [index_mf, index_partner_gcs]
        for ind_gc in indices_synapsed:
            edges_mf.append([ind_mf, ind_gc])
            #edges_mf.append([index, GCs.index(syn_partner), [GCs[GCs.index(syn_partner)].color] ]) #color?       
        #edges_mf.append(mf.color)
        edges_for_all_MFs.append(edges_mf)
    return edges_for_all_MFs"""

def extract_edges2(Cells): # generalized edges element structure: [Ind Cell, Ind Syn Partner], Type: <Class, 'List'>
    edges_for_all_Cells=[]
    for ind_cell, cell in enumerate(Cells):
        indices_synapsed=cell.synapse_partners
        edges_cell=[] 
        for ind_syn_partner in indices_synapsed:
            edges_cell.append([ind_cell, ind_syn_partner])
        edges_for_all_Cells.append(edges_cell)
    return edges_for_all_Cells

def print_dist_stats(distribution):
    print('{:<5}'.format('min:'), '{:>10}'.format(np.round(min(distribution),3)),\
          '{:<5}'.format('max:'), '{:>10}'.format(np.round(max(distribution),3)),\
          '{:<6}'.format('mean:'), '{:>10}'.format(np.round(np.mean(distribution),4)),\
          '{:<5}'.format('var:'), '{:>15}'.format(np.round(np.var(distribution),4)),\
          '{:<5}'.format('std:'), '{:>10}'.format(np.round(np.std(distribution),4))) 

def plot_distributions_together(data_to_plot):
    import matplotlib.pyplot as plt
    for data in data_to_plot:
        mean=np.round(np.mean(data[0]), 3)
        var =np.round(np.std (data[0]), 3)
        label=data[-1]+'u:'+str(mean)+' var:'+str(var)
        plt.plot(data[0], data[1], label=label)
    plt.title('Cumulative Distribution Comparison')
    plt.legend()
    plt.show()

def cumulative_distribution(distribution, print_dist=False):
    from statsmodels.distributions.empirical_distribution import ECDF
    ecdf = ECDF(distribution)    
    if print_dist:
        import matplotlib.pyplot as plt 
        plt.title('ECDF')
        plt.plot(ecdf.x, ecdf.y)
        plt.show()
    ##To avoid a bug in ECDF to produce nan or inf
    infinite = np.where(np.isfinite(ecdf.x)==0)
    x= np.delete(ecdf.x, infinite)
    y= np.delete(ecdf.y, infinite)
    return x, y
    


def connectivity_ratio_distribution(MFs, GCs, gc_color_map):
    ratios=[] 
    for ind_mf, mf in enumerate(MFs):
        score_early=0
        score_late=0
        indices_synapsed=mf.synapse_partners
        if len(indices_synapsed)==0: print('No synapsed partners for mf', ind_mf)
        for ind_gc in indices_synapsed:
            if GCs[ind_gc].color==gc_color_map[0]: score_early+=1
            elif GCs[ind_gc].color==gc_color_map[-1]: score_late+=1    
        ratio_mf= score_early-score_late        
        ratios.append(ratio_mf)
    return ratios

def Statistics_distribution(dist):
    mean=np.mean(dist)
    std=np.std(dist)
    #var=np.var(dist)
    return mean, std

def Statistics_Randomly_Connected_Cells(Num_MFs, Num_GCs, MF_Colormap, GC_Colormap, Num_Repeat=10):
    Averaging_mean_randomnet=[]
    Averaging_STD_randomnet=[]
    
    for _ in range(Num_Repeat):        
        Sample_MFs, Sample_GCs = randomly_conencted_sample_cells(Num_MFs, Num_GCs, MF_Colormap, GC_Colormap)
        Sample_ratio_dist=connectivity_ratio_distribution(Sample_MFs, Sample_GCs, GC_Colormap)
        mean, std= Statistics_distribution(Sample_ratio_dist)
        Averaging_mean_randomnet.append(mean)
        Averaging_STD_randomnet.append(std)
    random_mean, random_std= np.mean(Averaging_mean_randomnet), np.mean(Averaging_STD_randomnet)
    return random_mean, random_std

def shuffling(node_GC, node_MF, edges, gc_color_map, eval_shf_func=False): #shuffle the synapses between target GC groups and MF groups    
    
    total_early_gc_edges2=[]
    for mf in node_MF: # For selected MFs (sub) group
        early_gc_edges=[]
        for edge in [eg for eg in edges if eg[0]==mf]: # For whole edges with the mf            
            if gc_color_map[node_GC.index(edge[1])]==GC_colors[0]: #if synapsed with early GCs                
                early_gc_edges.append(edge)                               
        total_early_gc_edges2.append(early_gc_edges) #early GC edges list per MF
    edges=np.array(edges)        
    for edge_list in total_early_gc_edges2:
        if len(edge_list)>0:
            for edge in edge_list:
                edges=np.delete(edges, np.where(np.all(edges==edge, axis=1)), axis=0)
        
    np.random.shuffle(total_early_gc_edges2)     
    
    for ind, edge_list in enumerate(total_early_gc_edges2):
        if len(edge_list)>0:
            for edge in edge_list:                
                edge[0]='M%s'%ind
                #print(np.shape(edges))
                #print(np.shape([edge]))
                #edges=np.append(edges, [edge], axis=0)
       
    
    rebuilt_edges=[]    
    for ind, edge_list in enumerate(total_early_gc_edges2):
        if len(edge_list)>0:
            for edge in edge_list:                
                rebuilt_edges.append(edge)
    
    edges=np.concatenate((rebuilt_edges, edges), axis=0)    
    edges=np.ndarray.tolist(edges)
    return edges

def synaptic_partner_shuffling(MFs, GCs, GC_colormap, target_subgroup): #target_subgroups: 0=early, 2 =late
    GC_target_group = [gc for gc in GCs if gc.color==GC_colormap[target_subgroup]]
    target_edges = extract_edges2(GC_target_group) # [GC, MF] pairing
    #np.random.shuffle(target_edges)
    permuted_edges=np.random.permutation(target_edges)

    for ind, edge in enumerate(permuted_edges):
        ind_partner_MF=edge[1]
        GC_target_group[ind].synapse_partners=ind_partner_MF
    
    for mf in MFs: #Initialization
        mf.synapse_partners=[]

    for gc_ind, gc in enumerate(GCs):
        for syn_pt in gc.synapse_partners:
            MFs[syn_pt].synapse_partners.append(gc_ind)

    #return MFs, GCs
        

def Statistics_Shuffled_Cells(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap, Num_Repeat=10):
    Averaging_mean_shuffled=[]
    Averaging_STD_shuffled=[]
    
    edges1 = extract_edges2(Replica_MF_Objects)
    for _ in range(Num_Repeat):
        synaptic_partner_shuffling(Replica_MF_Objects, Replica_GC_Objects,  GC_Colormap, target_subgroup=0)
        shufflededges = extract_edges2(Replica_MF_Objects)
        ratio_distribution_shuffled = connectivity_ratio_distribution(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap)
        mean, std = Statistics_distribution(ratio_distribution_shuffled)
        print('m:', mean, 's:',std)
        Averaging_mean_shuffled.append(mean)
        Averaging_STD_shuffled.append(std)
    mean_shuffled, std_shuffled= np.mean(Averaging_mean_shuffled), np.mean(Averaging_STD_shuffled)   
    return mean_shuffled, std_shuffled

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