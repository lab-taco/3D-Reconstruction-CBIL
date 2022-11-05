from .Sample_Cell_Gen import *
import numpy as np
import sys

def basic_connectivity_statistics(Cells, Cell_type='GC'):
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
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
    mean=[]
    std=[]
    for data in data_to_plot:
        mean.append(np.round(np.mean(data[0]), 3))
        std.append(np.round(np.std (data[0]), 3))
        label=data[-1] #,+'u:'+str(mean)+' std:'+str(std)
        ax.plot(data[0], data[1], label=label)
    print('Broadness Difference:', std[1]-std[0])
    #plt.legend()
    lgd = ax.legend(bbox_to_anchor=(1, 0.5), loc="lower left")
    plt.title('Cumulative Distribution Comparison')
    #fig.savefig('SP65.png', dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.show()

def cumulative_distribution(distribution, label_cdf='CDF', print_dist=False, smoothing=True):
    if smoothing:
        distribution+=np.random.normal(0,1, size=len(distribution))
    from statsmodels.distributions.empirical_distribution import ECDF
    ecdf = ECDF(distribution)
           
    if print_dist:
        import matplotlib.pyplot as plt 
        plt.title('ECDF')
        plt.plot(ecdf.x, ecdf.y, label=label_cdf)
        plt.legend()
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
    return random_mean, random_std, Sample_ratio_dist


    #return MFs, GCs
def GC_synaptic_partner_exchange(MFs, GCs, GC_colormap, Size, D): 
    #target_subgroups: 0=early, 2 =late
    #Size = Size of exchanging cells for each group eg. Size 10: each 10 GCs from early and late
    #D= the number of exchange round; the degree of exchange
    for _ in range(D):
        GCs=np.array(GCs,dtype=object)
        Early= np.array([gc for gc in GCs if gc.color==GC_colormap[0]])
        Late = np.array([gc for gc in GCs if gc.color==GC_colormap[-1]])
        Ind_choice_from_Early = np.random.choice(np.arange(len(Early)), size=Size, replace=False)
        Ind_choice_from_Late = np.random.choice (np.arange(len(Late)), size=Size, replace=False)
        #print('Ind_choice_from_Early', Ind_choice_from_Early)
        #print('Ind_choice_from_Late', Ind_choice_from_Late)
        #print('E:', np.array([gc.synapse_partners for gc in Early[Ind_choice_from_Early]]) )
        #print('L:', np.array([gc.synapse_partners for gc in Late[Ind_choice_from_Late]]) )
        for i in range(Size):
            exchange = Early[Ind_choice_from_Early[i]].synapse_partners
            Early[Ind_choice_from_Early[i]].synapse_partners = \
                Late[Ind_choice_from_Late[i]].synapse_partners
            Late[Ind_choice_from_Late[i]].synapse_partners= exchange
        #print('exchanged')
        #print('E:', np.array([gc.synapse_partners for gc in Early[Ind_choice_from_Early]]) )
        #print('L:', np.array([gc.synapse_partners for gc in Late[Ind_choice_from_Late]]) )
        for mf in MFs: #Initialization
            mf.synapse_partners=[]

        for gc_ind, gc in enumerate(GCs):
            for syn_pt in gc.synapse_partners:
                MFs[syn_pt].synapse_partners.append(gc_ind)
        
def Edge_initialization(MFs, GCs, Edges):
    for ind, edges in enumerate(Edges):
        syn_partners =[]
        for ed in edges:
            syn_partners.append(ed[1])
        MFs[ind].synapse_partners=syn_partners

    for gc in GCs:
        gc.synapse_partners=[]

    for ind_mf, mf in enumerate(MFs):
        for syn_part in mf.synapse_partners:
            #print('syn_part', syn_part)
            GCs[syn_part].synapse_partners.append(ind_mf)


def Statistics_Shuffled_Cells(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap,\
                             exchange_size = 10, D_rounds=1, Num_Repeat=10):
    Averaging_mean_shuffled=[]
    Averaging_STD_shuffled=[]
    
    Initial_edges = extract_edges2(Replica_MF_Objects)
    initial_ratio= connectivity_ratio_distribution(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap)
    mean, std = Statistics_distribution(initial_ratio)
    #print('initial m:', mean, 'std', std)
    for _ in range(Num_Repeat):      
        GC_synaptic_partner_exchange(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap,\
                             exchange_size, D_rounds)
        #shufflededges = extract_edges2(Replica_MF_Objects)
        ratio_distribution_shuffled = connectivity_ratio_distribution(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap)
        mean, std = Statistics_distribution(ratio_distribution_shuffled)
        Averaging_mean_shuffled.append(mean)
        Averaging_STD_shuffled.append(std)

        Edge_initialization(Replica_MF_Objects, Replica_GC_Objects, Initial_edges)
        initial_ratio= connectivity_ratio_distribution(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap)
        mean, std = Statistics_distribution(initial_ratio)
        #print('initial m:', mean, 'std', std)

    #for i in range(Num_Repeat):
    #    print('Shuffled m:', Averaging_mean_shuffled[i], 'std', Averaging_STD_shuffled[i])
    mean_shuffled, std_shuffled= np.mean(Averaging_mean_shuffled), np.mean(Averaging_STD_shuffled)
    return ratio_distribution_shuffled, mean_shuffled, std_shuffled

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
