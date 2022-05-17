import numpy as np
import sys
from .Sample_Cell_Gen import *
from .Analysis_Connectivity import *
import matplotlib.pyplot as plt

def synaptic_partner_exchange(MF1, MF2, GC1, GC2, Size, D): 
    #target_subgroups: 0=early, 2 =lates
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

def overlapped_network(Num_MFs, Num_GCs, colormap_mf, colormap_gc, Size, D):
    
    MF1, GC1=randomly_conencted_sample_cells(int(Num_MFs/2), int(Num_GCs/2), colormap_mf, [colormap_gc[0]])
    MF2, GC2=randomly_conencted_sample_cells(int(Num_MFs/2), int(Num_GCs/2), colormap_mf, [colormap_gc[-1]]\
                                            , add_MF=int(Num_MFs/2), add_GC= int(Num_GCs/2))

    MFs = MF1+MF2
    GCs = GC1+GC2
    
    #Initial_edges = extract_edges2(MFs)
    initial_ratio= connectivity_ratio_distribution(MFs, GCs, colormap_gc)
    mean, std = Statistics_distribution(initial_ratio)
    #synaptic_partner_exchange(MF1, MF2, GC1, GC2, Size, D)
    #print('initial edges:', np.array(Initial_edges))
    print('initial m:', mean, 'std', std)
    for gc in GCs:
        if len(gc.synapse_partners)!=4: print(len(gc.synapse_partners))
    
    trend=[]
    stds=[]
    for i in range(1000):      
        GC_synaptic_partner_exchange(MFs, GCs, colormap_gc, Size, D)
        #shufflededges = extract_edges2(Replica_MF_Objects)
        ratio_distribution_shuffled = connectivity_ratio_distribution(MFs, GCs, colormap_gc)
        mean, std = Statistics_distribution(ratio_distribution_shuffled)
        #Swapped_edges = extract_edges2(MFs)
        print('Swapped m:', mean, 'std', std)
        stds.append(std)
        for gc in GCs:
            if len(gc.synapse_partners)!=4: print(len(gc.synapse_partners))
        trend.append([i, std])
        #print('Swapped edges:', np.array(Swapped_edges))
    #cumulative_distribution(np.array(trend).T, label_cdf='Overlapped', print_dist=True)
    #plt.plot(np.array(trend)[:,0], np.array(trend)[:1])
    plt.plot(np.array(trend)[:,0], np.array(trend)[:,1])
    print('np.mean(stds)', np.mean(stds))
    plt.axhline(y=np.mean(stds), color='r', linestyle='-')
    plt.show()