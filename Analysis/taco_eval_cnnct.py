import time
from degree_dist import *
import os
#DEFEAULT_edge_data='Volume_filling_155'
DEFEAULT_edge_data='tolerance_170_MFs_edges'
#DEFEAULT_edge_data='empty'
#DEFEAULT_edge_data='tolerance_10_MFs_edges13-02-2022-1843'
DEFEAULT_DRAW_DIST=True
DEFEAULT_DRAW_NET=False

def main(edge_data, draw_dist, draw_net):
    if edge_data=='empty':
        print('No edge data')
    else:
        edge_data=data_load(edge_data)
        print('edge data shape:', np.shape(edge_data))
        synapse_data=edge_data[:-1]
        num_GC=edge_data[-1][1]

        """
        Record edges information ----------------------------------------------------
        """ 
        ed_list=[]        
        for mf in synapse_data:
            for synapse in mf[:-1]:
                edge=['M%s'%synapse[0], 'G%s'%synapse[1]]
                ed_list.append(edge)

        """
        Record nodes information ----------------------------------------------------
        """ 
        nd_list1=[]
        nd_list2=[]
        for ind_gc in range(num_GC):
            nd_list1.append('G%s'%ind_gc)
        for ind_mf in range(len(synapse_data)):
            nd_list2.append('M%s'%ind_mf)

        """
        # Cells ----------------------------------------------------
        """ 
        print('# GCs:', len(nd_list1), '# MFs:', len(nd_list2))
        print('Early GC: %d,  Mid GC: %d,  Late GC: %d '\
                %(int(num_GC/3), int(num_GC/3), num_GC-2*(int(num_GC/3))))

        
        """
        Stats for whole MF edge ratio with shuffling for randomness benchmark ---------------------------------
        """ 
        '''
        mean, std = neuralnet_whole(nd_list1, nd_list2, ed_list,frequency_dist=draw_dist, \
            draw_net=draw_net, shuffle=True)
            #draw_net=draw_net)'''


        """
        Evalution of connectivity from edge data ----------------------------------------------------
        """ 
        neuralnet4(nd_list1, nd_list2, ed_list,  \
            MF_subgroup='early', frequency_dist=draw_dist, draw_net=draw_net)
        neuralnet4(nd_list1, nd_list2, ed_list,  \
            MF_subgroup='mid', frequency_dist=draw_dist, draw_net=draw_net)
        neuralnet4(nd_list1, nd_list2, ed_list,  \
            MF_subgroup='late', frequency_dist=draw_dist, draw_net=draw_net)
    

import argparse
if __name__ == '__main__':

    start_time = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--edge_data', type=str, dest='edge_data',
                        default=DEFEAULT_edge_data,
                        help='Data of edges that represent GC-MF network connectivity.')
    
    parser.add_argument('--draw_dist', type=int, dest='draw_dist',
                        default=DEFEAULT_DRAW_DIST,
                        help='Draw distribution of network connectivity\
                              input 1 for on')

    parser.add_argument('--draw_net', type=int, dest='draw_net',
                        default=DEFEAULT_DRAW_NET,
                        help='Draw network of given edges by using networkX library.\
                              input 1 for on')

    args = parser.parse_args()
    
    #main(args.Num_childs, args.draw_dist, args.draw_net)
    main(args.edge_data, args.draw_dist, args.draw_net)

    elapsed_time = time.time() - start_time
    print('Total elapsed time for evaluation:', elapsed_time)