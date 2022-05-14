#from Simulation.form_synapse import *
import imp
from Simulation.Data_Extraction import *

from Analysis.Analysis_Connectivity import *

from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt

DEFEAULT_GC_object_data= 'save_cell_objects_Num_parents10_time13-05-2022-1847_GCs'
DEFEAULT_MF_object_data= 'save_cell_objects_Num_parents10_time13-05-2022-1847_MFs'


def main(GC_data_name, MF_data_name, ANALYSIS_SPATIAL_DISTRIBUTION, ANALYSIS_CONNECTIVITY):

    print('-------Data Loading...--------------------------------------------------')
    GC_objects = np.ndarray.tolist(data_load(GC_data_name, DATA_PATH))
    MF_objects = np.ndarray.tolist(data_load(MF_data_name, DATA_PATH))
    print('GC_objects type:', type(GC_objects), 'data shape:', np.shape(GC_objects))
    print('MF_objects type:', type(MF_objects), 'data shape:', np.shape(MF_objects))


    if ANALYSIS_SPATIAL_DISTRIBUTION:
        print('------ANALYSIS_SPATIAL_DISTRIBUTION...-------------------------------')
        

    if ANALYSIS_CONNECTIVITY:
        print('------ANALYSIS_CONNECTIVITY...---------------------------------------')
        print('Statistics of the number of synaptic partners...')
        connectivity_statistics(GC_objects, Cell_type='GC')
        connectivity_statistics(MF_objects, Cell_type='MF')
        print('')

        print('Edge analysis...')
        Edges = extract_edges(MF_objects)       
        print('Edge type:', type(Edges), 'Edges shape:', np.shape(np.array(Edges, dtype=object))\
            , 'Edge[0] type:', type(Edges[0])) 
        #To use Numpy with the Lists with varying size of rows, use np.array(list, dtype=object)
        print('')




import argparse
import time
if __name__ == '__main__':

    start_time = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--GC_object_name', type=str, dest='data_name_GCs',
                        default=DEFEAULT_GC_object_data,
                        help='The saved GC Object data name')
    
    parser.add_argument('--MF_object_name', type=str, dest='data_name_MFs',
                        default=DEFEAULT_MF_object_data,
                        help='The saved MF Object data name')
    
    #boolen for argparse is tricky, So used 0,1 instead

    parser.add_argument('--Spatial_analysis', type=int, dest='ANALYSIS_SPATIAL_DISTRIBUTION',
                        default=1,
                        help='Running Spatial_analysis of cell soma distribution')
    parser.add_argument('--Connectivity_analysis', type=int, dest='ANALYSIS_CONNECTIVITY',
                        default=1,
                        help='Running Connectivity analysis of GC-MF synapses')    

    args = parser.parse_args()
    
    #main(args.Num_childs, args.draw_dist, args.draw_net)
    main(args.data_name_GCs, args.data_name_MFs, \
        args.ANALYSIS_SPATIAL_DISTRIBUTION, args.ANALYSIS_CONNECTIVITY)

    elapsed_time = time.time() - start_time
    print('Total elapsed time for analysis:', elapsed_time)  