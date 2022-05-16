#from Simulation.form_synapse import *
import imp
from Simulation.Data_Extraction import *

from Analysis.Analysis_Connectivity import *
from Analysis.Analysis_Spatial_dist import *
from Analysis.Slicing import *
from Simulation.Parameters import *
from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt

DEFEAULT_GC_object_data= 'save_cell_objects_Num_parents10_time13-05-2022-1847_GCs'
DEFEAULT_MF_object_data= 'save_cell_objects_Num_parents10_time13-05-2022-1847_MFs'

DATA_FOLDER='14-05-2022-1333'
PLOT_SPATIAL_DIST=False
PLOT_CUMU_DIST=False
def main(GC_data_name, MF_data_name, ANALYSIS_SPATIAL_DISTRIBUTION, ANALYSIS_CONNECTIVITY):

    print('-------Data Loading...----------------------------------------------------')
    GC_Objects, MF_Objects, GC_Colormap, MF_Colormap = load_all(DATA_FOLDER, DATA_PATH, Print=True)
    Num_GCs=len(GC_Objects)
    Num_MFs=len(MF_Objects)
    #len_synapse(GCs, MFs) -------------------------------------------
    
    if ANALYSIS_SPATIAL_DISTRIBUTION:
        Map_size_3D= [area_length, area_width, height_PCL]
        Map_size_2D= [area_length, height_PCL]
        Cell_radii=[radius_MFR, radius_GC]
        print('------ANALYSIS_SPATIAL_DISTRIBUTION...-------------------------------')
        #data_extraction_3d(GC_Objects, Map_size_2D, radius_GC, PLOTTING=True)
        #data_extraction_2d(GC_Objects, Map_size_2D, radius_GC, PLOTTING=True)
        
        MF_Captured, GC_Captured = Capture_Cells_at_Slice(MF_Objects, GC_Objects, \
                        Map_size_3D, Cell_radii, \
                        PLOTTING=PLOT_SPATIAL_DIST, Where_to_Slice='Random')
        Position_MFs=[[mf.body.pos.x, mf.body.pos.y] for mf in MF_Captured]
        Position_GCs=[[gc.body.pos.x, gc.body.pos.y] for gc in GC_Captured]
        #np.array(Position_MFs)[:, 0], np.array(Position_MFs)[:, 1]
        #np.array(Position_GCs)[:, 0], np.array(Position_GCs)[:, 1]
        random_dist = random_scattering(100, Map_size_2D, plotting=PLOT_SPATIAL_DIST)
        regular_dist = regular_scattering(num_grid=10, Map_size_2D=Map_size_2D, plotting=PLOT_SPATIAL_DIST)

    sys.exit()
    if ANALYSIS_CONNECTIVITY:
        print('------ANALYSIS_CONNECTIVITY...---------------------------------------')
        print('Statistics of the number of synaptic partners...')
        connectivity_statistics(GC_Objects, Cell_type='GC')
        connectivity_statistics(MF_Objects, Cell_type='MF')
        print('')

        #print('Edge analysis...')
        #Edges = extract_edges(MF_Objects)
        #Edges2 = extract_edges2(MF_Objects) #generalized
        #print('Edge type:', type(Edges), 'Edges shape:', np.shape(np.array(Edges, dtype=object))\
        #    , 'Edge[0] type:', type(Edges[0])) 
        #To use Numpy with the Lists with varying size of rows, use np.array(list, dtype=object)
        print('')

        print('------------Measuring Conenctivity Preference...----------------------------------------------')
        ratio_distribution = connectivity_ratio_distribution(MF_Objects, GC_Objects, GC_Colormap)
        mean_x, std_x= Statistics_distribution(ratio_distribution)
        print('Stats of the Reconstructed Network')
        print('mean_x:', mean_x, 'std_x:', std_x, 'var_x:', np.var(ratio_distribution))
        cumu_x, cumu_y = cumulative_distribution(ratio_distribution, print_dist=True)
        sys.exit()

        Repetition_for_stat=10
        random_mean, random_std = Statistics_Randomly_Connected_Cells(Num_MFs, Num_GCs, MF_Colormap, GC_Colormap, \
            Num_Repeat=Repetition_for_stat)
        print('Random Net Stat difference')
        print('mean difference:', mean_x-random_mean, 'std difference:', std_x-random_std )

        Replica_MF_Objects, Replica_GC_Objects = Replicate_Sample_Cells(MF_Objects, GC_Objects)
        #BEFOR SHUFFL
        ratio_distribution = connectivity_ratio_distribution(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap)
        mean_x, std_x= Statistics_distribution(ratio_distribution)
        print('Stats of the REPLICA Network before shuffling')
        print('mean_x:', mean_x, 'std_x:', std_x, 'var_x:', np.var(ratio_distribution))

        mean_shuffled, std_shuffled = Statistics_Shuffled_Cells(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap, Num_Repeat=10)      
        print('Shuffled Net Stat difference')
        #CDF_comparison_of_shuffling??
        print('mean difference:', mean_x-mean_shuffled, 'std difference:', std_x-std_shuffled )

        if PLOT_CUMU_DIST:
            cumu_x, cumu_y = cumulative_distribution(ratio_distribution, print_dist=False)
            sample_cumu_x, sample_cumu_y = cumulative_distribution(Sample_ratio_dist, print_dist=False)
            data_to_plot = [[cumu_x, cumu_y,               'Reconstructed']\
                          , [sample_cumu_x, sample_cumu_y, 'Rndmly-cnnctd']]
            plot_distributions_together(data_to_plot)
        
        



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
    print('-------Analysis ended---------------------------')
    print('Total elapsed time for analysis:', elapsed_time)  