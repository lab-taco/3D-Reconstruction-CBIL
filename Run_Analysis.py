#from Simulation.form_synapse import *
import imp
from Simulation.Data_Extraction import *

from Analysis.Analysis_Connectivity import *
from Analysis.Analysis_Spatial_dist import *
from Analysis.Slicing import *
from Analysis.Analysis_K_Function import *
from Simulation.Parameters import *
from Analysis.Helper_shapely import *
from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt

DEFEAULT_GC_object_data= 'save_cell_objects_Num_parents10_time13-05-2022-1847_GCs'
DEFEAULT_MF_object_data= 'save_cell_objects_Num_parents10_time13-05-2022-1847_MFs'

#DATA_FOLDER='14-05-2022-1333' #simple model
DATA_FOLDER ='16-05-2022-2348' # large population

PLOT_BASELINE=False
PLOT_SPATIAL_DIST=True
PLOT_CUMU_DIST=False
PLOT_CONVEX_HULL = False
K_ANALYSIS=True
NND_ANALYSIS=False
def main(GC_data_name, MF_data_name, ANALYSIS_SPATIAL_DISTRIBUTION, ANALYSIS_CONNECTIVITY):

    print('-------Data Loading...----------------------------------------------------')
    GC_Objects, MF_Objects, GC_Colormap, MF_Colormap = load_all(DATA_FOLDER, DATA_PATH, Print=True)
    Num_GCs=len(GC_Objects)
    Num_MFs=len(MF_Objects)
    #len_synapse(GCs, MFs) -------------------------------------------
    
    ANALYSIS_SPATIAL_DISTRIBUTION=True
    if ANALYSIS_SPATIAL_DISTRIBUTION:
        Map_size_3D= [area_length, area_width, height_PCL]
        Map_size_2D= [area_length, height_PCL]
        Cell_radii=[radius_MFR, radius_GC]        
        print('------ANALYSIS_SPATIAL_DISTRIBUTION...-------------------------------')
        #data_extraction_3d(GC_Objects, Map_size_2D, radius_GC, PLOTTING=True)
        #data_extraction_2d(GC_Objects, Map_size_2D, radius_GC, PLOTTING=True)
        View_3D_Dist(MF_Objects, GC_Objects)
        print('Density of Cells in Total Volume:', \
            (len(GC_Objects)+len(MF_Objects))/(area_length*area_width*height_PCL))

        MF_Captured, GC_Captured = Capture_Cells_at_Slice(MF_Objects, GC_Objects, \
                        Map_size_3D, Cell_radii, \
                        PLOTTING=PLOT_SPATIAL_DIST, Where_to_Slice='Random')
        #sys.exit()
        Position_MFs=np.asarray([[mf.body.pos.x, mf.body.pos.y] for mf in MF_Captured])
        Position_GCs=np.asarray([[gc.body.pos.x, gc.body.pos.y] for gc in GC_Captured])
        
        
        Num_cell=len(GC_Captured)
        Volume=Map_size_2D[0]*Map_size_2D[1]
        Density=Num_cell/Volume
        GRIDs=int(math.sqrt(Density*Volume))
        print('Base Cell Density:', Density, 'Volume',Volume, 'GRIDs', GRIDs)
        random_dist = random_scattering(Num_cell, Map_size_2D, plotting=PLOT_BASELINE)
        regular_dist = regular_scattering(num_grid=GRIDs, Map_size_2D=Map_size_2D, plotting=PLOT_BASELINE)
        if PLOT_BASELINE:
            plt.legend()
            plt.title('Baseline Cell Distributions')                                          
            plt.show()
            
        
        if NND_ANALYSIS:
            GC_NND, rnge_NNDgc = nearest_neighbor_distance_distribution(Position_GCs , Map_size_2D[0])
            MF_NND, rnge_NNDmf = nearest_neighbor_distance_distribution(Position_MFs , Map_size_2D[0])
            Rand_NND, rnge_NNDrnd = nearest_neighbor_distance_distribution(random_dist , Map_size_2D[0])
            Reglr_NND, rnge_NNDrglr = nearest_neighbor_distance_distribution(regular_dist , Map_size_2D[0])

            NND_List = [[GC_NND, rnge_NNDgc, 'GCs'], [MF_NND, rnge_NNDmf, 'MFs'], \
                        [Rand_NND, rnge_NNDrnd, 'Rand'], [Reglr_NND, rnge_NNDrglr, 'Reglr']]
            Plot_NND_all_together(NND_List)
            sys.exit()
            Plot_NND(GC_NND, range_NND, plotting_type='cdf', plot_label='GCs')
            Plot_NND(MF_NND, range_NND, plotting_type='cdf', plot_label='MFs')
            Plot_NND(Rand_NND, range_NND, plotting_type='cdf', plot_label='Rand')
            Plot_NND(Reglr_NND, range_NND, plotting_type='cdf', plot_label='Reglr')
            plt.title('Nearest Neighbor Distance, CDF')
            plt.legend()                                          
            plt.show()

            Plot_NND(GC_NND, range_NND, plotting_type='hist', plot_label='GCs')
            Plot_NND(MF_NND, range_NND, plotting_type='hist', plot_label='MFs')
            Plot_NND(Rand_NND, range_NND, plotting_type='hist', plot_label='Rand')
            Plot_NND(Reglr_NND, range_NND, plotting_type='hist', plot_label='Reglr')
            plt.title('Nearest Neighbor Distance, Hist')
            plt.legend()                                          
            plt.show()

        if K_ANALYSIS:
            Plot_Convex_Hull=True
            if Plot_Convex_Hull:
                MF_CvH=point_dist_to_convex_hull(Position_MFs, plotting=PLOT_CONVEX_HULL, plot_label='CvxH MFs')        
                GC_CvH=point_dist_to_convex_hull(Position_GCs, plotting=PLOT_CONVEX_HULL, plot_label='CvxH GCs')
                RD_CvH=point_dist_to_convex_hull(random_dist, plotting=PLOT_CONVEX_HULL, plot_label='CvxH Rand')
                RG_CvH=point_dist_to_convex_hull(regular_dist, plotting=PLOT_CONVEX_HULL, plot_label='CvxH Reglr')
                plt.scatter(Position_MFs[:, 0], Position_MFs[:, 1], label='MFs')
                plt.scatter(Position_GCs[:, 0], Position_GCs[:, 1], label='GCs')
                rd=np.array(random_dist)
                rg=np.array(regular_dist)
                plt.scatter(rd[:,0], rd[:,1], label='Rand dist., Num Cells:'+str(len(random_dist)))
                plt.scatter(rg[:,0], rg[:,1], label='Reglr dist., Num Cells:'+str(len(regular_dist)))
                
                plt.title('Convex Hulls')
                plt.legend()
                plt.show()   

            FUNCTION='L'
            plotting_each_Func=False            
            L_random  = my_K_func(random_dist, Map_size_2D, radius_GC, \
                                function_type=FUNCTION,graph=plotting_each_Func, return_L=True)
            L_regular = my_K_func(regular_dist, Map_size_2D, radius_GC, \
                                function_type=FUNCTION, graph=plotting_each_Func, return_L=True)
            L_MFs = my_K_func(Position_MFs, Map_size_2D, radius_GC, \
                                function_type=FUNCTION, graph=plotting_each_Func, return_L=True)
            L_GCs = my_K_func(Position_GCs, Map_size_2D, radius_GC, \
                                function_type=FUNCTION, graph=plotting_each_Func, return_L=True)
            plt.plot(L_random[:,0], L_random[:,1], color='k', label='random dist')
            plt.plot(L_regular[:,0], L_regular[:,1], color='y', label='regular dist')
            plt.plot(L_MFs[:,0], L_MFs[:,1], color='b', label='MFs')
            plt.plot(L_GCs[:,0], L_GCs[:,1], color='red', label='GCs')
            plt.plot(L_GCs[:,0], np.zeros(len(L_GCs[:,0])), color='c', ls=':', label=r'$L_{pois}$') 
            plt.title('Spatial Analysis using L function')
            plt.legend()
            plt.show()
            
        
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
                        default='DEFEAULT_GC_object_data',
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