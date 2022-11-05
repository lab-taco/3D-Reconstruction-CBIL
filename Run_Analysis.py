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
import pandas as pd
import openpyxl
from Analysis.Modularity import *

DEFEAULT_GC_object_data= 'save_cell_objects_Num_parents10_time13-05-2022-1847_GCs'
DEFEAULT_MF_object_data= 'save_cell_objects_Num_parents10_time13-05-2022-1847_MFs'

#DATA_FOLDER='14-05-2022-1333' #simple model
#DATA_FOLDER ='16-05-2022-2348' # large population
DATA_FOLDER='SP0'

PLOT_Connectivity_CDF=True
RAND_NET=False
SHUFFL=False

K_ANALYSIS=False
NND_ANALYSIS=False
PLOT_BASELINE=False
PLOT_CONVEX_HULL = False
PLOT_SPATIAL_DIST=True
def main(GC_data_name, MF_data_name, ANALYSIS_SPATIAL_DISTRIBUTION, ANALYSIS_CONNECTIVITY):

    print('-------Data Loading...----------------------------------------------------')
    GC_Objects, MF_Objects, GC_Colormap, MF_Colormap = load_all(DATA_FOLDER, DATA_PATH, Print=True)
    Num_GCs=len(GC_Objects)
    Num_MFs=len(MF_Objects)
    #len_synapse(GCs, MFs) -------------------------------------------
    
    print('Num GCs:', Num_GCs, 'Num MFs:', Num_MFs)
    #GC_synaptic_partner_exchange(MF_Objects, GC_Objects, GC_Colormap, Size=3, D=1)
    
    D=1
    #Num_MFs=10
    #overlapped_network(Num_MFs, Num_MFs*3, MF_Colormap, GC_Colormap, 50, D)
    #sys.exit()
    ANALYSIS_SPATIAL_DISTRIBUTION=False
    if ANALYSIS_SPATIAL_DISTRIBUTION:
        Map_size_3D= [area_length, area_width, height_PCL]
        Map_size_2D= [area_length, height_PCL]
        Cell_radii=[radius_MFR, radius_GC]        
        print('------ANALYSIS_SPATIAL_DISTRIBUTION...-------------------------------')
        #data_extraction_3d(GC_Objects, Map_size_2D, radius_GC, PLOTTING=True)
        #data_extraction_2d(GC_Objects, Map_size_2D, radius_GC, PLOTTING=True)
        #View_3D_Dist(MF_Objects, GC_Objects)
        print('Density of Cells in Total Volume:', \
            (len(GC_Objects)+len(MF_Objects))/(area_length*area_width*height_PCL))

        
        #MF_Captured, GC_Captured = Capture_Cells_at_Slice(MF_Objects, GC_Objects, \
        #                Map_size_3D, Cell_radii, \
        #                PLOTTING=PLOT_SPATIAL_DIST, Where_to_Slice='Random')

        MF_Captured, E_Captured, M_Captured, L_Captured = Capture_Cells_at_Slice_group_separation(MF_Objects, GC_Objects, \
                        GC_Colormap, Map_size_3D, Cell_radii, \
                        PLOTTING=PLOT_SPATIAL_DIST, Where_to_Slice='Random')

        GC_Captured =E_Captured+ M_Captured+ L_Captured
        #sys.exit()
        Position_MFs=np.asarray([[mf.body.pos.x, mf.body.pos.y] for mf in MF_Captured])
        Position_GCs=np.asarray([[gc.body.pos.x, gc.body.pos.y] for gc in GC_Captured])
        Position_E=np.asarray([[gc.body.pos.x, gc.body.pos.y] for gc in E_Captured])
        Position_M=np.asarray([[gc.body.pos.x, gc.body.pos.y] for gc in M_Captured])
        Position_L=np.asarray([[gc.body.pos.x, gc.body.pos.y] for gc in L_Captured])
        
        
        Num_cell=len(Position_E)
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
            EGC_NND, rnge_NNDE = nearest_neighbor_distance_distribution(Position_E , Map_size_2D[0])
            MGC_NND, rnge_NNDG = nearest_neighbor_distance_distribution(Position_M , Map_size_2D[0])
            LGC_NND, rnge_NNDL = nearest_neighbor_distance_distribution(Position_L , Map_size_2D[0])
            MF_NND, rnge_NNDmf = nearest_neighbor_distance_distribution(Position_MFs , Map_size_2D[0])
            Rand_NND, rnge_NNDrnd = nearest_neighbor_distance_distribution(random_dist , Map_size_2D[0])
            Reglr_NND, rnge_NNDrglr = nearest_neighbor_distance_distribution(regular_dist , Map_size_2D[0])

            NND_List = [[GC_NND, rnge_NNDgc, 'GCs'], [MF_NND, rnge_NNDmf, 'MFs'], \
                         [EGC_NND, rnge_NNDE, 'Early GCs'], [MGC_NND, rnge_NNDG, 'Mid GCs'], \
                            [LGC_NND, rnge_NNDL, 'Late GCs'], [MF_NND, rnge_NNDmf, 'MFs'], \
                        [Rand_NND, rnge_NNDrnd, 'Rand'], [Reglr_NND, rnge_NNDrglr, 'Reglr']]
            Plot_NND_all_together(NND_List)
            #sys.exit()
            """
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
            plt.show()"""

        if K_ANALYSIS:
            Plot_Convex_Hull=True
            if Plot_Convex_Hull:
                
                MF_CvH=point_dist_to_convex_hull(Position_MFs, plotting=PLOT_CONVEX_HULL, plot_label='CvxH MFs')        
                GC_CvH=point_dist_to_convex_hull(Position_GCs, plotting=PLOT_CONVEX_HULL, plot_label='CvxH GCs')
                E_CvH=point_dist_to_convex_hull(Position_E, plotting=PLOT_CONVEX_HULL, plot_label='CvxH GCs')
                M_CvH=point_dist_to_convex_hull(Position_M, plotting=PLOT_CONVEX_HULL, plot_label='CvxH GCs')
                L_CvH=point_dist_to_convex_hull(Position_L, plotting=PLOT_CONVEX_HULL, plot_label='CvxH GCs')
                RD_CvH=point_dist_to_convex_hull(random_dist, plotting=PLOT_CONVEX_HULL, plot_label='CvxH Rand')
                RG_CvH=point_dist_to_convex_hull(regular_dist, plotting=PLOT_CONVEX_HULL, plot_label='CvxH Reglr')
                plt.scatter(Position_MFs[:, 0], Position_MFs[:, 1], label='MFs')
                plt.scatter(Position_GCs[:, 0], Position_GCs[:, 1], label='GCs')
                plt.scatter(Position_E[:, 0], Position_E[:, 1], label='GCs E')
                plt.scatter(Position_M[:, 0], Position_M[:, 1], label='GCs M')
                plt.scatter(Position_L[:, 0], Position_L[:, 1], label='GCs L')
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
            L_GCsE = my_K_func(Position_E, Map_size_2D, radius_GC, \
                                function_type=FUNCTION, graph=plotting_each_Func, return_L=True)
            L_GCsM = my_K_func(Position_M, Map_size_2D, radius_GC, \
                                function_type=FUNCTION, graph=plotting_each_Func, return_L=True)
            L_GCsL = my_K_func(Position_L, Map_size_2D, radius_GC, \
                                function_type=FUNCTION, graph=plotting_each_Func, return_L=True)
            plt.plot(L_random[:,0], L_random[:,1], color='k', label='random dist')
            plt.plot(L_regular[:,0], L_regular[:,1], color='y', label='regular dist')
            plt.plot(L_MFs[:,0], L_MFs[:,1], color='b', label='MFs')
            plt.plot(L_GCs[:,0], L_GCs[:,1], color='red', label='GCs')
            plt.plot(L_GCsE[:,0], L_GCsE[:,1], label='GC E')
            plt.plot(L_GCsM[:,0], L_GCsM[:,1], label='GC M')
            plt.plot(L_GCsL[:,0], L_GCsL[:,1], label='GC L')
            plt.plot(L_GCs[:,0], np.zeros(len(L_GCs[:,0])), color='c', ls=':', label=r'$L_{pois}$') 
            plt.title('Spatial Analysis using L function')
            plt.legend()
            plt.show()
            
    if ANALYSIS_CONNECTIVITY:
        print('------ANALYSIS_CONNECTIVITY...---------------------------------------')
        print('Statistics of the number of synaptic partners...')
        basic_connectivity_statistics(GC_Objects, Cell_type='GC')
        basic_connectivity_statistics(MF_Objects, Cell_type='MF')
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
        data_mean_x, data_std_x= Statistics_distribution(ratio_distribution)
        print('Stats of the Reconstructed Network', DATA_FOLDER)
        print('mean_x:', data_mean_x, 'std_x:', data_std_x, 'var_x:', np.var(ratio_distribution))
        
        #cumu_x, cumu_y = cumulative_distribution(ratio_distribution, \
        #    label_cdf='GC-MF Ratio', print_dist=True)
        #sys.exit()

        if RAND_NET:
            Repetition_for_stat=100
            random_mean, random_std, Random_Sample_ratio_dist\
                 = Statistics_Randomly_Connected_Cells(Num_MFs, Num_GCs, \
                                                        MF_Colormap, GC_Colormap, \
                                                        Num_Repeat=Repetition_for_stat)
            print('Random Net Stat difference')
            print('mean:', random_mean, 'std:', random_std )
            print('mean difference:', data_mean_x-random_mean, 'std difference:', data_std_x-random_std )
        #sys.exit()
        if SHUFFL:
            Replica_MF_Objects, Replica_GC_Objects = Replicate_Sample_Cells(MF_Objects, GC_Objects)
            #BEFOR SHUFFL
            Replica_ratio_dist = connectivity_ratio_distribution(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap)
            repl_mean_x, repl_std_x= Statistics_distribution(Replica_ratio_dist)
            if not [repl_mean_x, repl_std_x] == [data_mean_x, data_std_x]:
                raise Exception('Replication failed')
            #print('Stats of the REPLICA Network before shuffling')
            #print('mean_x:', mean_x, 'std_x:', std_x, 'var_x:', np.var(ratio_distribution))
            #sys.exit()
            #AFTER SHUFFL
            SIZE_EXCHANGE = 10 
            D_ROUNDS=2
            ratio_distribution_shuffled, mean_shuffled, std_shuffled = \
                Statistics_Shuffled_Cells(Replica_MF_Objects, Replica_GC_Objects, GC_Colormap,\
                                 exchange_size = SIZE_EXCHANGE, D_rounds=D_ROUNDS, Num_Repeat=10)      
            print('Shuffled Net Stat difference')
            #CDF_comparison_of_shuffling??
            print('mean difference:', data_mean_x-mean_shuffled, 'std difference:', data_std_x-std_shuffled )

        if PLOT_Connectivity_CDF:
            cumu_x, cumu_y = cumulative_distribution(ratio_distribution, label_cdf='GC-MF Ratio', print_dist=False)

            Sample_MFs, Sample_GCs = randomly_conencted_sample_cells(Num_MFs, Num_GCs, MF_Colormap, GC_Colormap)
            Random_Sample_ratio_dist=connectivity_ratio_distribution(Sample_MFs, Sample_GCs, GC_Colormap)
            sample_cumu_x, sample_cumu_y = cumulative_distribution(Random_Sample_ratio_dist, print_dist=False)
            
            print('len cumu_x', len(cumu_x))
            #Shuffled_cumu_x, Shuffled_cumu_y = cumulative_distribution(ratio_distribution_shuffled, print_dist=False)
            data_to_plot = [[cumu_x, cumu_y,               'Reconstructed'],
                            [sample_cumu_x, sample_cumu_y, 'Rndmly-cnnctd']] #,
                            #[Shuffled_cumu_x, Shuffled_cumu_y, 'Shuffled']]

            plot_distributions_together(data_to_plot)
            sys.exit()
            print('np.shape(cumu_x)',np.shape(cumu_x), \
                'np.shape(cumu_y)',np.shape(cumu_y), type(cumu_x))
            Concat=np.vstack((cumu_x, cumu_y)).T
            print('np.concat(cumu_x)',np.shape(Concat))
            #print(Concat)
            #sys.exit()
            #df = pd.DataFrame(Concat, columns = ['x','y'])
            #df.to_excel('Connectivity_Reconstructed.xlsx')
            
        
        



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