#from Simulation.form_synapse import *
from Simulation.Data_Extraction import *

from Simulation.Parameters import *

from Analysis.Analysis_Connectivity import *
from Analysis.Analysis_Spatial_dist import *
from Analysis.Analysis_K_Function import *
from Analysis.Helper_shapely import *
from Analysis.Slicing import *
from Analysis.Modularity import *
from Analysis.Helper_GraphAnalysis import *

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import sys
from vpython import *


#DATA_FOLDER='16-05-2023-1711' #simple model
#DATA_FOLDER ='16-05-2023-1723' # large population

#DATA_FOLDER='SP100-Small'
#DATA_FOLDER='SP0'
#DATA_FOLDER='Biased/SP0'
#DATA_FOLDER='SP0_activityonly'
#DATA_FOLDER='Unbiased_SP0_old'
DATA_FOLDER='Unbaised2_SP75'


PLOT_Connectivity_CDF=True
RAND_NET=False
SHUFFL=False

K_ANALYSIS=False
NND_ANALYSIS=False
PLOT_BASELINE=False
PLOT_CONVEX_HULL = False
PLOT_SPATIAL_DIST=True


def main(ANALYSE_SPATIAL_DISTRIBUTION, ANALYSE_CONNECTIVITY):
    

    print('-------Data %s Loading...----------------------------------------------------'%DATA_FOLDER)
    GC_Objects, MF_Objects, GC_Colormap, MF_Colormap, Graph = load_all(DATA_FOLDER, DATA_PATH, Print=True)
    Num_GCs=len(GC_Objects)
    Num_MFs=len(MF_Objects)

    #Graph Data shapes
    #print(np.shape(Graph)[0])
    #print(Graph[0])
    #print(Graph[0][0])
    #print(type(Graph[0]), type(Graph[0][0]))
    #print(GC_Colormap[0], GC_Colormap[-1])
    #print(np.shape(GC_Colormap))
    #print(type(GC_Colormap))
    DATA_ANALYSIS=True

    if DATA_ANALYSIS==True:
    
    
        B = nx.Graph()

        Num_edges=0
        for ind_MF, Edges in enumerate(Graph):
            if len(Edges)==0:
                print('MF', ind_MF, "has no edges")
                continue
            #print('Edges', Edges)
            Num_edges+=len(Edges)
            #ind_MF=Edges[0][0]
            node_MF = "M%d"%ind_MF 
            B.add_node(node_MF, bipartite=0, color=MF_Objects[ind_MF].color, ind=ind_MF)

            for edge in Edges:
                ind_GC=edge[1]
                node_GC = "G%d"%ind_GC
                #node_GC = ind_GC

                #B.add_node(node_MF, bipartite=0, color='red')
                #B.add_node(node_GC, bipartite=1, color='blue')
                B.add_node(node_GC, bipartite=1, color=GC_Objects[ind_GC].color, ind=ind_GC)

                #print('edge', edge, 'type:', type((Ind_MF, Ind_GC)))
                #B.add_edges_from([(Ind_MF, Ind_GC)])
                B.add_edge(node_MF, node_GC)

            #print(B.number_of_edges())
            #print(B.edges())
            #break
        print('Graph constructed......... ->>>  B connected:', nx.is_connected(B))
        #print(B.number_of_edges())
        #print('Num_edges', Num_edges)



        Node_MFs, Node_GCs = bipartite.sets(B)
        print('len nodes - MF:', len(Node_MFs), 'GC:', len(Node_GCs))

        ##Nodes' Sorted Indices
        #Sorted_MFs=[]
        #for node_ind in [mf.ind for mf in Node_MFs]:
        #    for mf in Node_MFs:
        #        if mf['ind']==node_ind:
        #            Sorted_MFs.append(mf)
        #            break
        #
        #print(Sorted_MFs)
        #sys.exit()



        #print('node data type:', type(Node_MFs))


        print('One-mode Projection and Assortative coefficient calculation........')
        Graph_onemode_MF = bipartite.projected_graph(B, Node_MFs)
        #Assr_coeff_MFs = nx.degree_assortativity_coefficient(Graph_onemode_MF)
        #Assr_coeff_MFs_attribute = nx.attribute_assortativity_coefficient(Graph_onemode_MF, "ind")
        Assr_coeff_MFs_attribute = nx.numeric_assortativity_coefficient(Graph_onemode_MF, "ind")

        Graph_onemode_GC = bipartite.projected_graph(B, Node_GCs)
        #Assr_coeff_GCs = nx.degree_assortativity_coefficient(Graph_onemode_GC)
        #Assr_coeff_GCs_attribute = nx.attribute_assortativity_coefficient(Graph_onemode_GC, "ind")
        Assr_coeff_GCs_attribute = nx.numeric_assortativity_coefficient(Graph_onemode_GC, "ind")

        #print('Assortative Coefficient MF:',Assr_coeff_MFs, 'GC:', Assr_coeff_GCs)
        print('Assortative Coefficient Attribute MF:',Assr_coeff_MFs_attribute, 'GC:', Assr_coeff_GCs_attribute)

        """
        # GC-MF Bipartite graph
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        fig.suptitle(DATA_FOLDER)

        BM = bipartite.biadjacency_matrix(B, Node_MFs, Node_GCs)
        ax1.spy(BM, markersize=1)
        ax1.set_title('Biadjacency')

        A_MF = nx.adjacency_matrix(Graph_onemode_MF)
        ax2.spy(A_MF, markersize=1)
        ax2.set_title('MF onemode')


        A_GC = nx.adjacency_matrix(Graph_onemode_GC)
        ax3.spy(A_GC, markersize=1)
        ax3.set_title('GC onemode')

        plt.show()
        """
        """
        while True: rate(30)

        print('Drawing Networks......... ')
        pos = nx.spring_layout(Graph_onemode_MF)
        nx.draw_networkx_nodes(Graph_onemode_MF, pos, node_color="tab:red")
        nx.draw_networkx_labels(Graph_onemode_MF, pos)
        nx.draw_networkx_edges(Graph_onemode_MF, pos, Graph_onemode_MF.edges)
        plt.show()


        print('Drawing Networks......... (GC) ')
        pos = nx.spring_layout(Graph_onemode_GC)
        nx.draw_networkx_nodes(Graph_onemode_GC, pos, node_color="tab:blue")
        nx.draw_networkx_labels(Graph_onemode_GC, pos)
        nx.draw_networkx_edges(Graph_onemode_GC, pos, Graph_onemode_GC.edges)
        plt.show()
        #print(len(Node_MFs), len(Node_GCs))

        print('Drawing Networks......... (Bipartite) ')
        pos = nx.bipartite_layout(B, Node_MFs)
        nx.draw_networkx_nodes(B, pos)
        nx.draw_networkx_labels(B, pos)
        nx.draw_networkx_edges(B, pos, B.edges)
        nx.draw_networkx_nodes(B, pos, nodelist=Node_MFs, node_color="tab:red")
        plt.show()


        #B.add_nodes_from([1, 2, 3, 4], bipartite=0, color='red')
        #B.add_nodes_from(["a", "b", "c"], bipartite=1)
        #bottom_nodes, top_nodes = bipartite.sets(B)
        """
        
    #Num_GCs=10
    #Num_MFs=30
    #len_synapse(GCs, MFs) -------------------------------------------
    
    
    print('Num GCs:', Num_GCs, 'Num MFs:', Num_MFs)
    #GC_synaptic_partner_exchange(MF_Objects, GC_Objects, GC_Colormap, Size=3, D=1)
    
    
    #Num_MFs=10
    #overlapped_network(Num_MFs, Num_MFs*3, MF_Colormap, GC_Colormap, 50, D)
    
    ANALYSE_SPATIAL_DISTRIBUTION=False
    
    if ANALYSE_SPATIAL_DISTRIBUTION:
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
            
    if ANALYSE_CONNECTIVITY:
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

            print('np.shape(cumu_x)',np.shape(cumu_x), \
                'np.shape(cumu_y)',np.shape(cumu_y), type(cumu_x))
            Concat=np.vstack((cumu_x, cumu_y)).T
            print('np.concat(cumu_x)',np.shape(Concat))
            #print(Concat)
            #sys.exit()
            #df = pd.DataFrame(Concat, columns = ['x','y'])
            #df.to_excel('Connectivity_Reconstructed.xlsx')
            
            SAVE_Ratio_Dist=False
            if SAVE_Ratio_Dist:
                Name_ratio_dist='sample name'
                Data_ratio_dist=Concat
                data_save(Name_ratio_dist, Data_ratio_dist, DATA_PATH, DATA_FOLDER)
        



import argparse
import time
if __name__ == '__main__':

    start_time = time.time()
    
    parser = argparse.ArgumentParser()
    
    #boolen for argparse is tricky, So used 0,1 instead

    parser.add_argument('--Spatial_analysis', type=int, dest='ANALYSE_SPATIAL_DISTRIBUTION',
                        default=1,
                        help='Running Spatial_analysis of cell soma distribution')
    parser.add_argument('--Connectivity_analysis', type=int, dest='ANALYSE_CONNECTIVITY',
                        default=1,
                        help='Running Connectivity analysis of GC-MF synapses')    

    args = parser.parse_args()
    
    #main(args.Num_childs, args.draw_dist, args.draw_net)
    main(args.ANALYSE_SPATIAL_DISTRIBUTION, args.ANALYSE_CONNECTIVITY)

    elapsed_time = time.time() - start_time
    print('-------Analysis ended---------------------------')
    print('Total elapsed time for analysis:', elapsed_time)  