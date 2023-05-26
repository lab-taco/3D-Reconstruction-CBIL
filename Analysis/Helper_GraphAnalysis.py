from networkx.algorithms import bipartite
import networkx as nx
import sys



def GCL_Bipartite_Graph(GC_Objects, MF_Objects, Total_Edges, Plot=False): #Color not needed?
    B = nx.DiGraph() #For the order of nodes, create Digraph first, then change to unordered graph

    Num_edges=0
    for Edges_per_mf in Total_Edges:
        #print('Edges', Edges)
        Num_edges+=len(Edges_per_mf)
        ind_MF=Edges_per_mf[0][0]
        node_MF = "M%d"%ind_MF 
        B.add_node(node_MF, bipartite=0, color=MF_Objects[ind_MF], ind=ind_MF)

        for edge in Edges_per_mf:
            ind_GC=edge[1]
            node_GC = "G%d"%ind_GC
            #node_GC = ind_GC

            #B.add_node(node_MF, bipartite=0, color='red')
            #B.add_node(node_GC, bipartite=1, color='blue')
            B.add_node(node_GC, bipartite=1, color=GC_Objects[ind_GC], ind=ind_GC)
            
            #print('edge', edge, 'type:', type((Ind_MF, Ind_GC)))
            #B.add_edges_from([(Ind_MF, Ind_GC)])
            B.add_edge(node_MF, node_GC)
        
        #print(B.number_of_edges())
        #print(B.edges())
        #break
    Top_Nodes = {n for n, d in B.nodes(data=True) if d["bipartite"] == 0}
    Btm_Nodes = set(B) - Top_Nodes
    #Top_Nodes, Btm_Nodes = bipartite.sets(B)    #Does not work when the network not connected

    if Plot:
        Plot_Bipartite_Graph(B, Top_Nodes, Btm_Nodes)
    return B, Top_Nodes, Btm_Nodes

def Degree_Assortative_Mixing(Graph, Projection_Nodes, numeric_attribute="ind"):
    #One-mode projection
    Graph_onemode_projected= bipartite.projected_graph(Graph, Projection_Nodes)
    #Assortative coefficient
    Assr_coeff_MFs_attribute = nx.numeric_assortativity_coefficient(Graph_onemode_projected, numeric_attribute)
    return Assr_coeff_MFs_attribute



from .Sample_Cell_Gen import *
from Analysis.Analysis_Connectivity import extract_edges2
def Two_module_network(Num_MFs, Num_GCs, GC_Colormap, MF_Colormap):

    Module_size_MF, Module_size_GC = int(Num_MFs/2), int(Num_GCs/2) #Ind_MF2_start, Ind_GC2_start
    MF1, GC1=randomly_conencted_sample_cells(int(Num_MFs/2), int(Num_GCs/2), MF_Colormap, [GC_Colormap[0]])
    MF2, GC2=randomly_conencted_sample_cells(int(Num_MFs/2), int(Num_GCs/2), MF_Colormap, [GC_Colormap[-1]]\
                                            , add_MF=Module_size_MF, add_GC= Module_size_GC)

    MFs = MF1+MF2
    GCs = GC1+GC2
    MF_Edges= extract_edges2(MFs)
    
    return MFs, GCs, MF_Edges, Module_size_MF, Module_size_GC
    

def Two_module_edge_swapping(TM_B, Module_size_GC):
    # MF's bipartite=0 (left nodes)
    Copygraph=TM_B.copy()
    #Select a GC node in each module randomly
    ind_select_gc1, ind_select_gc2 = np.random.choice(Module_size_GC, 2)
    ind_select_gc2+=Module_size_GC
    select_gc1, select_gc2 = 'G%d'%ind_select_gc1, 'G%d'%ind_select_gc2
    print('Selected GCs:', select_gc1, select_gc2)
    

    #print(TM_B.edges)

    #Select an edge in each GC
    target_edges1 = [edge for edge in TM_B.edges if select_gc1 in edge]
    target_edges2 = [edge for edge in TM_B.edges if select_gc2 in edge]
    print('Target edges1:', target_edges1, len(target_edges1))
    print('Target edges2:', target_edges2, len(target_edges2))

    ind_select_edge1 = np.random.choice(len(target_edges1))
    ind_select_edge2 = np.random.choice(len(target_edges2))
    select_edge1, select_edge2 = target_edges1[ind_select_edge1], target_edges2[ind_select_edge2]
    print('Selected edges:', select_edge1, select_edge2)

    #target_mf1 = select_edge1[0]
    #target_mf2 = select_edge2[0]
    target_mf1 = [node for node in select_edge1 if node!=select_gc1][0]
    target_mf2 = [node for node in select_edge2 if node!=select_gc2][0]
    


    #Swap
    Swapped_edge1=(target_mf2, select_gc1)
    Swapped_edge2=(target_mf1, select_gc2)
    print('Swapped edges:', select_edge1, select_edge2)
    TM_B.remove_edges_from([select_edge1, select_edge2])
    TM_B.add_edges_from([Swapped_edge1, Swapped_edge2])
    #Check if two graph equal
    print(nx.is_isomorphic(Copygraph, TM_B))
    return TM_B





import matplotlib.pyplot as plt
def Plot_Bipartite_Graph(Graphtoplot, TopNodes, BtmNodes):
    print('Drawing Networks......... (Bipartite) ')

    print(TopNodes, BtmNodes)
    
    pos = dict()
    pos.update( (n, (1, Graphtoplot.nodes[n]["ind"])) for n in TopNodes ) # put nodes from TopNodes at x=1
    pos.update( (n, (2, Graphtoplot.nodes[n]["ind"])) for n in BtmNodes ) # put nodes from BottomNodes at x=2

    #pos = nx.bipartite_layout(Graph, TopNodes)
    nx.draw_networkx_nodes (Graphtoplot, pos)
    nx.draw_networkx_edges (Graphtoplot, pos, Graphtoplot.edges)
    nx.draw_networkx_labels(Graphtoplot, pos)
    nx.draw_networkx_nodes (Graphtoplot, pos, nodelist=TopNodes, node_color="tab:red")
    plt.show()

def graph_plotting():

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
