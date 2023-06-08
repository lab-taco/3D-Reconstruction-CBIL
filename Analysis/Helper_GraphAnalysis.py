from networkx.algorithms import bipartite
import networkx as nx
import sys



def GCL_Bipartite_Graph(GC_Objects, MF_Objects, Total_Edges, DiGraph=False, Plot=False, Module_separation=False): #Color not needed?
    if DiGraph: B = nx.DiGraph() #For the order of nodes, create Digraph first, then change to unordered graph, this is for node sequencing (Ass Coef will be same)
    else: B = nx.Graph() 

    Num_edges=0
    for ind_mfbyedge, Edges_per_mf in enumerate(Total_Edges):
        #print('Edges', Edges)
        #print(ind_mfbyedge,'-th edge:', len(Edges_per_mf))
        Num_edges+=len(Edges_per_mf)
        ind_MF=Edges_per_mf[0][0]
        node_MF = "M%d"%ind_MF
        mf_color = (MF_Objects[ind_MF].color.x, MF_Objects[ind_MF].color.y, MF_Objects[ind_MF].color.z)
        #B.add_node(node_MF, bipartite=0, color=list(MF_Objects[ind_MF].color), ind=ind_MF)
        B.add_node(node_MF, bipartite=0, color=mf_color, ind=ind_MF)

        for edge in Edges_per_mf:
            ind_GC=edge[1]
            node_GC = "G%d"%ind_GC
            #node_GC = ind_GC

            #B.add_node(node_MF, bipartite=0, color='red')
            #B.add_node(node_GC, bipartite=1, color='blue')
            #B.add_node(node_GC, bipartite=1, color=tuple(GC_Objects[ind_GC].color), ind=ind_GC)
            gc_color = (GC_Objects[ind_GC].color.x, GC_Objects[ind_GC].color.y, GC_Objects[ind_GC].color.z)
            B.add_node(node_GC, bipartite=1, color=gc_color, ind=ind_GC)
            
            #print('edge', edge, 'type:', type((Ind_MF, Ind_GC)))
            #B.add_edges_from([(Ind_MF, Ind_GC)])
            B.add_edge(node_MF, node_GC)
        
    #print('Num Total edges:', B.number_of_edges())
        #print(B.edges())
    Top_Nodes = {n for n, d in B.nodes(data=True) if d["bipartite"] == 0}
    Btm_Nodes = set(B) - Top_Nodes
    #Top_Nodes, Btm_Nodes = bipartite.sets(B)    #Does not work when the network not connected
    #print('Check nodes:', len(Top_Nodes), len(Btm_Nodes))

    if Plot: Plot_Bipartite_Graph(B, Top_Nodes, Btm_Nodes, Module_separation)

    if DiGraph: 
        print('Network Constructed, changing to undirected.....')
        B=B.to_undirected()
    return B, Top_Nodes, Btm_Nodes

def Degree_Assortative_Mixing(Graph, Projection_Nodes, Mode_numeric=True, numeric_attribute="ind", static_attribute="color"):
    #One-mode projection
    Graph_onemode_projected= bipartite.projected_graph(Graph, Projection_Nodes)
    #Assortative coefficient
    #print('average_degree_connectivity:', nx.average_degree_connectivity(Graph_onemode_projected))
    if Mode_numeric: 
        return nx.numeric_assortativity_coefficient(Graph_onemode_projected, numeric_attribute)
    else: 
        return nx.attribute_assortativity_coefficient(Graph_onemode_projected, static_attribute)



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
    

def Two_module_edge_swapping(TM_B, Module_size_GC, Print_Analysis=False):
    if Print_Analysis: Copygraph=TM_B.copy()

    # MF's bipartite=0 (left nodes)

    #Swap
    """
    1. To pick edge until the swap candidates are unique
    while True:    
        select_edge1, select_edge2, Swapped_edge1, Swapped_edge2 = Pick_edges_to_swap_Twomodule(TM_B, Module_size_GC, Print_Analysis)
        Swapped_edge1_duplicate, Swapped_edge2_duplicate = Swapped_edge1 in TM_B.edges, Swapped_edge2 in TM_B.edges
        if not Swapped_edge1_duplicate and Swapped_edge2_duplicate:
            break
    TM_B.remove_edges_from([select_edge1, select_edge2])
    TM_B.add_edges_from([Swapped_edge1, Swapped_edge2])
    """

    """
    2. To skip swapping if the swap candidate is duplicated
    """
    select_edge1, select_edge2, Swapped_edge1, Swapped_edge2 = Pick_edges_to_swap_Twomodule(TM_B, Module_size_GC, Print_Analysis)
    Swapped_edge1_duplicate, Swapped_edge2_duplicate = Swapped_edge1 in TM_B.edges, Swapped_edge2 in TM_B.edges
    if not Swapped_edge1_duplicate:
        TM_B.remove_edges_from([select_edge1])
        TM_B.add_edges_from([Swapped_edge1])
    if not Swapped_edge2_duplicate:
        TM_B.remove_edges_from([select_edge2])
        TM_B.add_edges_from([Swapped_edge2])

    #if [Swapped_edge1, Swapped_edge2] not in TM_B.edges:
    #    print("target edges1 After swap:", [edge for edge in TM_B.edges if select_gc1 in edge])
    #    print("target edges2 After swap:", [edge for edge in TM_B.edges if select_gc2 in edge])
    #    raise Exception ('Swapped edges not added')
    #Check if two graph equal
    if Print_Analysis:
        print("is_isomorphic:", nx.is_isomorphic(Copygraph, TM_B))
        #print("target edges1 After swap:", [edge for edge in TM_B.edges if select_gc1 in edge])
        #print("target edges2 After swap:", [edge for edge in TM_B.edges if select_gc2 in edge])
    return TM_B

def Pick_edges_to_swap_Twomodule(TM_B, Module_size_GC, Print_Analysis):
    # MF's bipartite=0 (left nodes)
    #Select a GC node in each module randomly
    ind_select_gc1, ind_select_gc2 = np.random.choice(Module_size_GC, 2)
    ind_select_gc2+=Module_size_GC
    select_gc1, select_gc2 = 'G%d'%ind_select_gc1, 'G%d'%ind_select_gc2

    GC1_exist, GC2_exist = select_gc1 in TM_B.nodes, select_gc2 in TM_B.nodes
    if not (GC1_exist and GC2_exist):
        raise Exception ('GC1_exist:', GC1_exist, "GC2_exist:", GC2_exist, "Selected GCs do not exist, During Two_module_edge_swapping")

    #print(TM_B.edges)

    #Select an edge in each GC module
    target_edges1 = [edge for edge in TM_B.edges if select_gc1 in edge]
    target_edges2 = [edge for edge in TM_B.edges if select_gc2 in edge]
        
    if not (len(target_edges1) and len(target_edges2)):
        raise Exception ('Len GC1 edges:', len(target_edges1), "Len GC2 edges:", len(target_edges2), "GC edges are not same, maybe check if not 4")
    
    ind_select_edge1 = np.random.choice(len(target_edges1))
    ind_select_edge2 = np.random.choice(len(target_edges2))
    select_edge1, select_edge2 = target_edges1[ind_select_edge1], target_edges2[ind_select_edge2]

    target_mf1 = select_edge1[0] if select_edge1[0]!=select_gc1 else select_edge1[1]
    target_mf2 = select_edge2[0] if select_edge2[0]!=select_gc2 else select_edge2[1]

    MF1_exist, MF2_exist = target_mf1 in TM_B.nodes, target_mf2 in TM_B.nodes
    if not (MF1_exist and MF2_exist):
        raise Exception ('MF1_exist:', MF1_exist, "MF2_exist:", MF2_exist, "Selected MFs do not exist, During Two_module_edge_swapping")
    
    Swapped_edge1=(target_mf2, select_gc1)
    Swapped_edge2=(target_mf1, select_gc2)

    if Print_Analysis:
        print('Selected GCs:', select_gc1, select_gc2)
        print('Target edges1:', target_edges1, len(target_edges1))
        print('Target edges2:', target_edges2, len(target_edges2))
        print('Selected edges:', select_edge1, select_edge2)
        print('Swapped edges:', Swapped_edge1, Swapped_edge2)
    
    return select_edge1, select_edge2, Swapped_edge1, Swapped_edge2
    


import matplotlib.pyplot as plt
def Plot_Bipartite_Graph(Graphtoplot, TopNodes, BtmNodes, Module_separation=False):
    print('Drawing Networks......... (Bipartite) ')

    #print(TopNodes, BtmNodes)
    
    pos = dict()
    if not Module_separation: 
        pos.update( (n, (1, Graphtoplot.nodes[n]["ind"])) for n in TopNodes ) # put nodes from TopNodes at x=1
    else: 
        Separate_space=5
        pos.update( (n, (1, Graphtoplot.nodes[n]["ind"])) for n in TopNodes if Graphtoplot.nodes[n]["ind"] < len(TopNodes)/2)
        pos.update( (n, (1, Graphtoplot.nodes[n]["ind"] + Separate_space)) for n in TopNodes if Graphtoplot.nodes[n]["ind"] >= len(TopNodes)/2)

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
