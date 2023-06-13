import numpy as np

class Sample_cells_with_synapses: #Only have synaptic partner indices?
    def __init__(self, cell_type, color_):
        self.synapse_partners=[] #as indices
        self.color=color_

        if cell_type=='MF':
            self.cell_type='MF'
        elif cell_type=='GC':
            self.cell_type='GC'

def coloring_by_indexing(color_list, cells):
    color_count = len(color_list)
    cell_count = len(cells)
    
    for i, attribute in enumerate(color_list):
        attribute_proportion = cell_count // color_count
        attribute_population = cells[i * attribute_proportion : (i+1) * attribute_proportion]
        
        if i == color_count - 1:
            attribute_population.extend(cells[(i+1) * attribute_proportion:])
        
        for cell in attribute_population:
            cell.color==attribute
    

def randomly_conencted_sample_cells(Num_MFs, Num_GCs, colormap_mf, colormap_gc, \
                                    add_MF=0, add_GC=0):    
    Sample_MF_Objects=[]
    for _ in range(Num_MFs):
        #sample_mf=Sample_cells_with_synapses('MF', np.random.choice(colormap_mf))
        sample_mf=Sample_cells_with_synapses('MF', colormap_mf[0])
        Sample_MF_Objects.append(sample_mf)

    Sample_GC_Objects=[]
    for ind_gc in range(Num_GCs):
        #sample_gc=Sample_cells_with_synapses('GC', np.random.choice(colormap_gc))
        sample_gc=Sample_cells_with_synapses('GC', colormap_gc[0])
        synaptic_partners=np.random.choice(Num_MFs, size=4, replace=False)
        sample_gc.synapse_partners=sorted(synaptic_partners+add_MF)
        Sample_GC_Objects.append(sample_gc)
        for ind_sp in synaptic_partners:
            Sample_MF_Objects[ind_sp].synapse_partners.append(ind_gc+add_GC)
    
    coloring_by_indexing(colormap_mf, Sample_MF_Objects)
    coloring_by_indexing(colormap_gc, Sample_GC_Objects)
    return Sample_MF_Objects, Sample_GC_Objects

def Cell_Gen_by_given_edges(edges, Num_GCs, colormap_mf, colormap_gc):
    Sample_GC_Objects=[]
    for ind in range(Num_GCs):
        sample_gc=Sample_cells_with_synapses('GC', colormap_gc[ind])
        Sample_GC_Objects.append(sample_gc)

    Sample_MF_Objects=[]
    Num_MFs=len(edges)
    for ind_mf in range(Num_MFs):
        sample_mf=Sample_cells_with_synapses('MF', colormap_mf[ind_mf])
        sample_mf.synapse_partners=edges[ind_mf]        
        Sample_MF_Objects.append(sample_mf)
        for ind_sp in sample_mf.synaptic_partners:
            Sample_GC_Objects[ind_sp].synapse_partners.append(ind_mf)
        
    Sample_GC_Objects=[] 

def Replicate_Sample_Cells(MFs, GCs):
    Sample_GC_Objects=[]
    for gc in GCs:
        sample_gc=Sample_cells_with_synapses('GC', gc.color)
        sample_gc.synapse_partners=gc.synapse_partners
        Sample_GC_Objects.append(sample_gc)

    Sample_MF_Objects=[]
    for mf in MFs:
        sample_mf=Sample_cells_with_synapses('MF', mf.color)
        sample_mf.synapse_partners=mf.synapse_partners
        Sample_MF_Objects.append(sample_mf)
        
    return Sample_MF_Objects, Sample_GC_Objects