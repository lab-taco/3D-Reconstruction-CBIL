class Sample_cells_with_synapses: 
    def __init__(self, cell_type, color_):
        self.synapse_partners=[] #as indices
        self.color=color_

        if cell_type=='MFR':
            self.cell_type='MFR'
        elif cell_type=='GC':
            self.cell_type='GC'

def randomly_conencted_sample_cells(Num_MFs, Num_GCs, colormap_mf, colormap_gc):    
    Sample_MF_Objects=[]
    for _ in range(Num_MFs):
        sample_mf=Sample_cells_with_synapses('MF', np.random.choice(colormap_mf))        
        Sample_MF_Objects.append(sample_mf)
        
    Sample_GC_Objects=[]    
    for ind_gc in range(Num_GCs):
        sample_gc=Sample_cells_with_synapses('GC', np.random.choice(colormap_gc))
        synaptic_partners=np.random.choice(Num_MFs, size=4, replace=False)
        sample_gc.synapse_partners=sorted(synaptic_partners)
        Sample_GC_Objects.append(sample_gc)
        for ind_sp in synaptic_partners:
            Sample_MF_Objects[ind_sp].synapse_partners.append(ind_gc)
   
    return Sample_MF_Objects, Sample_GC_Objects