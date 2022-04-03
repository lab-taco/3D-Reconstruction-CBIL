from .Parameters import *
from vpython import *
import numpy as np
import math

from .Cell_Body_Controls import check_out_of_volume
from .Cell_Body_Controls import GC_migration_complete

color_list = [color.red, color.orange, color.yellow, color.green, color.blue, \
                hat(vector(46,43,95)), hat(vector(139, 0, 255)), color.black]

Mig_Early=hat(vector(192,192,192))
Mig_Mid=hat(vector(135, 121, 78))
Mig_Late=hat(vector(109, 135, 100))
MF_colors = [Mig_Early, Mig_Mid, Mig_Late]

def colormap(Num_Variation, Cell_Type, Plot_Graph=False):
    import seaborn as sns
    #145=green, 20=red
    if Cell_Type=='MF':
        palette = sns.dark_palette((260, 75, 60), n_colors=Num_Variation, input="husl")        
    elif Cell_Type=='GC': 
        palette = sns.diverging_palette(145, 20, s=60, n=Num_Variation, as_cmap=False)
        

    if Plot_Graph:
        import matplotlib.pyplot as plt
        sns.palplot(palette)
        #plt.plot(palette)
        plt.show()

    color_map_for_vpython=[]
    for i in range(len(palette)):
        color_i=hat(vector(palette[i][0], palette[i][1],palette[i][2]))
        #color_i=vector(palette[i][0], palette[i][1],palette[i][2])
        color_map_for_vpython.append(color_i)
    return color_map_for_vpython


def random_sample_2d(length, width, location, cell_type='MF'):        
    length=np.random.choice(length)
    width=np.random.choice(width)
    height=location
    coord=vector(length,height,width)
    return coord

def sample_MFs_for_random_position(set_MFs, max_MF,area_length,area_width, height_PCL, static=False):
    MFs=set_MFs
    MF_sample_position=[]    
    for i in range(max_MF):
        if static==True: # for sample test
            MF_sample_position.append(height_PCL/2)
        else:  # for general experiments
            MF_sample_position.append(np.random.choice(height_PCL))
    #MF_sample_position.sort(reverse=True)
    for i in range(max_MF):        
        cell_position=random_sample_2d(area_length,area_width,MF_sample_position[i])
        if i<max_MF/3:
            mf_color=MF_colors[0]
        elif i<max_MF*2/3:
            mf_color=MF_colors[1]
        else:
            mf_color=MF_colors[2]
        MFs.append(Cells('MFR', height_PCL, init_position=cell_position, color_=mf_color ))    
    return MFs


def positions_MFs_childs(Position_parents, Num_childs):
    positions_childs=[]
    while not len(positions_childs)==Num_childs:
        r=np.random.exponential(100)
        theta=np.random.choice(360)* np.pi/180
        phi=np.random.choice(360)* np.pi/180
        #in Radian
        coordinate_converted=\
            vector(r*np.sin(theta)*np.cos(phi)\
                ,  r*np.sin(theta)*np.sin(phi)\
                ,  r*cos(theta))
        if check_out_of_volume(Position_parents+coordinate_converted)==False:
            positions_childs.append(Position_parents+coordinate_converted)      

    return positions_childs


def sample_MFs_clusters(List_MFs, Num_parents, Num_childs, max_MF, Mig_Timing_Variation,\
                         area_length,area_width, height_PCL,\
                         static=False, vpython=True):
    MFs=List_MFs
    MF_sample_position=[]    
    for i in range(Num_parents):
        if static==True: # for sample test
            MF_sample_position.append(height_PCL/2)
        else:  # for general experiments
            MF_sample_position.append(np.random.choice(height_PCL))
    #MF_sample_position.sort(reverse=True)
    colormap_for_mig_time=colormap(Mig_Timing_Variation, Cell_Type='MF')
    #print('len(colormap_for_mig_time)',len(colormap_for_mig_time))
    if max_MF>=Mig_Timing_Variation: 
        buf_size = max_MF//(Mig_Timing_Variation)
        if max_MF%buf_size!=0: buf_size+=1
    else: raise Exception("max_MF<Mig_Timing_Variation")    
    mig_timing_ind=0
    for i in range(Num_parents):
        cell_position=random_sample_2d(area_length,area_width,MF_sample_position[i])
        MFs.append(Cells('MFR', init_position=cell_position, \
            color_=colormap_for_mig_time[mig_timing_ind], vpython=vpython))        
        if len(MFs)%buf_size==0: 
            #print('len(MFs)',len(MFs), 'mig_timing_ind', mig_timing_ind)
            if mig_timing_ind<Mig_Timing_Variation:mig_timing_ind+=1

        positions_childs=positions_MFs_childs(cell_position, Num_childs)
        for j in positions_childs:
            MFs.append(Cells('MFR', init_position=j, \
                color_=colormap_for_mig_time[mig_timing_ind], vpython=vpython))
            if len(MFs)%buf_size==0: 
                #print('len(MFs)',len(MFs), 'mig_timing_ind', mig_timing_ind)
                if mig_timing_ind<Mig_Timing_Variation: mig_timing_ind+=1

    return MFs, colormap_for_mig_time

def num_cells_per_MF_clusters(MFs):
    early, mid, late =0, 0, 0
    for i in MFs:
        if i.rosette.color==MF_colors[0]:
            early+=1
        elif i.rosette.color==MF_colors[1]:
            mid+=1
        elif i.rosette.color==MF_colors[2]:
            late+=1
    print('# MFs per migration timing:','early:', early, 'mid',mid,'late',late)


"""
def GC_coloring(len_GCs, Max_GC, time_division=1):
    #if len_GCs< int(Max_GC/3):
    if len_GCs< Max_GC/3:
        return color_list[0]  #red
    elif len_GCs < Max_GC*2/3:
        return color_list[-1]  #black
    else:
        return color_list[3] #green

def GC_Single_coloring(len_GCs, Max_GC, time_division=1):
    if len_GCs< Max_GC*0.3:
        return color_list[0]
    elif len_GCs> Max_GC*0.7:
        return color_list[3]
    else:
        return color_list[-1]"""


class non_v_cell:
    def __init__(self, pos, radius, color):
        self.pos=pos
        self.radius=radius
        self.color=color

class Cells: 
    def __init__(self, cell_type, init_position, color_=color.white, vpython=False):
        #self.height_PCL=height_PCL
        self.init_position=init_position
        self.color=color_
        self.flag_arrived_final_destination=False
        self.synapse_partners=[]
        self.vpytyhon=vpython
        if cell_type=='MFR':
            self.cell_type='MFR'
            self.radius=radius_MFR
            if vpython==True:
                self.rosette=sphere(pos=init_position, radius=self.radius, color=self.color)
            else:
                self.rosette=non_v_cell(pos=init_position, radius=self.radius, color=self.color)
            self.superposition=False
            self.activity_level=0
            #self.above_PCL=False        
            #self.radiation1=sphere(pos=init_position, radius=radiation1, color=color.magenta, opacity=0.3)
            #self.radiation2=sphere(pos=init_position, radius=radiation2, color=color.magenta, opacity=0.2)
            #self.radiation3=sphere(pos=init_position, radius=radiation3, color=color.magenta, opacity=0.15)
            #self.body=compound([self.rosette, self.radiation3, self.radiation2, self.radiation1])
            self.body=self.rosette
        elif cell_type=='GC':
            self.gauge=0.01
            self.cell_type='GC'
            self.radius=radius_GC
            self.superposition=False            
            if vpython==True:
                self.body=ellipsoid(pos=init_position, length=self.radius, \
                                                    height=self.radius*4, \
                                                    width =self.radius, color=self.color)
            else:
                self.body=non_v_cell(pos=init_position, radius=self.radius, color=self.color)

                                                
        elif cell_type=='GC_sample':
            self.gauge=0
            labelpos=init_position+vector(0,0,10)
            self.display=label(pos=labelpos, text=self.gauge)
            self.cell_type='GC_sample'
            self.radius=radius_GC
            #self.destination=0
            self.superposition=False            
            if vpython==True:
                self.body=ellipsoid(pos=init_position, length=self.radius, \
                                                    height=self.radius*4, \
                                                    width =self.radius, color=self.color)
            else:
                self.body=non_v_cell(pos=init_position, radius=self.radius, color=self.color)
        else: 
            raise Exception('Unknown cell types')
    def getter(self, attribute):
        t= id(attribute)
        return t
        
    def moveup(self):
        self.body.pos.y+= 0.12*height_PCL/(24*7*time_division)
        # 0.15% of MFs arrive at PCL by P0 
        # While another 0.15% are supposed to be arrived at PCL till P7
        # -> The time that proportional population of randomly sampled MFs w.r.t. the upper height of PCL*0.15 take to arrive PCL for 7 days =24*7 hours

    
    def movedown(self, current_depth_IGL=height_PCL):        
        if self.cell_type=='GC':
            if not self.flag_arrived_final_destination:
                if self.body.pos.y>1e-3+current_depth_IGL:  #precision 0.001
                    self.body.pos.y-=height_PCL/(10*time_division)
                    if self.body.pos.y<=1e-3+current_depth_IGL:
                        self.body.pos.y=current_depth_IGL
                        GC_migration_complete(self, self.vpytyhon)
                        self.flag_arrived_final_destination=True
            #if self.flag_arrived_final_destination and not self.body.pos.y<self.radius:
            #    self.body.pos.y-=(height_PCL*0.20)/(24*7)
        elif self.cell_type=='GC_sample': #Only difference is the existence of the gauge-display
            if not self.flag_arrived_final_destination:
                if self.body.pos.y>1e-3+current_depth_IGL:  #precision 0.001
                    self.body.pos.y-=height_PCL/(10*time_division)
                    self.display.pos.y-=height_PCL/(10*time_division)
                    if self.body.pos.y<=1e-3+current_depth_IGL:
                        self.body.pos.y=current_depth_IGL
                        GC_migration_complete(self, self.vpytyhon)
                        self.flag_arrived_final_destination=True
        
        elif self.cell_type=='MFR':            
            self.body.pos.y-=(height_PCL*0.20)/(24*7*time_division)



''' 
#obsolete functions
def GC_destination(current_depth_IGL):
    strata_subdivision=['top', 'mid', 'bottom']
    prob_destination=np.random.choice(strata_subdivision, p=[0.06,0.12,0.82])        
    if prob_destination=='top':
        transit_distance=np.random.choice(int(current_depth_IGL*0.44)+1)
    elif prob_destination=='mid':
        transit_distance=np.random.choice(int(current_depth_IGL*0.44)+1)+current_depth_IGL*0.44
    elif prob_destination=='bottom':
        transit_distance=np.random.choice(int(current_depth_IGL*0.12)+1)+current_depth_IGL*0.88        
    if transit_distance<radius_GC:
        transit_distance=radius_GC
    #print('transit_distance: ', transit_distance)        
    return transit_distance
    
'''