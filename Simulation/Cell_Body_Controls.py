from .Parameters import *
from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt 
import sys
import time
sys.setrecursionlimit(10**6)  #function recursion limit, default was 1000


"""
def config():
    height_PCL=radius_GC*42
    area_length=height_PCL
    area_width=height_PCL

    radiation1=radius_MFR*2.0*scale_radiation_range
    radiation2=radius_MFR*5.0*scale_radiation_range
    radiation3=radius_MFR*8.0*scale_radiation_range

    intensity_radiation1=0.01/time_division
    intensity_radiation2=0.005/time_division
    intensity_radiation3=0.001/time_division
    return height_PCL, area_length, area_width\
            ,radiation1, radiation2, radiation3\
            ,intensity_radiation1, intensity_radiation2, intensity_radiation3"""

#area_length=300
#area_width=300
#height_PCL=300

#area_length=85
#area_width=17
#height_PCL=113

#area_length=10
#area_width=10
#height_PCL=50


def after_migration_coloring(subgroup, allcell):
    for cell in allcell:
        if cell in subgroup:
            cell.body.color=color.green
        else:
            cell.body.color=color.red
     
'''

def flood_check(collided_GC_position, ricochet_vector):    
    floodcheck_y = collided_GC_position.body.pos.y+ricochet_vector.y
    floodcheck_x = collided_GC_position.body.pos.x+ricochet_vector.x
    floodcheck_z = collided_GC_position.body.pos.z+ricochet_vector.z
    if floodcheck_y > height_PCL-radius_GC or floodcheck_y<0+radius_GC:
        x_flood=True
    if floodcheck_x > area_length-radius_GC or floodcheck_x<0+radius_GC:
        y_flood=True
    if floodcheck_z > area_width-radius_GC or floodcheck_z<0+radius_GC:
        z_flood=True
    else: flood=None
    return flood
'''

def collision_check(S1, S2):
    distance=math.sqrt((S1.body.pos.x-S2.body.pos.x)**2\
                      +(S1.body.pos.y-S2.body.pos.y)**2\
                      +(S1.body.pos.z-S2.body.pos.z)**2)
    if distance+1e-4<(S1.radius+S2.radius)*(1-tolerance):
        S1.superposition, S2.superposition= True, True
        return True, distance        
    else: 
        return False, None

def ricochet(cell, collided_cell):
    distance=math.sqrt((cell.body.pos.x-collided_cell.body.pos.x)**2\
                    +(cell.body.pos.y-collided_cell.body.pos.y)**2\
                    +(cell.body.pos.z-collided_cell.body.pos.z)**2)
    necessary_distance=(cell.radius+collided_cell.radius)*(1-tolerance)
    if distance>necessary_distance:
        print('collided pos:',collided_cell.body.pos,'cell pos:',cell.body.pos,\
               'distance:',distance)
        print('cell type:',cell.cell_type, 'collided_cell type:', collided_cell.cell_type,\
                'necessary distance:', necessary_distance)
        #    'norm_v:', norm_v,'distance:',distance)
        raise Exception('No need ricohet: enough distance')
    elif distance==0:
        ricochet_vector=vector(necessary_distance,necessary_distance,necessary_distance)
    #collided_cell.body.color=color.cyan        
    else:
        v=vector((collided_cell.body.pos.x-cell.body.pos.x)\
            ,(collided_cell.body.pos.y-cell.body.pos.y)\
            ,(collided_cell.body.pos.z-cell.body.pos.z))
        #norm_v=distance=math.sqrt( (v.x)**2+(v.y)**2+(v.z)**2 )        
        new_norm=necessary_distance-distance
        ricochet_vector=new_norm*(v/distance)
        #if too close to layer division
    '''
    flood_check=collided_cell.body.pos.y+ricochet_vector.y
    if flood_check > height_PCL-radius_cell or flood_check < 0+radius_cell:        
        xz_direction_vector=vector(collided_cell.body.pos.x-cell.body.pos.x,0,collided_cell.body.pos.z-cell.body.pos.z)\
                        /sqrt((collided_cell.body.pos.x-cell.body.pos.x)**2+(collided_cell.body.pos.z-cell.body.pos.z)**2)
        if flood_check> height_PCL-radius_cell: 
            y_vector=vector(0, height_PCL-collided_cell.body.pos.y-radius_cell, 0)
        elif flood_check< 0+radius_cell: 
            y_vector=vector(0, 0-collided_cell.body.pos.y+radius_cell, 0)
        ricochet_vector=xz_direction_vector+y_vector
        #collided_cell.body.color=color.magenta
    '''
    collided_cell.body.pos+=ricochet_vector
    #collided_cell.final_position=collided_cell.body.pos
    


def collision_control(new_cell, collided_cells, cells):
    #print('length collided cells:', len(collided_cells))
    for col_cell in [collided_cell for collided_cell in collided_cells]:
        ricochet(new_cell, col_cell) #ricochet off each collided cell from the target cell
    new_cell.superposition=False

    for col_cell in [collided_cell for collided_cell in collided_cells]:        
        superposition_check(col_cell, cells) # same procedure for each ricocheted cells        


def superposition_check(new_cell, cells):
    #if new_cell.cell_type=='GC' and not new_cell.body.__class__.__name__=='sphere':
    #    new_cell.body.color=color.cyan
    #    for i in range(10):
    #        sleep(1)
    #        new_cell.body.pos.z+=3
    #    raise Exception('incomplete migration')        
    
    collided_cells=[] 
    #collision detection loop
    for target_cell in \
        [cell for cell in cells if cell!=new_cell and cell.flag_arrived_final_destination]:  
        collision, _ = collision_check(new_cell, target_cell)
        if collision==True:            
            collided_cells.append(target_cell)
    
    
    if len(collided_cells)>0: # If collision exist, collision control 
        collision_control(new_cell, collided_cells, cells)    
    new_cell.superposition=False

def final_collision_check(all_cells):
    print('final collision checking...')
    remnant_count=0
    sp_count=0    
    for cell in all_cells:
        if cell.superposition:  sp_count+=1        
    print('# superposition remain:', sp_count)
    for cell_checking in all_cells:        
        for target_cell in [cell for cell in all_cells if cell!=cell_checking]:
            final_col, distance =collision_check(cell_checking,target_cell)
            if final_col==True:
                remnant_count+=1               
                
                target_cell.body.color=color.cyan
                cell_checking.body.color=color.magenta
                '''
                print('Collision still remain btwn',all_cells.index(cell_checking),\
                    'and',all_cells.index(target_cell),'distance:',distance)'''
                
                final_col=False
    print('# collision remain', remnant_count)  


def check_out_of_volume(position):
    if position.x>area_length or position.x<0 \
        or position.y>height_PCL or position.y<0 \
        or position.z>area_width or position.z<0:
        return True
    else:
        return False



def flooded(cell):
    cell.body.visible=False
    cell.body=sphere(pos=cell.body.pos, radius=cell.radius, color=cell.color)


def all_superposition_check(cells):       
    for cell in cells:
        if cell.flag_arrived_final_destination:
            superposition_check(cell, cells)

    

def Molecule_absorption(GCs, MFs, collision_check_cells, vpython):
    for GC in [cell for cell in GCs if not cell.flag_arrived_final_destination and cell.gauge<1]:
        for MF in MFs:
            distance=math.sqrt((GC.body.pos.x-MF.rosette.pos.x)**2\
                              +(GC.body.pos.y-MF.rosette.pos.y)**2\
                              +(GC.body.pos.z-MF.rosette.pos.z)**2)
            if distance<radius_GC+radiation1:
                #GC.gauge+=intensity_radiation1
                GC.gauge+=MF.activity_level*intensity_radiation1
            elif distance<radius_GC+radiation2:
                #GC.gauge+=intensity_radiation2-
                GC.gauge+=MF.activity_level*intensity_radiation2
            elif distance<radius_GC+radiation3:
                #GC.gauge+=intensity_radiation3
                GC.gauge+=MF.activity_level*intensity_radiation3
            
            if GC.gauge>=1:
                GC.gauge=1
        if GC.cell_type=='GC_sample':
            GC.display.text=GC.gauge
        #print('Gauge:',GC.gauge)        

def Stop_or_not(GCs, collision_check_cells, vpython):
    for GC in [cell for cell in GCs if not cell.flag_arrived_final_destination]:
        s_o_n = ['stop', 'not']
        draw=np.random.choice(s_o_n, p=[GC.gauge, 1-GC.gauge])
        if draw=='stop':
            GC.flag_arrived_final_destination=True
            GC_migration_complete(GC, vpython)
            superposition_check(GC, collision_check_cells)

def GC_migration_complete(cell, vpython):
    if vpython==True:
        cell.body.visible=False
        cell.body=sphere(pos=cell.body.pos, radius=cell.radius, color=cell.color)
        
#def GC_migrates(cells, current_depth_IGL, time_steps, collision_check_duration):       
def Cell_migratedown(Cells, current_depth_IGL, time_steps):
    for cell in [c for c in Cells if not c.flag_arrived_final_destination]:
        cell.gauge=GC_Migration_Speed*(1/(cell.body.pos.y)) # think about the function curve and gc location >> Currently, migration end does not depend on activity, but only dependant on the location of migrating GC
        #cell.gauge=1
        cell.movedown(current_depth_IGL)               
        #if len(Cells)>1 and gc.flag_arrived_final_destination: \
        if cell.flag_arrived_final_destination: \
            #and time_steps%collision_check_duration==0:
            superposition_check(cell, Cells)


def MF_growup(MFs):    
    for cell in MFs:                       
        cell.moveup()

def MFs_migration_complete(MFs):
    for mf in MFs:
        mf.flag_arrived_final_destination=True
        
def MFs_repulsed(MFs):       
    for cell in MFs:               
        cell.movedown()