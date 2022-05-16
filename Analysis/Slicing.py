from vpython import *
import numpy as np
from matplotlib import pyplot as plt
"""
        Show 2d slicing through vpython----------------------------------------------------
""" 

def Capture_Cells_at_Slice(MF_Objects, GC_Objects, Map_size_3D, Cell_radii, \
                        PLOTTING=True, Where_to_Slice='Random'):
    len_x_axis, len_y_axis, len_z_axis= Map_size_3D
    if Where_to_Slice=='Random': 
        Slicing_Point=np.random.choice(len_y_axis)
    radius_MFR, radius_GC = Cell_radii
    Visible_range=radius_GC*3

    MF_Captured=[]
    GC_Captured=[]
    All_cells=MF_Objects+GC_Objects
    for cell in All_cells:
        if cell.body.pos.x+cell.radius<=Slicing_Point+Visible_range\
        and cell.body.pos.x-cell.radius>=Slicing_Point-Visible_range:
            if cell in MF_Objects: MF_Captured.append(cell)
            elif cell in GC_Objects: GC_Captured.append(cell)

    if PLOTTING:
        if len(MF_Captured+GC_Captured)>0:
            Position_MFs=[[mf.body.pos.x, mf.body.pos.y] for mf in MF_Captured]
            Position_GCs=[[gc.body.pos.x, gc.body.pos.y] for gc in GC_Captured]
            
            plt.scatter(np.array(Position_MFs)[:, 0], np.array(Position_MFs)[:, 1], label='MFs')
            plt.scatter(np.array(Position_GCs)[:, 0], np.array(Position_GCs)[:, 1], label='GCs')
            plt.legend(title='Captured Cells')
            plt.title('Visualized Slice in 2D')
            plt.show()
        else:
            print('There is no cells in the slice, \
                try adjust visible points or increase cell populaion size')

    return MF_Captured, GC_Captured

def Show_slice(normal_cells, injected_cells):
    radius_scale=1
    cr = shapes.circle(radius=1*radius_scale, thickness=0.1)
    cr_hole = shapes.circle(radius=0.9*radius_scale)
    cr2 = shapes.circle(radius=1*radius_scale)

    points=[]
    for i.body.pos in normal_cells:
        ext1=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr_hole, color=color.black)
        ext2=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr, color=color.red)
        ext=compound([ext1, ext2])
        points.append(ext)
    print('Red cell done')
    for i in injected_cells:
        ext=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr2, color=color.green)
        points.append(ext)
    print('Green cell done')