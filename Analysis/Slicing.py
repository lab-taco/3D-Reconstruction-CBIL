from vpython import *
import numpy as np
from matplotlib import pyplot as plt
import sys
"""
        Show 2d slicing through vpython----------------------------------------------------
""" 
def Capture_Cells_at_Slice(MF_Objects, GC_Objects, Map_size_3D, Cell_radii, \
                        PLOTTING=True, Where_to_Slice='Random'):
    len_x_axis, len_y_axis, len_z_axis= Map_size_3D #X, Y, Z = area length, height PCL, area width
    radius_MFR, radius_GC = Cell_radii
    Visible_range=radius_MFR*8
    Visible_volume = Visible_range*2*len_x_axis*len_y_axis

    if Where_to_Slice=='Random': # Slice in the center half
        Slicing_Point=np.random.choice(int(len_z_axis/2))+int(len_z_axis/2)

    MF_Captured=[]
    GC_Captured=[]
    All_cells=MF_Objects+GC_Objects
    for cell in All_cells:
        if cell.body.pos.z+cell.radius<=Slicing_Point+Visible_range\
        and cell.body.pos.z-cell.radius>=Slicing_Point-Visible_range:
            if cell in MF_Objects: MF_Captured.append(cell)
            elif cell in GC_Objects: GC_Captured.append(cell)

    
    Num_Captured = len(MF_Captured)+len(GC_Captured)
    print('Density Captured by Slicing:', Num_Captured/Visible_volume)
    if PLOTTING:
        if len(MF_Captured+GC_Captured)>0:
            Position_MFs=[[mf.body.pos.x, mf.body.pos.y] for mf in MF_Captured]
            Position_GCs=[[gc.body.pos.x, gc.body.pos.y] for gc in GC_Captured]
            
            plt.scatter(np.array(Position_MFs)[:, 0], np.array(Position_MFs)[:, 1], label='MFs')
            plt.scatter(np.array(Position_GCs)[:, 0], np.array(Position_GCs)[:, 1], label='GCs')
            plt.legend(title='Captured Cells')
            plt.title('Sliced surface')
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


def View_3D_Dist(MFs, GCs):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    scatter_MFs = ax.scatter([mf.body.pos.x for mf in MFs] \
                            , [mf.body.pos.y for mf in MFs] \
                            , [mf.body.pos.z for mf in MFs], color='b', marker='o', label='MFs')
    scatter_GCs =ax.scatter([gc.body.pos.x for gc in GCs]\
                            , [gc.body.pos.y for gc in GCs]\
                            , [gc.body.pos.z for gc in GCs], color='r', marker='^', label='GCs')
    plt.title('3D Distribution')
    ax.set_xlabel('X', fontsize = 10)
    ax.set_ylabel('Y', fontsize = 10)
    ax.set_zlabel('Z', fontsize = 10)
    ax.legend([scatter_MFs, scatter_GCs], ['MFs', 'GCs'])
    plt.show()