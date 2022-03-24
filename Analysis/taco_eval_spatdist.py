from vpython import *
import numpy as np
import matplotlib.pyplot as plt 
import time
from K_func import *


#from hanging_threads import start_monitoring
#start_monitoring(seconds_frozen=10, test_interval=100)

DEFEAULT_point_data='empty'
DEFEAULT_SHOW_SLICE=True
DEFEAULT_Eval_K=False



def main(point_data, Show_slice, Eval_K):
    point_data='not empty'
    if point_data=='empty':
        print('No edge data')
    else:
        point_data='tolerance0.1_160'

        normal_cells_data=point_data+'_normal_MFs'
        normal_cells = data_load(normal_cells_data)
        normal_cells[:,2]=np.round(normal_cells[:,2])

        injected_cells_data=point_data+'_injected_MFs'
        injected_cells = data_load(injected_cells_data)
        injected_cells[:,2]=np.round(injected_cells[:,2])

        print('Num total cells', len(injected_cells)+len(normal_cells))
        print('Inected cells', len(injected_cells))
        print('Normal_cells', len(normal_cells))


        """
        Show 2d slicing through vpython----------------------------------------------------
        """ 
        if Show_slice:
            radius_scale=1
            cr = shapes.circle(radius=1*radius_scale, thickness=0.1)
            cr_hole = shapes.circle(radius=0.9*radius_scale)
            cr2 = shapes.circle(radius=1*radius_scale)

            points=[]
            for i in normal_cells:
                ext1=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr_hole, color=color.black)
                ext2=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr, color=color.red)
                ext=compound([ext1, ext2])
                points.append(ext)
            print('Red cell done')
            for i in injected_cells:
                ext=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr2, color=color.green)
                points.append(ext)
            print('Green cell done')
                

        """
        Evalution of spatial distribution of cells through K function----------------------------------------------------
        """ 
        if Eval_K:
            total_data=np.concatenate((normal_cells_data, injected_cells_data))
            total_data=total_data[:,:2]
            list_del=[]
            for i in enumerate(total_data):
                if i[1][0]>area_length or i[1][1]>height_PCL:
                    list_del.append(i[0])
            squared_data=np.delete(point_data, list_del,0)
            plt.scatter(point_data[:,0], point_data[:,1])
            #plt.scatter(point_data[:,0], point_data[:,1])
            #plt.scatter(squared_data[:,0], squared_data[:,1])

            plt.title('Data remained outlier removal')
            plt.show()

            my_K_func(squared_data)
    
        

import argparse
if __name__ == '__main__':

    start_time = time.time()
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--point_data', type=str, dest='point_data',
                        #default=DEFEAULT_point_data,
                        help='Data of cell points that is obtained from a slice of IGL volume')
    
    parser.add_argument('--draw_dist', type=int, dest='show_slice',
                        default=DEFEAULT_SHOW_SLICE,
                        help='Show visualization of 2D slice of cells')

    parser.add_argument('--draw_net', type=int, dest='Eval_K',
                        default=DEFEAULT_Eval_K,
                        help='Show spatial distribution analysis through K function.')

    args = parser.parse_args()
    
    main(args.point_data, args.show_slice, args.Eval_K)

    elapsed_time = time.time() - start_time
    print('Total elapsed time for evaluation:', elapsed_time)