import numpy as np
from matplotlib import pyplot as plt



def data_extraction_3d(cells):
    
    print('extracting data...')
    inspection_range=np.linspace(0+1, area_width-1, area_width+1-2)

    total_inspection=[]
    for ins_plain in inspection_range:
        inspect=[]
        for cell in cells:
            distance= cell.body.pos.z - ins_plain #distance from center
            if abs(distance)< radius_GC:
            #distance=cell.body.pos.z
            #if cell.body.pos.z>0.5*area_width-radius_GC:
                cell_coordinate=[ round(cell.body.pos.x, 3)\
                                 ,round(cell.body.pos.y, 3)\
                                 ,round(distance, 3) ]
                #print('Cell coordinate', cell_coordinate)
                inspect.append(cell_coordinate)
                #inspect=np.append(inspect, cell_coordinate, axis=0)
        inspect=np.asarray(inspect)
        total_inspection.append(inspect)

    #for ind, val in enumerate(total_inspection):
    #    print(ind,'-th')
    #    print(val)


    #if FLAG_save:
    #if True:
    #    #data_save(data_name, inspect)
    #    data_save(data_name, total_inspection)
    #print('close the plot to proceed further...')

    '''
    for ind, val in enumerate(total_inspection):
        if len(val)>0:    
            plt.scatter(val[:,0], val[:,1])
        else:
            print('There is no cells for',ind,'-th', data_name)
        plt.show()'''

def data_extraction_2d(cells, data_name, FLAG_save):
    inspect=[]
    for cell in cells:        
        #distance= cell.body.pos.z - 0.5*area_width    #distance from center
        #if abs(distance)< radius_GC:
        if cell.body.pos.z>0.5*area_width-radius_GC:
            inspect.append([round(cell.body.pos.x, 8)\
                            , round(cell.body.pos.y,8)])
    inspect=np.asarray(inspect)

    if FLAG_save:
        data_save(data_name, inspect)
    
    plt.scatter(inspect[:,0], inspect[:,1])
    plt.show()

import math
def Euclidean_dist(x2, y2, x1, y1):
    return math.sqrt((x2-x1)**2+(y2-y1)**2)
    

def sample_xy_coor_random(domain_x, domain_y):
    x = np.random.random_integers(0, domain_x)
    y = np.random.random_integers(0, domain_y)
    return x, y
def random_scattering(num_cells, Map_size_2D, plotting=False, exclusion_zone=False, append_arr=[]):
    N=num_cells
    x_axis, y_axis = Map_size_2D
    arr=[]
    if len(append_arr)>0:
        arr=np.ndarray.tolist(append_arr)

    for _ in range(0, N):
        x, y = sample_xy_coor_random(x_axis, y_axis)
        if exclusion_zone and len(arr)>0:
            exclusion_range=2.5
            excls_count=0
            max_loop=100
            while(excls_count != len(arr)):                                
                for cell in arr:
                    distance=Euclidean_dist(cell[0], cell[1], x, y)                    
                    if distance < exclusion_range:                         
                        break
                    else:
                        excls_count+=1
                if excls_count != len(arr):
                    x, y = sample_xy_coor_random(x_axis, y_axis)
                    excls_count=0
                max_loop+=1
                if max_loop>10000:
                    raise Exception('Too many cells for regular distribution')
            #print('excls count', excls_count, 'num points:', len(arr))
                            
        arr.append([x, y])
        if False:        
            arr2=np.asarray(arr)
            plt.scatter(arr2[:,0], arr2[:,1])
            plt.xlim(-1,43)
            plt.ylim(-1,43)
            plt.show()
    arr=np.asarray(arr)
    if plotting:        
        plt.scatter(arr[:,0], arr[:,1])
        plt.xlim(-1,42)
        plt.ylim(-1,42)
        plt.title('Random Spatial Sampling')
        plt.legend(title='Num Cells:{}'.format(str(num_cells)))
        plt.show()   
    
    return arr

def regular_scattering(num_grid, Map_size_2D, plotting=False):    
    arr=[]
    x_axis, y_axis = Map_size_2D
    for i in range(0, num_grid):
        for j in range(0, num_grid):
            arr.append([x_axis*i/num_grid, y_axis*j/num_grid])
    
    if plotting:
        print(np.shape(arr))
        arr=np.asarray(arr)
        plt.scatter(arr[:,0], arr[:,1])
        plt.xlim(-1,x_axis)
        plt.ylim(-1,y_axis)
        plt.title('Regular Spatial Sampling')
        plt.legend(title='Num Grid:{}'.format(str(num_grid)))
        plt.show()
    return arr