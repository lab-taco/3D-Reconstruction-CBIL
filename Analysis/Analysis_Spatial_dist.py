from cProfile import label
import numpy as np
from matplotlib import pyplot as plt


def data_extraction_3d(cells, Map_size_2D, radius_GC, PLOTTING=False):
    x_axis, y_axis = Map_size_2D
    print('extracting data...')
    inspection_range=np.linspace(0+1, x_axis-1, x_axis+1-2)

    total_inspection=[]
    for ins_plain in inspection_range:
        inspect=[]
        for cell in cells:
            distance= cell.body.pos.z - ins_plain #distance from center
            if abs(distance)< radius_GC:
            #distance=cell.body.pos.z
            #if cell.body.pos.z>0.5*x_axis-radius_GC:
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
def data_extraction_2d(cells, Map_size_2D, radius_GC, PLOTTING=False):
    x_axis, y_axis = Map_size_2D
    inspect=[]
    for cell in cells:        
        #distance= cell.body.pos.z - 0.5*x_axis    #distance from center
        #if abs(distance)< radius_GC:
        if cell.body.pos.z>0.5*x_axis-radius_GC:
            inspect.append([round(cell.body.pos.x, 8)\
                            , round(cell.body.pos.y,8)])
    inspect=np.asarray(inspect)

    if PLOTTING:
        plt.scatter(inspect[:,0], inspect[:,1])
        plt.show()

import math
def Euclidean_dist(x2, y2, x1, y1):
    return math.sqrt((x2-x1)**2+(y2-y1)**2)
def sample_xy_coor_random(domain_x, domain_y):
    x = np.random.random_integers(0, domain_x)
    y = np.random.random_integers(0, domain_y)
    return x, y
def random_scattering(num_cells, Map_size_2D, plotting=False, \
                    exclusion_zone=False, append_arr=[]):
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
        plt.scatter(arr[:,0], arr[:,1], label='Random dist., Num Cells:'+str(num_cells))
        #plt.xlim(-1,42)
        #plt.ylim(-1,42)
        #plt.title('Random Spatial Sampling, Num Cells:{}'.format(str(num_cells)))
        ##plt.legend(title='Num Cells:{}'.format(str(num_cells)))
        #plt.show()   
    
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
        plt.scatter(arr[:,0], arr[:,1], label='Regular dist., , Num Cells:'+str(len(arr)))
        #plt.xlim(-1,x_axis)
        #plt.ylim(-1,y_axis)
        #plt.title('Regular Spatial Sampling, Num Grid:{}'.format(str(num_grid)))
        ##plt.legend(title='Num Grid:{}'.format(str(num_grid)))
        #plt.show()
    return arr

def cells_in_scope(origin, data, distance):
    data_in_scope = [cell for cell in data if not np.all(cell==origin)\
                        and cell[0]<origin[0]+distance and cell[0]>origin[0]-distance\
                        and cell[1]<origin[1]+distance and cell[1]<origin[1]-distance]
    return data_in_scope

#def nearest_neighbor_dist(data, range_inspect, plotting_type='hist', plot_label=''):
def nearest_neighbor_distance_distribution(data, range_inspect):
    nearest_distance_distribution=[]
    for origin_cell in data:
        for range_of_scope in range(1, range_inspect):
            inspect = cells_in_scope(origin_cell, data, range_of_scope)
            if len(inspect)>0:
                #print(len(inspect), 'cells are caugt at range', range_of_scope)
                distance_distribution=[]
                for cell in inspect:
                    distance=Euclidean_dist(cell[0], cell[1], origin_cell[0], origin_cell[1])                    
                    distance_distribution.append(distance)
                the_nearest_dist = np.amin(distance_distribution)
                nearest_distance_distribution.append(the_nearest_dist)
                #break
    return nearest_distance_distribution, range_inspect
    #mean = np.round(np.mean(nearest_distance_distribution),3)
    #std  = np.round(np.std(nearest_distance_distribution), 3)
    #from statsmodels.distributions.empirical_distribution import ECDF
    #ecdf = ECDF(nearest_distance_distribution)
    
    #Label= plot_label+' \u03BC{} \u03C3:{}'.format(str(mean) ,str(std))
    #if plotting_type=='hist':
    #    bin_range=np.arange(0, range_inspect)
    #    plt.hist(nearest_distance_distribution, bins=bin_range, \
    #        density=False, align='mid', label=Label)
    #elif plotting_type=='cdf':
    #    from statsmodels.distributions.empirical_distribution import ECDF
    #    ecdf = ECDF(nearest_distance_distribution)
    #    plt.plot(ecdf.x, ecdf.y, label=Label)
    #    
        
    


def Plot_NND(NND, Range_inspect, plotting_type='Histo', plot_label='NND'):
    mean = np.round(np.mean(NND),3)
    std  = np.round(np.std (NND),3)
    Label= plot_label+' \u03BC{} \u03C3:{}'.format(str(mean) ,str(std))
    
    if plotting_type=='Histo':
        bin_range=np.arange(0, Range_inspect)
        plt.hist(NND, bins=bin_range, \
            density=False, align='mid', label=Label)
    elif plotting_type=='CDF':
        from statsmodels.distributions.empirical_distribution import ECDF
        ecdf = ECDF(NND)
        plt.plot(ecdf.x, ecdf.y, label=Label)

def Plot_NND_all_together(NND_List):
    Graph_Type= ['CDF', 'Histo']    
    for graph_type in Graph_Type:
        for nnd_list in NND_List:
            nnd, range_NND, label =nnd_list
            print('Given nnd:', type(nnd), np.shape(nnd), 'label:', label)
            Plot_NND(nnd, range_NND, plotting_type=graph_type, plot_label=label)
        plt.title('Nearest Neighbor Distance, {}'.format(graph_type))
        plt.legend()                                          
        plt.show()


def mean_std_rounded(distribution):
    mean = np.round(np.mean(distribution),3)
    std  = np.round(np.std(distribution), 3)
    return mean, std