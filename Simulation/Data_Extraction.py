# Data management
import matplotlib.pyplot as plt 
from .Parameters import *
import os
import numpy as np

def data_save(data_name, data, PATH):
    #cur_path = os.path.dirname(os.path.realpath(__file__))
    directory = PATH +'/data'
    if not os.path.exists(directory):
        os.makedirs(directory)

    save_as=directory+'/'+data_name
    np.save(save_as, data)
    print('Data', data_name,'Saved at:',PATH)
    
def data_load(data_name):
    cur_path = os.path.dirname(os.path.realpath(__file__))
    directory = cur_path+'/data'
    #load= np.load(directory+'/'+ data_name+'.npy')
    load= np.load(directory+'/'+ data_name+'.npy', allow_pickle=True)
    return load


import pickle

def pickle(dataname, data):
    pickle.dump(data, open(dataname+".pkl", "w"))

def unpickle(dataname):
    pickle.load(open(dataname+".pkl"))


#def data_extraction_3d(cells, data_name, FLAG_save):
def data_extraction_3d(cells, data_name):
    
    print('extracting data', data_name,'...')
    inspection_range=np.linspace(0+1, area_width-1, area_width+1-2)

    total_inspection=[]
    for ins_plain in inspection_range:

        #inspect=np.array([])
        inspect=[]
        for cell in cells:        
            #distance= cell.body.pos.z - 0.5*area_width    #distance from center
            distance= cell.body.pos.z - ins_plain
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
    if True:
        #data_save(data_name, inspect)
        data_save(data_name, total_inspection)
    print('close the plot to proceed further...')

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
