# Data management
import matplotlib.pyplot as plt 
from .Parameters import *
import os
import numpy as np

def data_save(data_name, data, PATH, dir_name):

    #cur_path = os.path.dirname(os.path.realpath(__file__))
    directory = PATH +'/data/'+dir_name
    if not os.path.exists(directory):
        os.makedirs(directory)

    save_as=directory+'/'+data_name
    np.save(save_as, data)
    print('Data', data_name,'Saved at:', directory)
    
def data_load(data_name, PATH):
    #load= np.load(directory+'/'+ data_name+'.npy')
    load= np.load(PATH+'/'+ data_name, allow_pickle=True)
    return load

def load_all(folder_name, PATH, Print=True):
    Location=PATH+'/data/'+folder_name    
    contents = os.listdir(Location)
    for c in contents:
        if 'GC_Objects' in c:
            GC_Objects = np.ndarray.tolist(data_load(c, Location))
        elif 'MF_Objects' in c:
            MF_Objects = np.ndarray.tolist(data_load(c, Location))
        elif 'GC_Colormap' in c:
            GC_Colormap = np.ndarray.tolist(data_load(c, Location))
        elif 'MF_Colormap' in c:
            MF_Colormap = np.ndarray.tolist(data_load(c, Location))
        elif 'Graph' in c:
            #Graph = np.ndarray.tolist(data_load(c, Location))
            Graph = data_load(c, Location)
    
    if Print:
        print('Data Loaded from:', Location)
        type_and_shape('GC_Objects', GC_Objects)
        type_and_shape('MF_Objects', MF_Objects)
        type_and_shape('GC_Colormap', GC_Colormap)
        type_and_shape('MF_Colormap', MF_Colormap)
        type_and_shape('Graph', Graph)
    return GC_Objects, MF_Objects, GC_Colormap, MF_Colormap, Graph

def load_coefficient(folder_name, PATH, DATA_NAME):
    Location=PATH+'/data/'+folder_name    
    contents = os.listdir(Location)
    for c in contents:
        if c==DATA_NAME:
            data=np.ndarray.tolist(data_load(c, Location))

    #print(np.shape(data), type(data))
    
    t_range =data[0]
    List_Assr_coeff_MF=data[1]
    List_Assr_coeff_GC=data[2]
    
    return t_range, List_Assr_coeff_MF, List_Assr_coeff_GC



def type_and_shape(Name, data):
    print(str(Name), 'type:', type(data), 'data shape:', np.shape(data))


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
