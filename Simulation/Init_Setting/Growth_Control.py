
from ..Parameters import *
import numpy as np
import matplotlib.pyplot as plt 

def exponential_growth(graph=False):
    time_range=simulation_time
    #dt = np.linspace(0, 0.95*time_range, int(0.95*time_range))
    dt = np.linspace(0, time_range, time_range+1)
    #arr=np.zeros(time_range - int(len(dt)))    
    
    
    max_num=height_PCL
    lable_y='IGL Depth'
     
    y= (1-np.exp(-0.013*dt/time_division))
    #print('time range:', time_range, 'sp dt:', np.shape(dt), 'sp arr', np.shape(arr), 'sp y', np.shape(y))
    #y=np.concatenate((arr,y))
    present_value=height_PCL* y #Note that without compensation, present_cell[-1] < max_num    
    if present_value[-1]!=height_PCL:
        present_value[-1]=height_PCL

    if graph:
        t=np.linspace(0, time_range, time_range+1)  
        plt.plot(t, present_value)
        plt.axvline(x=14*24*time_division)
        plt.axvline(x=7*24*time_division)
        #plt.plot(dt, max_num-present_value) 
        plt.xlabel("Time lapsed (h)")
        plt.ylabel(lable_y)
        plt.show()    

    return max_num-present_value

def linear_growth(Max_GC, graph=False): 
    time_range=24*20*time_division
    peak_period = np.linspace(24*7*time_division, 24*14*time_division, 24*7*time_division+1)
    y= Max_GC*(peak_period-24*7*time_division)/(24*7*time_division)
    before_peak=np.zeros(24*7*time_division)
    after_peak=Max_GC*np.ones(24*6*time_division)
    y_concat=np.concatenate([before_peak,y,after_peak])

    #print('before', np.shape(before_peak))
    #print('peak', np.shape(peak_period))
    #print('after', np.shape(after_peak))
    #print('concat', np.shape(y))
    #print('due num', 24*20*time_division)
    
    if graph:
        dt = np.linspace(0, time_range, time_range+1)        
        plt.plot(dt, y_concat)
        plt.axvline(x=14*24*time_division)
        plt.axvline(x=7*24*time_division)
        #plt.plot(peak_period, y) 
        #plt.plot(dt, max_num-present_value) 
        plt.xlabel("Time lapsed (h)")
        plt.ylabel('Linear growth')
        plt.show()    
    return y_concat

def logistic_growth(factor, max_value, graph=False):
    #time_range=24*20 #P0~P19
    time_range=simulation_time
    dt = np.linspace(-time_range/2, (time_range/2), time_range)  
    max_num=max_value
    
    if factor=='IGL':        
        lable_y='IGL Depth'
    if factor=='GC':
        #max_num=3*max_value  #cuz ratio GC:MF is constrained as 3:1 
        lable_y='Num GCs'
    
    z = 1/(1 + np.exp(-0.03*dt/time_division)/time_division)    #logistic function    
    present_value=z*max_num #Note that without compensation, present_cell[-1] < max_num
    #while(present_value[-1]<max_num): #compensation for the loss from exp calculation.
    #    present_value+=1    
    if graph:
        plt.axvline(x=14*24*time_division)
        plt.axvline(x=7*24*time_division)
        plt.plot(dt+time_range/2, present_value) 
        plt.xlabel("Time lapsed (h)")
        plt.ylabel(lable_y)
        plt.show()
    if factor=='IGL':
        present_value=max_value-present_value
    if factor=='GC':
        present_value[-1]=max_num

    
    return present_value
