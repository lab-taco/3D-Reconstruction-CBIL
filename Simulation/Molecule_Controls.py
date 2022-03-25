#from Parameters import *

from vpython import *
import numpy as np

def exponentiated_exponential_function(x, x_var=0, time_division=1, lambda_=3.5, alpha_=5.5, sigma=100, k=100, graph=False):
    #with sigma variation    
    # Function range = 0.5*sigma
    if x<=x_var*time_division:
        return 0
    eef=k*alpha_*lambda_\
        *np.exp(-lambda_*(x-x_var*time_division)/(sigma*time_division))\
        *((1-np.exp(-lambda_*(x-x_var*time_division)/(sigma*time_division)))**(alpha_-1))    
    #if graph:
    #    plt.plot(x, eef)
    #    plt.show()
    return eef

def MF_activity_coordinator(time_step, time_division):
    t=200
    interval=50
    Early_MF_activity=exponentiated_exponential_function(time_step, \
                x_var=t-interval, time_division=time_division)+1
    Mid_MF_activity=exponentiated_exponential_function(time_step, \
                x_var=t, time_division=time_division )+1
    Late_MF_activity=exponentiated_exponential_function(time_step, \
                x_var=t+interval, time_division=time_division)+1
    MF_activities=[Early_MF_activity, Mid_MF_activity, Late_MF_activity]
    return MF_activities


def MF_activity_coordinator_separatedactivities(time_step, time_division):
    t=220
    k= 140
    s=60
    Early_MF_activity=exponentiated_exponential_function(time_step, \
                x_var=t-100, time_division=time_division, sigma=s, k=k)+1
    Mid_MF_activity=exponentiated_exponential_function(time_step, \
                x_var=t, time_division=time_division, sigma=s/3, k=k)+1
    Late_MF_activity=exponentiated_exponential_function(time_step, \
                x_var=t+50, time_division=time_division, sigma=s, k=k)+1
    MF_activities=[Early_MF_activity, Mid_MF_activity, Late_MF_activity]
    return MF_activities
