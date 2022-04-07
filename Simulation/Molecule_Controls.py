#from Parameters import *

from time import time
from vpython import *
import numpy as np

def exponentiated_exponential_function(x, x_var=0, time_division=1, lambda_=3.5, alpha_=5.5, sigma=100, k=100, graph=False):
    #sigma = variation, lambda = peak value
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
    k= 140
    s=60
    t=200
    activity_curve=exponentiated_exponential_function(time_step, \
                x_var=t, time_division=time_division, sigma=s )+1
    return activity_curve


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



def sign_curv_function(time_steps, curv_width, t_division, \
                        t_shift=0, peak_val=10, superpose=0, bias=1):
    val_func=bias*np.ones(time_steps)
    len_curve=curv_width*t_division+superpose
    curve_at_time=peak_val*np.sin((np.pi*(np.arange(len_curve))/len_curve))
    val_func[t_shift:t_shift+len_curve]+=curve_at_time
    return val_func

def curvs_generations(num_curvs, t_division, PK=10, SP=100, draw=False):
    #num_curvs: Mig_Timing_Variation, SP: Superposing degree, PK= Peak Value, 
    P7=24*7
    P6=24*6
    P14=24*14
    t=np.arange(24*20*t_division)
    cw=P7//num_curvs
    #print('len t', len(t), 'len curv', cw, 'P6', P6*t_division, 'P14',P14*t_division)
    curvs=[]
    for curv in range(num_curvs):    
        act_curve = sign_curv_function(len(t), cw, t_division,\
                 t_shift=(P7+cw*curv)*t_division, peak_val=PK, superpose=SP)
        #print('curve', curv, 'at', (P6+cw*curv)*t_division)
        curvs.append(act_curve)
        
    if draw:
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText
        for cv in curvs:
            plt.plot(t, cv)
        plt.axvline(x=P7*t_division, color='k', linestyle='--',label='P7')
        plt.axvline(x=P14*t_division, color='k', linestyle='--',label='P14')        
        plt.xlabel('time')
        plt.ylabel('activity')
        plt.text(len(t)*0.75, min(cv)+1, "Superposing=%d"%SP)
        plt.text(len(t)*0.1, min(cv), "bias=%d"%min(cv))
        plt.legend()
        plt.show()
    return curvs