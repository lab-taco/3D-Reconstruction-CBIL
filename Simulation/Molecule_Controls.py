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
    len_curve=curv_width*t_division+superpose*5
    curve_at_time=peak_val*np.sin((np.pi*(np.arange(len_curve))/len_curve))
    val_func[t_shift:t_shift+len_curve]+=curve_at_time
    return val_func

import scipy.stats as stats
import math

def gaussian_curves_activity(x, duration, mu_init, t_shift, Superposition=0, time_division=1):
    mu = mu_init+t_shift
    curve_width= duration/6 # 6 sigma; width of bell curve = peak time duration of each curve
    #variance = 1, sigma=math.sqrt(variance)=1
    sigma =1*time_division*curve_width# do higher than 0
    #sigma2=1*time_division*curve_width*(1-Superposition/100)+1*curve_width*(Superposition/100)*duration # width term
    sigma2=1*time_division*curve_width*(1-Superposition/200)+1*curve_width*(Superposition/200)*duration # width term
    if sigma2<=0: raise Exception("sigma2 value is to low")

    term1 = 1/(sigma*math.sqrt(2*math.pi))
    term2 = np.exp(-0.5*((x-mu)/sigma2)**2)    
    g_curve=term1*term2
    return g_curve



def curvs_generations(num_curvs, t_division, PK=10, SP=100, draw=False):
    #num_curvs: Mig_Timing_Variation, SP: Superposing degree, PK= Peak Value, 
    P7=24*7*t_division
    P14=24*14*t_division
    Early_start=10
    
       
    peak_time_duration=P7//(num_curvs-1)
    #cw=peak_time_duration/6 #Curve Width
    #print('len t', len(t), 'len curv', cw, 'P6', P6*t_division, 'P14',P14*t_division)
    t=np.arange(24*20*t_division) 
    curvs=[]
    for i in range(num_curvs):
        #time_shift=((P7-Early_start)+cw*ind_curv*(1-SP/100))*t_division+(P14-P7)/10*SP/100        
        #act_curve = sign_curv_function(len(t), cw, t_division, t_shift=int(time_shift), peak_val=PK, superpose=SP)
        act_curve = gaussian_curves_activity(t, peak_time_duration, P7, i*peak_time_duration, \
                                                Superposition=SP, time_division=t_division)
        #print('curve', curv, 'at', (P6+cw*curv)*t_division)
        curvs.append(act_curve)
        
    if draw:
        import matplotlib.pyplot as plt
        from matplotlib.offsetbox import AnchoredText
        for cv in curvs:
            plt.plot(t, cv)
        plt.axvline(x=P7, color='k', linestyle='--',label='P7')
        plt.axvline(x=P14, color='k', linestyle='--',label='P14')        
        plt.xlabel('time')
        plt.ylabel('activity')
        #plt.text(len(t)*0.75, min(cv)+1, "Superposing=%d"%SP)
        #plt.text(len(t)*0.1, min(cv), "bias=%d"%min(cv))
        plt.title("Superposing=%d"%SP)
        plt.legend()
        plt.show()
    return curvs