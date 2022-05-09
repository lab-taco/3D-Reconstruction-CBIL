
from vpython import *
import numpy as np
import math
import sys
#from function import*
from .Parameters import *

def distance(S1, S2):
    distance=math.sqrt((S1.body.pos.x-S2.body.pos.x)**2\
                      +(S1.body.pos.y-S2.body.pos.y)**2\
                      +(S1.body.pos.z-S2.body.pos.z)**2)
    return distance

def normalization(valueset): #Luce's Choice Axiom
    return np.array(valueset)/np.array(valueset).sum()

def softmax(valueset): # Exponential
    return np.exp(valueset)/np.exp(valueset).sum()
    

# Distance/Activity dependent
def Form_Synapse(GCs, MFs):
    for GC in [cell for cell in GCs if cell.flag_arrived_final_destination\
                                     and len(cell.synapse_partners)==0]:
        distance_values=[] #distribution of distance        
        activity_values=[]
        INDICES_MFs=np.arange(len(MFs))
        for mf in MFs:  
            distance_values.append(distance(GC, mf))
            activity_values.append(mf.activity_level)            
            #sorted_distance_values=distance_values.sort(reverse=True)
        Contact_rate_by_distance = softmax(distance_values)
        #distance_values=np.array(distance_values)
        
        #print('distance_values',distance_values)
        #print('np.log(softmax(distance_values)', np.log(softmax(distance_values)))

        IND_initial_contacts=np.random.choice(INDICES_MFs, 10, \
                        replace=False, p=Contact_rate_by_distance) #draw 10 mfs by distance
        #print('Contact_rate_by_distance', Contact_rate_by_distance)        
        #print('IND_initial_contacts',IND_initial_contacts)
        #print('activity_values', activity_values)
        
        activity_values=np.array(activity_values)
        activity_of_the_contacted = np.take(activity_values, IND_initial_contacts)
        #print('activity_of_the_contacted', activity_of_the_contacted)
        #print('normalization_contacted', normalization(activity_of_the_contacted))

        IND_activity_draws=np.random.choice(IND_initial_contacts, 4, \
                        replace=False, p=softmax(activity_of_the_contacted)) #draw 4 mfs by activity        
        synapsed_MFs=[MFs[ind] for ind in IND_activity_draws]
        #print('IND_synapsed_MFs', IND_activity_draws)
        #print(synapsed_MFs)

        #synapsed_MFs=[MFs(np.where(distance_values==distance)) for distance in draws]
        #Selection of MFs that form actual synapse
        
        for mf in synapsed_MFs:
            GC.synapse_partners.append(mf)
            mf.synapse_partners.append(GC)

def Initial_Contacts(GCs, MFs):
    for GC in [cell for cell in GCs if cell.flag_arrived_final_destination\
                                     and len(cell.synapse_partners)==0]:
        distance_values=[] #distribution of distance        
        activity_values=[]
        INDICES_MFs=np.arange(len(MFs))
        for mf in MFs:  
            distance_values.append(distance(GC, mf))
            activity_values.append(mf.activity_level)            
            #sorted_distance_values=distance_values.sort(reverse=True)
        Contact_rate_by_distance = softmax(distance_values)

        IND_initial_contacts=np.random.choice(INDICES_MFs, 10, \
                        replace=False, p=Contact_rate_by_distance) #draw 10 mfs by distance    
        synapsed_MFs=[MFs[ind] for ind in IND_initial_contacts]
        
        for mf in synapsed_MFs:
            GC.synapse_partners.append(mf)
            mf.synapse_partners.append(GC)
""" 
# Distance/Activity dependent
def Form_Synapse(GCs, MFs, MF_activities):
    for GC in [cell for cell in GCs if cell.flag_arrived_final_destination\
                                     and len(cell.synapse_partners)==0]:
        distance_values=[] #distribution of distance
        relative_activity= []
        for mf in MFs:  
            distance_values.append(distance(GC, mf))
            relative_activity.append(mf.activity_level/MF_activities)
            #print('mf.activity_level value', mf.activity_level)
            #relative_activity.append(softmax(np.array(relative_activity)))
            #sorted_distance_values=distance_values.sort(reverse=True)
        distance_values=np.array(distance_values)

        #for i in MF_activities:
        #    print('MF activity', i)
        #print('Distances:', distance_values)
        #print('Distances softmax:', np.around(softmax(distance_values), decimals=3))
        #print('Relative activity:', relative_activity)     

        val = relative_activity/distance_values
        #print('val \n', val)

        #print('val/np.sum(val) \n', val/np.sum(val))
        #print('Prob', val/np.sum(val))
        #print('prob', negative_softmax(distance_values/relative_activity))
        #print('prob', softmax(relative_activity/distance_values))
        draws=np.random.choice(distance_values, 4, \
                        replace=False, p=val/np.sum(val)) #draw 10 mfs by distance                        
                        #replace=False, p=negative_softmax(distance_values)) #draw 10 mfs by distance

        #print('draws',draws)
        distance_values=np.ndarray.tolist(distance_values)
        synapsed_MFs=[]
        for distance_val in draws:
            ind=distance_values.index(distance_val)
            synapsed_MFs.append(MFs[ind])
        #print(synapsed_MFs)

        #synapsed_MFs=[MFs(np.where(distance_values==distance)) for distance in draws]
        #Selection of MFs that form actual synapse
        
        for mf in synapsed_MFs:
            GC.synapse_partners.append(mf)
            mf.synapse_partners.append(GC)              
"""

def num_synapse(GCs,MFs):
    num_synapse_GC=0
    num_synapse_MF=0
    No_synapse_GC=0
    for ind, gc in enumerate(GCs):
        num_synapse_GC+=len(gc.synapse_partners)
        if len(gc.synapse_partners)<4:
            print(ind,'th GC len(gc.synapse_partners)',len(gc.synapse_partners))
    for mf in MFs:
        #print(MFs.index(mf),'th',len(mf.synapse_partners))
        num_synapse_MF+=len(mf.synapse_partners)

    print('# synapse GC:',num_synapse_GC, '# GC:', len(GCs), '#synapse/GC:',num_synapse_GC/len(GCs))
    print('# synapse MF:',num_synapse_MF, '# MF:', len(MFs), '#synapse/MF:',num_synapse_MF/len(MFs))

def len_synapse(GCs, MFs):

    sum_avg=0 #AVG. length of dendrites for individual GCs
    k=0
    for gc in GCs:
        k+=1
        sum_length=0
        if len(gc.synapse_partners)>0:
            for mf in gc.synapse_partners:
                sum_length+=distance(gc, mf) #each dendrite's length            
            sum_avg+=sum_length/len(gc.synapse_partners) #avg_length per gc
        else:
            print(k,'-th GC has no synaptic partner',\
                    'migration finished?:',gc.flag_arrived_final_destination,\
                    'collision?:', gc.superposition)
        
    avg_len_synapse= sum_avg/len(GCs)   #total AVG. length of GC dendrites
    print('Average length of GC dendrites:',avg_len_synapse)
        

