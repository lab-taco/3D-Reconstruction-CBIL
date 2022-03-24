import time
from vpython import *
import numpy as np
"""
Blinking
"""
def blinking(cells, lamda, Time='Fixed'):    
    print('Blinking')
    beta=1/lamda
    while True:
        delta_t = np.random.exponential(beta)
        time.sleep(delta_t)
        #print('dt:', delta_t)
        #firing_cell=np.random.choice(cells, 10)
        firing_cell=np.random.choice(cells)
        firing_cell.rosette.color=color.cyan
        time.sleep(0.05)
        firing_cell.rosette.color=firing_cell.color


def blinking2(cells, lamda, color):        
    beta=1/lamda
    delta_t = np.random.exponential(beta)
    time.sleep(delta_t)
    #print('dt:', delta_t)
    firing_cell=np.random.choice(cells, int(lamda))
    #firing_cell=np.random.choice(cells)
    for i in firing_cell:
        i.rosette.color=color
    time.sleep(0.05)
    for i in firing_cell:
        i.rosette.color=i.color


def Blink(MFs):
    g1=[mf for mf in MFs if mf.body.color==Mig_Early]
    g2=[mf for mf in MFs if mf.body.color==Mig_Mid]
    g3=[mf for mf in MFs if mf.body.color==Mig_Late]
    #blinking(MFs, 100)
    t=0
    while True:    
        print('Time step', t)
        blinking2(g1, exponentiated_exponential_function(t)+1, color=color.cyan)
        blinking2(g2, exponentiated_exponential_function(t, alpha_=35)+1, color=color.magenta)
        blinking2(g3, exponentiated_exponential_function(t, alpha_=155)+1, color=color.blue)
        if t<40: t+=10
        if t>=40: t+=1
        if t>=300: 
            print('end of blink')
            sys.exit()


def intensity_mapping(MFs):
    print('generating intensity map')
    start_time = time.time()

    map_=np.zeros([height_PCL,area_length, area_width])

    max_range=int(radius_GC+radiation3)
    k=0
    for mf in MFs:        
        print(k,'-th mf')
        k+=1

        origin=mf.rosette.pos
        for x in range(-max_range, max_range+1):
            for y in range(-max_range, max_range+1):
                for z in range(-max_range, max_range+1):
                    distance=math.sqrt((x-origin.x)**2\
                                      +(y-origin.y)**2\
                                      +(z-origin.z)**2)
                    if distance<radius_GC+radiation1:
                        map_[x][y][z]+=intensity_radiation1
                    elif distance<radius_GC+radiation2:
                        map_[x][y][z]+=intensity_radiation2
                    elif distance<radius_GC+radiation3:
                        map_[x][y][z]+=intensity_radiation3
    print(type(map_))
    print(map_.shape)
    elapsed_time = time.time() - start_time
    print('elapsed time:', elapsed_time)
    return map_
            
            