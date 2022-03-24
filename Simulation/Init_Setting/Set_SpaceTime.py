#from function import *
#from ..Parameters import *
from vpython import *


class time_box:
    def __init__(self, pos, text):
        self.pos=pos
        self.text=text

class class_time_display:
    def __init__(self,area_length, height_PCL, time_division, time_exhibit='Embryonic stage',\
                 counter=0, vpython=True):
        #self.time_exhibit=time_exhibit
        self.counter=counter  #-(24*7* time_division)
        self.day=0
        self.start=False
        if vpython:
            self.display=label(pos=vector(area_length*2,height_PCL/2,0), text=time_exhibit)
        else:
            self.display=time_box(pos=vector(area_length*2,height_PCL/2,0), text=time_exhibit)

        self.time_division=time_division

    def time_count(self):
        if self.start==False:
            self.start=True
        else:
            self.counter+=1
            if (self.counter%(24*self.time_division)) ==0:
                self.day+=1
                self.display.text=('P%d'%self.day)
                #if self.day<15:
                #print('Day %d, h %d'%(self.day, self.counter/self.time_division))
                print('P%d'%self.day)
            if self.counter==20*24*self.time_division:
                print('Simulation time end, after processing....')


class layer_box:
    def __init__(self, pos):
        self.pos=pos

class moving_layer:
    def __init__(self, area_length, height_PCL, area_width, vpython=True):
        length=area_length/2
        width=area_width/2
        if vpython:
            self.expanding_border=box(pos=vector(length,height_PCL,width),\
                                size=vector(area_length,0.1,area_width),\
                                color=color.magenta)
        else:
            self.expanding_border=layer_box(pos=vector(length,height_PCL,width))
        
    def expand(self, depth):
        self.expanding_border.pos.y=depth   
        #self.expanding_border.pos.y-=height_PCL/(24*20)   


def layer_division(area_length, height_PCL, area_width, height_WM=0, labeling=False, vpython=True):
    length=area_length/2
    width=area_width/2
    if vpython:
        Intersection1=box(pos=vector(length,height_PCL,width), size=vector(area_length,0.1,area_width), color=color.blue)  ## btwn PCL & IGL    
        Intersection2=box(pos=vector(length,height_WM ,width), size=vector(area_length,0.1,area_width))  ## btwn IGL & WM
    else:
        Intersection1=layer_box(pos=vector(length,height_PCL,width))  ## btwn PCL & IGL    
        Intersection2=layer_box(pos=vector(length,height_WM ,width))  ## btwn IGL & WM


    if labeling==True and vpython==True:
        label(pos=vector(length,height_PCL,width), text='The bottom of PCL')
        label(pos=vector(length,height_WM,width), text='The bottom of IGL')

