from shapely.geometry import MultiPoint
from shapely.geometry import Point
from shapely.geometry import Polygon
import numpy as np
import math
from matplotlib import pyplot as plt
from astropy.stats import RipleysKEstimator

def my_K_func(data, Map_size_2D, cell_radius, function_type = 'K', graph=False, \
                return_L=False, rt_clut=False, cutoff=0, geo=False):
    x_axis, y_axis = Map_size_2D
    #for data_point in data:
    #    if data_point[0]>x_axis:
    #        raise Exception('migration misleading')
    #    if data_point[1]>y_axis:
    #        raise Exception('migration misleading')
    
    Convex_hull_ofdatapoints = MultiPoint(data).convex_hull
    #area=x_axis*(y_axis/2)
    area=Convex_hull_ofdatapoints.area
    #print(data)
    if len(data)==0: return 0
    lambda_inverse =  area/len(data) #density inverse

    r_min=cell_radius*2
    #r_max=sqrt(area/2)
    r_max=int(x_axis/2)
    r = np.linspace(r_min, r_max, r_max-r_min+1) #radii of inspection

    k=[]
    l=[]
    poi=[]
    print('calculating K function...')
    for t in r:        
        total_sum=0
        #for i in [point for point in list(enumerate(data))]: #enumerate(data) = index, value
        for ind, point in enumerate(data): #enumerate(data) = index, value
            
            #c2 = Point(point[0], point[1]).buffer(t)
            #intersection = cvh.intersection(c2)   
            #compensation = c2.area/intersection.area
            #print('compensation difference', compensation, compensation2)
            #print('c',c, c2.area, 'difference:',c-c2.area)
            #print('intersections', intersec, intersection.area, 'difference:',intersec-intersection.area, )

            #print('\n{:<6}'.format('intrsctn 1:'), '{:>6}'.format(np.round(intersec,3))                   ,\
            #        '{:<6}'.format('circle   1:'), '{:>6}'.format(np.round(c,3))                          ,\
            #        '{:<6}'.format('cmpnstn  1:'), '{:>6}'.format(np.round(compensation,3))               ,\
            #      '\n{:<6}'.format('intrsctn 2:'), '{:>6}'.format(np.round(intersection.area,3))          ,\
            #        '{:<6}'.format('circle   2:'), '{:>6}'.format(np.round(c2.area,3))                    ,\
            #        '{:<6}'.format('cmpnstn  2:'), '{:>6}'.format(np.round(compensation2,3))              ,\
            #      '\n{:<6}'.format('intrsctn d:'), '{:>6}'.format(np.round(intersec-intersection.area, 3)),\
            #        '{:<6}'.format('circle   d:'), '{:>6}'.format(np.round(c-c2.area,3))                  ,\
            #        '{:<6}'.format('cmpnstn  d:'), '{:>6}'.format(np.round(compensation-compensation2,3)) )
            
            #if intersec-intersection.area>2:
            #    g = gpd.GeoSeries([cvh])
            #    g.plot()
            #    plt.show()
            #    
            #    g = gpd.GeoSeries([cvh, c2])
            #    g.plot()
            #    plt.show()
            ##compensation=1 if no compensation
            #if ind==10:sys.exit()
            
            indicator_sum=0
            #for j in [point for point in list(enumerate(data)) if not point[0]==ind]:
            for ind2, point2 in [e for e in enumerate(data) if not e[0]==ind]:
                distance=np.sqrt((point[0]-point2[0])**2+(point[1]-point2[1])**2)
                #print('distance:',distance)
                if distance <= t:
                    indicator=1
                else:
                    indicator=0             
                indicator_sum+=indicator
            if indicator_sum>0:
                compensation_rate=edge_correction_convexhull(point ,t, Convex_hull_ofdatapoints, cell_radius/5)
                total_sum+=indicator_sum*compensation_rate            
        prob=total_sum/(len(data))
        k_t = lambda_inverse*prob

        if k_t<0:
            raise Exception('negative total sum for K Function of', t, 'radius')
        #homgen_k = np.math.pi*(t**2)
        #print('homogenous: %f, actual k: %f'% (homgen_k, k_t))
        k.append([t, k_t])
        #l.append([t, np.sqrt(k_t/math.pi)-t])
        l.append([t, np.sqrt(k_t/math.pi)-t])
        #l_poi=np.sqrt(Kest.poisson(r)/math.pi)
        
    k=np.asarray(k)
    l=np.asarray(l)

    #print('Spatial Calculation finished')
    
    Kest = RipleysKEstimator(area=area, x_max=x_axis, y_max=y_axis, x_min=-10, y_min=-10)
    #z = np.random.uniform(low=0, high=100, size=np.shape(data))
    
    
    if function_type=='K' and graph:
        print('Plotting K function')
        plt.plot(r, Kest.poisson(r), color='green', ls=':', label=r'$K_{pois}$')
        #same with plt.plot(r, math.pi*r**2)
        plt.plot(r, Kest(data, r, mode='ripley'), \
                color='black', label=r'$K_{ripley}$')    
        plt.plot(k[:,0], k[:,1], label=r'$K_{My_K}$') #K function

        #plt.plot(r, Kest(data=data, radii=r, mode='none'), color='red', ls='--', label=r'$K_{un}$')
        #plt.plot(r, Kest(data=data, radii=r, mode='translation'), color='black', ls=':', label=r'$K_{trans}$')
        #plt.plot(r, Kest(data=data, radii=r, mode='ohser'), color='blue', ls='-.', label=r'$K_{ohser}$')
        #plt.plot(r, Kest(data=data, radii=r, mode='var-width'), color='green',label=r'$K_{var-width}$')
        
    if function_type=='L' and graph:
        #print('Plotting L function')
        l_poi=np.sqrt(Kest.poisson(r)/math.pi)
        plt.plot(r, l_poi-r, color='green', ls=':', label=r'$L_{pois}$')    
        plt.plot(r, Kest.Lfunction(data,r), color='black', label=r'$L_{imported}$')
        plt.plot(l[:,0], l[:,1], label=r'$My L$') #L function    
    if graph: 
        plt.title('Spatial Analysis using '+function_type+'function')
        plt.legend()
        plt.show()

    if rt_clut:
        l_poi=np.sqrt(Kest.poisson(r)/math.pi)
        difference=l[:,1]-l_poi
        #print(difference)
        degree_of_clustering=np.round(np.mean(difference),3)
        print('The degree of clustering:', degree_of_clustering)        
        return degree_of_clustering
    if return_L:
        return l  #<< L not 1


import numpy as np
import math

def area_segment(k, r): # distance from the origin to intersection, radius
    #if k>r: print('k:', k, 'r:',r)
    #print('k:', k, 'r:',r)
    sin_alpha=math.sqrt(r**2-k**2)/r    
    angle=2*np.arcsin(sin_alpha)        
    tri=0.5*(r**2)*np.sin(angle)
    sect=0.5*(r**2)*angle
    #print('sector:', sect, 'tri:',tri)
    area=sect-tri
    return area
def area_segment_corner(u, v, r, decision_code=None): # relative distance from the point to the intersections, radius        
    y_D =math.sqrt(r**2-u**2)
    R=0.5*(y_D-v)*(math.sqrt(r**2-v**2)-u)
    angle=np.arcsin(y_D/r)-np.arcsin(v/r)
    #angle=(np.arcsin(y_D/r)-np.arcsin(v/r))*180/pi
    tri=0.5*(r**2)*np.sin(angle)
    sect=0.5*(r**2)*angle        
    #print('sector:', sect, 'tri:',tri)
    Y=sect-tri
    
    if R<0 or Y<0:
        print('r:',r, 'u:',u,'v:',v, 'y_D:', y_D)   
        print('tri:',tri,'sector:',sect, 'angle:', angle)
        print('R:',R, 'Y:',Y)
        print('code:', decision_code)
        raise Exception('negative RY')
    area=Y+R
    return area

def Euclidean_dist(x2, y2, x1, y1):
    return math.sqrt((x2-x1)**2+(y2-y1)**2)


def edge_correction_convexhull(data_point ,t, Convex_hull_ofdatapoints, simplification_tolerance):
    circle=Point(data_point).buffer(t)
    #simple_circle=circle.simplify(simplification_tolerance)
    
    #plt.plot(*Convex_hull_ofdatapoints.exterior.xy, label='conv hull')
    #plt.plot(*simple_circle.exterior.xy, label='circle')
    #plt.legend()
    #plt.show()
    computational_tolerance=math.pi*(t**2)-circle.area #min precision in circle area estimation
    if circle.intersects(Convex_hull_ofdatapoints):
        intersections = circle.intersection(Convex_hull_ofdatapoints).area
        #print('calc:', math.pi*(t**2), 'circle.area', circle.area) 
        #return math.pi*(t**2)/(math.pi*(t**2)-intersections)
        #print(circle.area, intersections)
        if intersections<=computational_tolerance: return 1
        else: return circle.area/intersections
    elif intersections<0:
        raise Exception('intersections:', intersections, 'Wrong calculation for edge correction')
    else: return 0
