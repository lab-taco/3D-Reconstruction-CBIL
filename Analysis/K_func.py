import numpy as np
from math import *
from matplotlib import pyplot as plt
#from astropy.stats import RipleysKEstimator
#from function import *

def Rieman_sum(t, dist, radius):    #intersection, intersection to circumference, radius    
    num_fraction=dist-t+1   #each 1
    delta=(dist-t)/num_fraction
    range_=np.linspace(t, dist, num_fraction)        
    y=np.sqrt(radius**2 - range_**2)
    area=2*np.sum(delta*y)
    return area
    
def area_segment(k, r): # distance from the origin to intersection, radius
    #if k>r: print('k:', k, 'r:',r)
    #print('k:', k, 'r:',r)
    sin_alpha=sqrt(r**2-k**2)/r    
    angle=2*np.arcsin(sin_alpha)        
    tri=0.5*(r**2)*np.sin(angle)
    sect=0.5*(r**2)*angle
    #print('sector:', sect, 'tri:',tri)
    area=sect-tri
    return area
def area_segment_corner(u, v, r, decision_code=None): # relative distance from the point to the intersections, radius        
    y_D =sqrt(r**2-u**2)
    R=0.5*(y_D-v)*(sqrt(r**2-v**2)-u)
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
    return sqrt((x2-x1)**2+(y2-y1)**2)

def edge_correction(x_i, y_i, radius, height_PCL=height_PCL, area_length=area_length):
    decision_code =[False, False, False, False] 
    bottom=0
    ceiling=height_PCL

    #The decision_code represents x_over, y_over, x_below, y_below
    if area_length < x_i+radius:decision_code[0]=True  #x_over=True
    if ceiling < y_i+radius: decision_code[1]=True  #y_over=True
    if 0 > x_i-radius:          decision_code[2]=True  #x_below=True
    if bottom > y_i-radius:          decision_code[3]=True  #y_below =True

    #if   decision_code==[0,0,0,0]: return 1,pi*(radius**2), pi*(radius**2)
    if   decision_code==[0,0,0,0]: return 1
    elif decision_code==[1,0,0,0]: area_s= area_segment(area_length-x_i, radius)
    elif decision_code==[0,1,0,0]: area_s= area_segment(ceiling-y_i, radius)
    elif decision_code==[0,0,1,0]: area_s= area_segment(x_i, radius)
    elif decision_code==[0,0,0,1]: area_s= area_segment(y_i, radius) 
    #or 2 segments whether a corner is superposed or not
    elif decision_code==[1,1,0,0]:
        area_s=area_segment(area_length-x_i, radius)+area_segment(ceiling-y_i, radius)
        if Euclidean_dist(x_i, y_i, area_length, ceiling)<radius:
            superposed=area_segment_corner(area_length-x_i, ceiling-y_i, radius, decision_code)
            area_s-=superposed
    elif decision_code==[0,1,1,0]:
        area_s=area_segment(ceiling-y_i, radius)+area_segment(x_i, radius)
        if Euclidean_dist(x_i, y_i, 0, ceiling)<radius:
            superposed=area_segment_corner(x_i, ceiling-y_i, radius, decision_code) #decision_code
            area_s-=superposed                
    elif decision_code==[0,0,1,1]:
        area_s=area_segment(x_i, radius)+area_segment(y_i, radius)
        if Euclidean_dist(x_i, y_i, 0, bottom)<radius:
            superposed=area_segment_corner(x_i, y_i, radius, decision_code)
            area_s-=superposed
    elif decision_code==[1,0,0,1]:
        area_s=area_segment(area_length-x_i, radius)+area_segment(y_i, radius)
        if Euclidean_dist(x_i, y_i, area_length, bottom)<radius:
            superposed=area_segment_corner(area_length-x_i, y_i, radius, decision_code)
            area_s-=superposed
    #or 2 segments with no corner
    elif decision_code==[1,0,1,0]:
        area_s= area_segment(area_length-x_i, radius)+area_segment(x_i, radius)
    elif decision_code==[0,1,0,1]:
        area_s= area_segment(ceiling-y_i, radius)+area_segment(y_i, radius) 
    #or 3 segments with 2 corners
    elif decision_code==[1,1,1,0]:
        area_s= area_segment(area_length-x_i, radius)\
                +area_segment(ceiling-y_i, radius)\
                +area_segment(x_i, radius)
        if Euclidean_dist(x_i, y_i, area_length, ceiling)<radius:
            superposed=area_segment_corner(area_length-x_i, ceiling-y_i, radius, decision_code)
            area_s-=superposed
        if Euclidean_dist(x_i, y_i, 0, ceiling)<radius:
            superposed=area_segment_corner(x_i, ceiling-y_i, radius, decision_code) #decision_code
            area_s-=superposed
    elif decision_code==[1,1,0,1]:
        area_s= area_segment(area_length-x_i, radius)\
                +area_segment(ceiling-y_i, radius)\
                +area_segment(y_i, radius)
        if Euclidean_dist(x_i, y_i, area_length, ceiling)<radius:
            superposed=area_segment_corner(area_length-x_i, ceiling-y_i, radius, decision_code)
            area_s-=superposed
        if Euclidean_dist(x_i, y_i, area_length, 0)<radius:
            superposed=area_segment_corner(area_length-x_i, y_i, radius, decision_code)
            area_s-=superposed
    elif decision_code==[1,0,1,1]:
        area_s= area_segment(area_length-x_i, radius)\
                +area_segment(x_i, radius)\
                +area_segment(y_i, radius)
        if Euclidean_dist(x_i, y_i, 0, 0)<radius:
            superposed=area_segment_corner(x_i, y_i, radius, decision_code)
            area_s-=superposed
        if Euclidean_dist(x_i, y_i, area_length, 0)<radius:
            superposed=area_segment_corner(area_length-x_i, y_i, radius, decision_code)
            area_s-=superposed
    elif decision_code==[0,1,1,1]:
        area_s= area_segment(ceiling-y_i, radius)\
                +area_segment(x_i, radius)\
                +area_segment(y_i, radius)
        if Euclidean_dist(x_i, y_i, 0, ceiling)<radius:
            superposed=area_segment_corner(x_i, ceiling-y_i, radius, decision_code) 
            area_s-=superposed
        if Euclidean_dist(x_i, y_i, 0, 0)<radius:
            superposed=area_segment_corner(x_i, y_i, radius, decision_code)
            area_s-=superposed
    #or 4 segments
    elif decision_code==[1,1,1,1]:
        area_s= area_segment(area_length-x_i, radius)\
                + area_segment(ceiling-y_i, radius)\
                + area_segment(x_i, radius)\
                + area_segment(y_i, radius)
        if Euclidean_dist(x_i, y_i, area_length, ceiling)<radius:
            superposed=area_segment_corner(area_length-x_i, ceiling-y_i, radius, decision_code)
            area_s-=superposed
        if Euclidean_dist(x_i, y_i, 0, ceiling)<radius:
            superposed=area_segment_corner(x_i, ceiling-y_i, radius, decision_code)
            area_s-=superposed
        if Euclidean_dist(x_i, y_i, 0, 0)<radius:
            superposed=area_segment_corner(x_i, y_i, radius, decision_code)
            area_s-=superposed
        if Euclidean_dist(x_i, y_i, area_length, 0)<radius:
            superposed=area_segment_corner(area_length-x_i, y_i, radius, decision_code)
            area_s-=superposed
            
    else:            
        raise Exception('Wrong code:', decision_code, 'x:',x_i, 'y:', y_i, 'r', radius)
    #print('weight:', pi*(radius**2)/area_s)
    if area_s<0:
        print('area_s:', area_s, 'code:', decision_code)
        raise Exception('Wrong calculation')
    #return pi*(radius**2)/(pi*(radius**2)-area_s), pi*(radius**2), area_s
    return pi*(radius**2)/(pi*(radius**2)-area_s)

'''
#from K_func import *
import geopandas as gpd
import geoplot'''
from shapely.geometry import MultiPoint
from shapely.geometry import Point
def edge_correction2(point):
    pass

def my_K_func(data, graph=False, rt_L=False, rt_clut=False, cutoff=0, geo=False):
    for data_point in data:
        if data_point[0]>area_length:
            raise Exception('migration misleading')
        if data_point[1]>height_PCL:
            raise Exception('migration misleading')
    
    mp=MultiPoint(data)
    cvh = mp.convex_hull

    area=area_length*(height_PCL/2)
    #print(data)
    if len(data)==0: return 0
    lambda_inverse =  area/len(data) #density inverse

    r_min=radius_GC*2
    #r_max=sqrt(area/2)
    r_max=area_length-12
    r = np.linspace(r_min, r_max, r_max-r_min+1) #radius

    k=[]
    l=[]
    poi=[]
    #print('calculating K function...')
    for t in r:        
        total_sum=0
        #for i in [point for point in list(enumerate(data))]: #enumerate(data) = index, value
        for ind, point in enumerate(data): #enumerate(data) = index, value
            compensation=edge_correction(point[0],point[1] ,t, height_PCL=height_PCL-cutoff)
            
            #compensation, c, intersec = edge_correction(point[0],point[1] ,t, height_PCL=height_PCL-cutoff)
            
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

            if compensation<0:
                print('negative compensation:', compensation)
                raise Exception('negative compensation')
            indicator_sum=0
            #for j in [point for point in list(enumerate(data)) if not point[0]==ind]:
            for ind2, point2 in [e for e in enumerate(data) if not e[0]==ind]:
                distance=np.sqrt((point[0]-point2[0])**2+(point[1]-point2[1])**2)
                #print('distance:',distance)
                if distance < t:
                    indicator=1
                else:
                    indicator=0             
                indicator_sum+=indicator              
            total_sum+=indicator_sum*compensation            
        prob=total_sum/(len(data))
        k_t = lambda_inverse*prob

        if k_t<0:
            print('total sum:', total_sum)
            raise Exception('negative total sum')
        #homgen_k = np.pi*(t**2)
        #print('homogenous: %f, actual k: %f'% (homgen_k, k_t))
        k.append([t,k_t])
        l.append([t,np.sqrt(k_t/pi)-t])
        
    k=np.asarray(k)
    l=np.asarray(l)

    #print('Spatial Calculation finished')
    
    Kest = RipleysKEstimator(area=area, x_max=area_length, y_max=height_PCL, x_min=-10, y_min=-10)
    #z = np.random.uniform(low=0, high=100, size=np.shape(data))
    
    
    show_k=False
    show_l=True
    if show_k and graph:
        print('Plotting K function')
        plt.plot(r, Kest.poisson(r), color='green', ls=':', label=r'$K_{pois}$')
        #same with plt.plot(r, pi*r**2)
        #plt.plot(r, Kest(data=data, radii=r, mode='ripley'), color='black', label=r'$K_{z_ ripley}$')    
        plt.plot(k[:,0], k[:,1]) #K function

        #plt.plot(r, Kest(data=data, radii=r, mode='none'), color='red', ls='--', label=r'$K_{un}$')
        #plt.plot(r, Kest(data=data, radii=r, mode='translation'), color='black', ls=':', label=r'$K_{trans}$')
        #plt.plot(r, Kest(data=data, radii=r, mode='ohser'), color='blue', ls='-.', label=r'$K_{ohser}$')
        #plt.plot(r, Kest(data=data, radii=r, mode='var-width'), color='green',label=r'$K_{var-width}$')
        if graph: plt.show()
    if show_l and graph:
        #print('Plotting L function')
        l_poi=np.sqrt(Kest.poisson(r)/pi)
        plt.plot(r, l_poi-r, color='green', ls=':', label=r'$L_{pois}$')    
        #plt.plot(r, Kest.Lfunction(data,r), color='black', label=r'$L_{data}$')
        plt.plot(l[:,0], l[:,1]) #L function    
        if graph: plt.show()

    if rt_clut:
        l_poi=np.sqrt(Kest.poisson(r)/pi)
        difference=l[:,1]-l_poi
        #print(difference)
        degree_of_clustering=np.round(np.mean(difference),3)
        print('The degree of clustering:', degree_of_clustering)        
        return degree_of_clustering
    if rt_L:
        return l[:, 1]

#def squaring(data, height_range=[0, height_PCL]):

def f_cutoff(data, cutoff_val):    
    for points in data:
        #print('eq',points[1],'-',cutoff_val,'=',points[1]-cutoff_val)
        points[1]-=cutoff_val
    return data


def squaring(data, ceiling=height_PCL, bottom=0):
    
    list_del=[]
    for i in enumerate(data):
        #if i[1][0]>area_length or i[1][1]>height_range[-1]:
        if i[1][0]>area_length or i[1][1]>ceiling:
            list_del.append(i[0])
        #if i[1][0] <0  or i[1][1]<height_range[0]:
        if i[1][0] <0  or i[1][1]<bottom:
            list_del.append(i[0])
    return np.delete(data, list_del,0)

def z_rounding(data):
    for subset in data:
        if len(subset)>0: subset[:,2]=np.round(subset[:,2])
    return data

def printnum(data):        
    print('Early GCs', len(data[0]))
    print('Mid GCs', len(data[1]))
    print('Late GCs', len(data[2]))
    print('Num total cells', len(data[0])+len(data[1])+len(data[2]))
    print('MFs', len(data[-1]))   

def total_GCs_concat(data):
    ttl_data=[]
    for sub_set in data:
        if len(sub_set)>0: ttl_data.append(sub_set)
    return np.concatenate(ttl_data)

def cells_in_scope(origin, data, distance):
    data_in_scope = [cell for cell in data if not np.all(cell==origin)\
                                            and cell[0]<origin[0]+distance and cell[0]>origin[0]-distance\
                                            and cell[1]<origin[1]+distance and cell[1]<origin[1]-distance]
    return data_in_scope


def nearest_neighbor_dist(data, graph=False, rt_moments=False, rt_cum=False,  title=''):    
    nearest_distance_distribution=[]
    for origin_cell in data:
        for range_of_scope in range(1, height_PCL):
            inspect = cells_in_scope(origin_cell, data, range_of_scope)
            if len(inspect)>0:
                #print(len(inspect), 'cells are caugt at range', range_of_scope)
                distance_distribution=[]
                for cell in inspect:
                    distance=Euclidean_dist(cell[0], cell[1], origin_cell[0], origin_cell[1])                    
                    distance_distribution.append(distance)
                the_nearest_dist = np.amin(distance_distribution)
                nearest_distance_distribution.append(the_nearest_dist)
                break
    mean = np.round(np.mean(nearest_distance_distribution),3)
    std  = np.round(np.std(nearest_distance_distribution), 3)
    from statsmodels.distributions.empirical_distribution import ECDF
    ecdf = ECDF(nearest_distance_distribution)
    if rt_cum:
        return ecdf.x, ecdf.y

    if graph:
        bin_range=np.arange(0, height_PCL)
        plt.hist(nearest_distance_distribution, bins=bin_range, density=False, align='mid')
        plt.title(title)
        plt.legend(title='Mean:{}  Std:{}'.format(str(mean) ,str(std)), loc="upper right")
        #plt.legend(title='Mean:{}  Std:{}'.format(str(mean) ,str(std)), loc="upper right", fontsize='xx-large')
        plt.show()

        from statsmodels.distributions.empirical_distribution import ECDF
        ecdf = ECDF(nearest_distance_distribution)        
        plt.plot(ecdf.x, ecdf.y)
        plt.show()

        

    
    if rt_moments:
        return mean, std


#def recon_2d_slice(stained_GCs, other_GCs, mfs):
def recon_2d_slice(load, load2,load3, load4):

    radius_scale=1
    cr = shapes.circle(radius=1*radius_scale, thickness=0.1)
    cr_hole = shapes.circle(radius=0.9*radius_scale)
    cr2 = shapes.circle(radius=1*radius_scale)
    cr3 = shapes.circle(radius=0.5*radius_scale)


    points=[]
    for i in load:
        ext=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr2, color=color.green)
        points.append(ext)
    print('Injected GCs done')

    for i in load2:
        ext1=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr_hole, color=color.black)
        ext2=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr, color=color.red)
        ext=compound([ext1, ext2])
        points.append(ext)
    print('Normal GC1s done')
    for i in load3:
        ext1=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr_hole, color=color.black)
        ext2=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr, color=color.red)
        ext=compound([ext1, ext2])
        points.append(ext)
    print('Normal GC2s done')

    for i in load4:
        ext=extrusion(path=[vec(i[0],i[1],i[2]), vec(i[0],i[1], i[2]+0.1)], shape=cr3, color=color.red)
        points.append(ext)
    print('MFs done')    

                

def regular_scattering(length, Map_size):
    leng=length
    arr=[]
    for i in range(0, leng):
        for j in range(0, leng):
            arr.append([42*i/leng, 42*j/leng])
    
    print(np.shape(arr))
    arr=np.asarray(arr)
    plt.scatter(arr[:,0], arr[:,1])
    plt.xlim(-1,42)
    plt.ylim(-1,42)
    plt.show()
    return arr

def random_scattering(num_cells, graph=False, exclusion_zone=False, append_arr=[]):
    N=num_cells
    arr=[]
    if len(append_arr)>0:
        arr=np.ndarray.tolist(append_arr)

    for i in range(0, N):            
        x = np.random.random_integers(0, 42)
        y = np.random.random_integers(0, 42)        
        if exclusion_zone and len(arr)>0:
            exclusion_range=2.5
            excls_count=0
            max_loop=100
            while(excls_count != len(arr)):                                
                for cell in arr:
                    distance=Euclidean_dist(cell[0], cell[1], x, y)                    
                    if distance < exclusion_range:                         
                        break
                    else:
                        excls_count+=1
                if excls_count != len(arr):                    
                    x = np.random.random_integers(0, 42)
                    y = np.random.random_integers(0, 42)
                    excls_count=0
                max_loop+=1
                if max_loop>10000:
                    raise Exception('Too many cells for regular distribution')
            #print('excls count', excls_count, 'num points:', len(arr))
                            
        arr.append([x, y])
        if False:        
            arr2=np.asarray(arr)
            plt.scatter(arr2[:,0], arr2[:,1])
            plt.xlim(-1,43)
            plt.ylim(-1,43)
            plt.show()
    arr=np.asarray(arr)
    if graph:        
        plt.scatter(arr[:,0], arr[:,1])
        plt.xlim(-1,42)
        plt.ylim(-1,42)
        plt.show()   
    
    return arr

def clustering_with_shuffle(arr, arr2):
    size_of_scope=10
    for i in range(200): #pick 100 pairs and swap to cluster
        indices1=np.random.choice(arr.shape[0],1)
        indices2=np.random.choice(arr2.shape[0],1)        
        pair=np.concatenate((arr[indices1], arr2[indices2]))

        num_inboundcell_for_pair=[]
        for target_cells in pair:
            numcell_inbound=0
            for cell in arr2:                
                if Euclidean_dist(target_cells[1], target_cells[0],cell[1], cell[0])<size_of_scope:
                    numcell_inbound+=1                    
            num_inboundcell_for_pair.append(numcell_inbound)
        if num_inboundcell_for_pair[0]>num_inboundcell_for_pair[1]:            
            arr =np.delete(arr, np.where(np.all(arr == pair[0], axis=1) ), axis=0)
            arr2=np.delete(arr2,np.where(np.all(arr2== pair[1], axis=1) ), axis=0)            
            
            arr  = np.append(arr , [pair[1]], axis=0)
            arr2 = np.append(arr2, [pair[0]],  axis=0)

    return arr, arr2

        
