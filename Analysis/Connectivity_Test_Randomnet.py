from degree_dist import *
from K_func import *
import numpy as np
from matplotlib import pyplot as plt
import os
import sys

num_MFs=40
#num_GC=4
num_GC=num_MFs*3
draw_dist=False
draw_net=False


node_GC, node_MF, ed_list = random_networks(num_GC, num_MFs)
#node_GC, node_MF, ed_list= input_processing('tolerance_170_MFs_edges')

#gc_color_map = migration_timing(node_GC, GC_colors)
#ed_list = shuffling(node_GC, node_MF, ed_list, gc_color_map)

ratio5=neuralnet4(node_GC, node_MF, ed_list,  \
    MF_subgroup='total', plot_frequency_dist=draw_dist, draw_net=draw_net,\
         shuffle=False, rt_ratio=True)
ratio5=neuralnet4(node_GC, node_MF, ed_list,  \
    MF_subgroup='total', plot_frequency_dist=draw_dist, draw_net=draw_net,\
         shuffle=True, rt_ratio=True)

#sys.exit()


dfference=broadness_difference(node_GC, node_MF, ed_list)
print(dfference)

CDF_comparison_of_shuffling(node_GC, node_MF, ed_list)

sys.exit()
    



ratio_list=[ratio1, ratio2, ratio3, ratio4]
                     
cdf_list=[]
for ratio in ratio_list:
    pdf1, pdf2, cdf1, cdf2 = plot_degree_dist(ratio, value_return2=True)
    cdf_list.append([cdf1, cdf2])
    plt.hist(pdf2, bins=pdf1, density=False, align='mid')
    plt.show()

log_range_= np.arange(-5, 5)
for cdf in cdf_list:
    #plt.plot(log_range_, cdf[1])
    plt.plot(cdf[0], cdf[1])
plt.legend(["Early", "Mid", "Late", "Shuffled"], loc="lower right")
plt.show()

sys.exit()


#neuralnet3(nd_list1, nd_list2, ed_list,  \
#    MF_subgroup='early', frequency_dist=draw_dist, draw_net=draw_net)
e=[]
m=[]
l=[]
for i in range(100):
    early, mid, late = neuralnet3(nd_list1, nd_list2, ed_list,  \
        MF_subgroup='early', frequency_dist=draw_dist, draw_net=draw_net)
    e.append(early)
    m.append(mid)
    l.append(late)

print('Mean #GC per MF  Early:',np.mean(e), 'Mid:',np.mean(m), 'Late:', np.mean(l))


#print('Early GC: %d,  Mid GC: %d,  Late GC: %d '%(int(num_GC/3), int(num_GC/3), num_GC-2*(int(num_GC/3))))
