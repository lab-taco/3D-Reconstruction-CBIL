import sys
import os
import numpy as np
from Simulation.Data_Extraction import load_coefficient
import matplotlib.pyplot as plt

DATA_FOLDER='Coefficient'
DATA_PATH = os.getcwd()
DATA_NAME_1= "Coefficient_overSwap1Large.npy"
DATA_NAME_2= "Coefficient_overSwap2Large.npy"
DATA_NAME_3= "Coefficient_overSwapmethod4large.npy"
DATA_NAME_3= "Coefficient_overSwapmethod4largeFullswap.npy"

#t_range, List_Assr_coeff_MF1, List_Assr_coeff_GC1 = load_coefficient(DATA_FOLDER, DATA_PATH, DATA_NAME_1)
#t_range, List_Assr_coeff_MF2, List_Assr_coeff_GC2 = load_coefficient(DATA_FOLDER, DATA_PATH, DATA_NAME_2)
t_range, List_Assr_coeff_MF3, List_Assr_coeff_GC3 = load_coefficient(DATA_FOLDER, DATA_PATH, DATA_NAME_3)
t_range=np.array(t_range)
t_range/=3000*2

plt.plot(t_range, List_Assr_coeff_MF3, 'r', label='MF')
plt.plot(t_range, List_Assr_coeff_GC3, 'b', label='GC')
plt.xlabel("Proportion of Swapped")
plt.ylabel("Coefficient")
plt.title("Coefficient over edge swapping")

"""
plt.plot(t_range, List_Assr_coeff_MF1, label='MF1')
plt.plot(t_range, List_Assr_coeff_GC1, label='GC1')
plt.plot(t_range, List_Assr_coeff_MF2, label='MF2')
plt.plot(t_range, List_Assr_coeff_GC2, label='GC2')
plt.xlabel("Num Swap")
plt.ylabel("Coefficient")
plt.title("Coefficient over edge swapping")
"""

#plt.plot(List_Assr_coeff_MF1, List_Assr_coeff_MF2, 'r', label='MF')
#plt.plot(List_Assr_coeff_GC1, List_Assr_coeff_GC2, 'b', label='GC')
#plt.xlabel("Num pickup")
#plt.ylabel("Num swap")
#plt.title("Pickup-swap rate")


plt.legend()
plt.show()

sys.exit()
data_set = (('A1', 'B50'), ('C75', 'A99'))

# Check if at least one element satisfies the condition
has_element = any(x.startswith('A') and int(x[1:]) < 100 
                  for pair in data_set 
                    for x in pair)

# Print the result
if has_element:
    print("At least one element starts with 'A' and has a number less than 100.")
else:
    print("No element satisfies the condition.")


print(not True)