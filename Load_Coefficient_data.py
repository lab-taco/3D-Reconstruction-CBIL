import sys
import os
from Simulation.Data_Extraction import load_coefficient
import matplotlib.pyplot as plt

DATA_FOLDER='Coefficient'
DATA_PATH = os.getcwd()

t_range, List_Assr_coeff_MF, List_Assr_coeff_GC = load_coefficient(DATA_FOLDER, DATA_PATH)


plt.plot(t_range, List_Assr_coeff_MF, 'r', label='MF')
plt.plot(t_range, List_Assr_coeff_GC, 'b', label='GC')
plt.xlabel("Num Swap")
plt.ylabel("Coefficient")
plt.title("Coefficient over edge swapping")
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