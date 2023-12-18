import sys
import os.path


import matplotlib.pyplot as plt
with open("vacancy_coord.txt", 'r') as f:
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y1 = [float(line.split()[1]) for line in lines]
    y2 = [float(line.split()[2]) for line in lines]
    y3 = [float(line.split()[3]) for line in lines]

with open("solute_coord.txt", 'r') as f:
    lines = f.readlines()
    #x = [float(line.split()[0]) for line in lines]
    y4 = [float(line.split()[1]) for line in lines]
    y5 = [float(line.split()[2]) for line in lines]
    y6 = [float(line.split()[3]) for line in lines]

plt.plot(x,y1,"r--",x,y2,"g--",x,y3,"b--",
	x,y4,"r^",x,y5,"g^",x,y6,"b^")
#plt.title("disSV vs time")
plt.xlim(0,6e6)
plt.ylim(0,0.5e6)
plt.show()