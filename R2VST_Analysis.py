import sys
import os.path

#fh addition
r2output=open("r2output.txt","w")
vacancy_coord=open("vacancy_coord.txt","w")

if len(sys.argv) != 2:                               #fh: if you are not giving an argument, print this
    print("Input history.dis and output time evolution of r^2 and r1*r2")
    exit("Usage: %s <history.dis>" % sys.argv[0])

# global variables
vbra= ((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5));    #fh: bravais lattice coordinates
iltcp= []                                                       #fh: custom index?
x= []                                                           #fh: x,y,z coordinates of atom/specie
y= []
z= []
# global variables

def cal_dis2(i, j, k):                                          #fh: calculate distance from origin
    x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0]
    y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1]
    z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2]

    return x*x + y*y + z*z

def cal_dot(i1, j1, k1, i2, j2, k2):                            #fh: calculate dot product
    x1= i1*vbra[0][0] + j1*vbra[1][0] + k1*vbra[2][0]
    y1= i1*vbra[0][1] + j1*vbra[1][1] + k1*vbra[2][1]
    z1= i1*vbra[0][2] + j1*vbra[1][2] + k1*vbra[2][2]
    x2= i2*vbra[0][0] + j2*vbra[1][0] + k2*vbra[2][0]
    y2= i2*vbra[0][1] + j2*vbra[1][1] + k2*vbra[2][1]
    z2= i2*vbra[0][2] + j2*vbra[1][2] + k2*vbra[2][2]
    #print(x1,y1,z1,x2,y2,z2)
    return x1*x2 + y1*y2 + z1*z2

vcc= []
sol= [0, 0, 0]
#with open(sys.argv[1], 'r') as IN: # Read his.dis
with open("history.dis", 'r') as IN: # Read his.dis
    l= 0; N= 99999; 
    dis2D= 0; dis2S= 0; disDOT= 0 
    for line in IN:
        l +=1
        
        if   l==1: 
            N= int(line)                                #fh: if it is first line, gather number of species
            #fh: interesting point: the N for number of species is flexible for each data entry
            #print("number of things: ",N)
        elif l==2: (T, step, time1_, time2_, time3_) = line.split()    #fh: if it is second line, gather temp, step, time, and time2
        else:
            (tp, ltcp, x, y, z) = line.split()                       #fh: if it is third line, gather type and coordinates
            if l==3 and int(tp)==-1: exit("Error: first atom is not a defect")   #fh: if first atom is a B atom
            if l >3 and int(tp)!=-1: exit("Error: a not-first atom is not a solute")   #fh: if anything after first atom is not a B atom
                                                   #fh: this implies that he wants vacancy first followed by B atoms
 
            x= int(x)                            #fh: gather coordinates of te specie
            y= int(y)
            z= int(z)
            if l==3:                            #fh: if you have a vacancy, calculate its distance and add to this variable dis2D
                dis2D += cal_dis2(x, y, z)                 #also add its coordinates to a vcc array
                vcc= [x, y, z]
                vacancy_coord.write(str(step))
                vacancy_coord.write(" ")
                vacancy_coord.write(str(x))
                vacancy_coord.write(" ")
                vacancy_coord.write(str(y))
                vacancy_coord.write(" ")
                vacancy_coord.write(str(z))
                vacancy_coord.write(" ")
                vacancy_coord.write(str(x*x+y*y+z*z))
                vacancy_coord.write("\n")
            else:
                dis2S  += cal_dis2(x, y, z)     #fh: if you have a B atom, calculate its distance and add to this varialbe dis2S
                sol[0] += x; sol[1] += y; sol[2] += z      #also add its x, y, and z coordinates to sol array

        if (l-2)==N:                                 #fh: if you have reached to last specie/ line in that time step
            #print(sol[0],sol[1],sol[2])
            disSV= cal_dot(sol[0], sol[1], sol[2], vcc[0], vcc[1], vcc[2])       #fh: calculate dot product of B atom and vacancy
            disSS= cal_dot(sol[0], sol[1], sol[2], sol[0], sol[1], sol[2])       #fh: dot product of B atom with B atom
            #print(disSV)
            if N==1: 
                print(step, float(time1_), dis2D, "END")
                #r2output.write(step,float(time1_),dis2D,"END")
            else:    
                #print(step, float(time1_), dis2D, dis2S/(N-1), disSV, disSS)
                #print(dis2S,"\n")
                r2output.write(step)
                r2output.write(" ")
                r2output.write((time1_))
                r2output.write(" ")
                r2output.write(str(dis2D))
                r2output.write(" ")
                r2output.write(str(dis2S/(N-1)))
                r2output.write(" ")
                r2output.write(str(disSV))
                r2output.write(" ")
                r2output.write(str(disSS))
                r2output.write(" ")
                r2output.write("\n")
            
            l= 0; dis2D= 0; dis2S= 0; disDOT= 0; vcc= []; sol= [0, 0, 0]          #fh: reset variables defined earlier
            
#frank addition
r2output.close()
vacancy_coord.close()


