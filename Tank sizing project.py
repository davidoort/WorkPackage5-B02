import math
import numpy as np

#std parameters (lengths in [mm])
vol = 1.244070691*(1000)**3      # mm^3
#t1 = 0.0013
t2 = 0.341 #mm
p = 1. #MPa
m = 1967.3 #dry+fuel= 1337+630.3 = 1967.3 kg
g = 9.81
rho = 4430.0/(1000**3) # kg/mm^3
gforce = 6
#material properties
E = 113.8*1000    #MPa
v = 0.342

def calculator(R,t1):
    #calculated parameters
    A = 2 * math.pi * R *t1 #mm^2
    I = math.pi*(R**3)*t1   #mm^4
    sphere = (4./3.)*math.pi*(R**3)   #mm^3
    L = 2*R + (vol-sphere)/(math.pi*(R**2)) #mm

    #launch loads
    stress = m*g*gforce/A #MPa = N/mm^2

    #formulas concerning buckling critical stresses

    k = 2*(12*(L**4)*(1-(v**2))/((math.pi**4)*(R**2)*(t1**2)))**0.5 #[-]
    Q = (p/E)*((R/t1)**2) #[-]

    column = (math.pi**2) * E * I / (A* (L**2)) #same units as stress (MPa)
    shell = (1.983-0.983*(math.e**(-23.14*Q)))*k*(math.pi**2)*E*((t1/L)**2)/(12*(1-v**2)) #MPa

    #Safety Margins
    #mscol = column/stress - 1
    #msshell = shell/stress - 1

    #Mass tank
    tank_mass = rho*(2*math.pi*R*t1*(L-2*R)+4*math.pi*(R**2)*t2)  #in kg


    #Vibrations (natural frequencies) by approximating the tank as a beam

    f_a=0.25*math.sqrt(A*E/(tank_mass*L/1000)) #>35Hz  #L has to be in meters here

    f_l = 0.560*math.sqrt(E*I/(tank_mass*L**3/1000))     #>12Hz
    
    return [stress,min(column,shell),tank_mass,f_a,f_l,L]




#printing values
#print column, "column"
#print shell, "shell"
#print stress
#print mscol, "margin column"
#print msshell, "margin shell"

#------------------------------------------------
#                ITERATION
#------------------------------------------------

#iteration range (in mm!)
rmin=500. #500. maximum length of the tank is 2000 mm
rmax=667. # DO NOT CHANGE THIS, this is the maximum value of R that holds the whole volume in a sphere 

t1min=0.1 #0.1
t1max=2.

points = 10

listR = np.linspace(rmin,rmax,points) # or alternatively np.arange(start,end,step)  
listt1= np.linspace(t1min,t1max,points)

stress_applied = []
stress_cr = []
mass = []
potential_t1 = []
potential_R = []
#ANALYTICAL WAY: DISCARD all combinations that have stress_cr < stress_applied or natural frequencies below vibration frequencies and then find the index of the combination with lowest mass

for R in listR:
    for t1 in listt1:
        calculations = calculator(R,t1)
        tank_mass = calculations[2]
        actual_stress = calculations[0]
        critical_stress = calculations[1]
        f_l = calculations[4]
        f_a = calculations[3]

        #MAKE CHECKS
        if  f_l >12 and f_a > 35 and critical_stress>actual_stress:
            potential_R.append(R)
            potential_t1.append(t1)
            mass.append(tank_mass)
            stress_applied.append(actual_stress)
            stress_cr.append(critical_stress)

pos = mass.index(min(mass))
print "Optimal thickness [mm]: ",potential_t1[pos]
print "Optimal radius [mm]",potential_R[pos]
print "L: ",calculator(potential_R[pos],potential_t1[pos])[5]

#VISUAL WAY: plots


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


fig = plt.figure()
ax = fig.gca(projection='3d')


# Plot the surface.
pltstress_applied = ax.plot_surface(potential_R, potential_t1, stress_applied, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
pltstress_cr =ax.plot_surface(potential_R, potential_t1, stress_cr, cmap=cm.coolwarm,
                      linewidth=0, antialiased=False)

# Customize the axes.
ax.set_xlim(rmin, rmax)
ax.set_ylim(t1min, t1max)
#ax.set_zlim(min(min(),min()), 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(pltstress_applied, shrink=0.5, aspect=5)
fig.colorbar(pltstress_cr, shrink=0.5, aspect=5)

plt.show()











