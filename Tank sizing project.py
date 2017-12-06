from math import *
import matplotlib.pyplot as plt
#std parameters
vol = 1.244070691
R = 0.6
t1 = 0.0013
t2 = 0.341 #mm
p = 1000.
m = 1967.3 #dry+fuel= 1337+630.3 = 1967.3
g = 9.81
rho = 4.43 # g/mm^3k
gforce = 6.5
#material properties
E = 27.*10**9
v = 0.3

#calculated parameters
A = 2 * pi * R *t1
I = pi*(R**3)*t1
L = 1.5

#launch loads
stress = m*g*gforce/A

#formulas concerning buckling
k = 2*(12*(L**4)*(1-(v**2))/((pi**4)*(R**2)*(t1**2)))**0.5

Q = (p/E)*((R/t1)**2)

column = ((pi**2) * E * I / (A* (L**2)))

shell = ((1.983-0.983*(e**(-23.14*Q)))*k*(pi**2)*E*((t1/L)**2)/(12*(1-v**2)))


mscol = column/stress -1
msshell = shell/stress - 1

tank_mass = rho*(2*pi*R*t1*h_cil+4*pi*R**2*t2)


#printing values
print column, "column"
print shell, "shell"
print stress
print mscol, "margin column"
print msshell, "margin shell"





