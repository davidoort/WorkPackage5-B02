from math import *
import matplotlib.pyplot as plt

L = 1.
v = 0.3
R = 0.5
t1 = 0.0001
E = 27.*10**9
I = 123456.
A = 0.9
p = 10000.


k = 2*sqrt(12*(L**4)*(1-(v**2))/((pi**4)*(R**2)*(t1**2)))
Q = (p/E)*((R/t1)**2)
stresscr = ((pi**2) * E * I / (A* (L**2)))
stressbuckl = ((1.983-0.983*(e**(-23.14*Q)))*k*(pi**2)*E*((t1/L)**2)/(12*(1-v**2)))

print stresscr, "MPa"
print stressbuckl, "MPa"


