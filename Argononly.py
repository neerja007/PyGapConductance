import cantera as ct
import math

print('Temperature and Pressure obtained from filmdrop : ')
t1 = float(input())
p1 = float(input())

gas = ct.Solution('Argon.xml')
gas. TP = t1, p1
d= float(2*(((gas.volume_mole/(6.023*10**23))*(3/(4*3.14)))**(1/3)))
lamda=float((.0821*t1)/(math.sqrt(2)*math.pi*d*d*p1*6.023*10**23))

k= gas.thermal_conductivity
Cp= gas.cp
v= gas.viscosity
Pr= (v*Cp/k)
Cv= gas.cv
gamma= Cp/Cv
#alpha of Ar- .85

G_jump= float(2*((2-.85)/.85)*(gamma/(gamma+1))*(k/(Pr*Cp))*lamda)
#assuming surface roughness to be 1

G= float(G_jump+1)

h = float(k/G)
print(h)









