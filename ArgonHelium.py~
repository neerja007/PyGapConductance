print("Enter Helium and Xenon concentration")
He_conc = float(input())
Ar_conc = float(input())

import cantera as ct
import math

print('Temperature and Pressure obtained from filmdrop : ')
t1 = float(input())
p1 = float(input())

print('Enter radius of fuel rod and temperature at the surface and temperature at center of fuel rod');
rf = float(input())
ts = float(input())
tc = float(input())

def cal_domeconductance(rf,ts):

  print('Enter radius at which temperature should be found (Mid point of fuel rod had r = 0)')
  rp = float(input())
  T = float( ts + ((tc-ts)*(1-(rp/rf)**2)))
  print T 
  
def cal_gapconductance(p1,t1):

  Ar = ct.Solution('Argon.xml')
  Ar. TP = t1, p1
  d= float(2*(((Ar.volume_mole/(6.023*10**23))*(3/(4*3.14)))**(1/3)))
  lamda_Ar=float((.0821*t1)/(math.sqrt(2)*math.pi*d*d*p1*6.023*10**23))

  He = ct.Solution('Silane.xml')
  He.X = 'HE:1'
  He. TP = t1, p1
  d= float(2*(((He.volume_mole/(6.023*10**23))*(3/(4*3.14)))**(1/3)))
  lamda_He=float((.0821*t1)/(math.sqrt(2)*math.pi*d*d*p1*6.023*10**23))

  #Argon calculation
  k_Ar= float(Ar.thermal_conductivity)
  Cp_Ar= float(Ar.cp)
  v_Ar= float(Ar.viscosity)
  Pr_Ar=float( (v_Ar*Cp_Ar/k_Ar))
  Cv_Ar= float(Ar.cv)
  S_Ar = float(Ar.viscosity*(273.15+Ar.T)/(273.15**3/2))
  #alpha of Ar- .85

  #Helium calculation
  k_He= float(He.thermal_conductivity)
  Cp_He= float(He.cp)
  v_He= float(He.viscosity)
  Pr_He= float((v_He*Cp_He/k_He))
  Cv_He= float(He.cv)
  S_He = float(He.viscosity*(273.15+He.T)/(273.15**3/2))
  #alpha of He- .25
  #lamda at 0C and 1 bar = 1.73*10^-5

  gamma = float(((He_conc*Cp_He)+(Ar_conc*Cp_Ar))/( (He_conc*Cv_He)+(Ar_conc*Cv_Ar)))

  alpha= float((((He_conc*.25/math.sqrt(He.mean_molecular_weight))+(Ar_conc*.25/math.sqrt(Ar.mean_molecular_weight)))/
               ( (He_conc/math.sqrt(He.mean_molecular_weight)) +(Ar_conc/math.sqrt(Ar.mean_molecular_weight)) )))

  S =  float((He_conc*S_He) + (Ar_conc*S_Ar))
               
  L = float( (1.73*10**-5)*(1+ (S/273))*He.T/( (He.P* (1+(S/He.T))*273)))

  K = float(He_conc*k_He + Ar_conc*k_Ar)      

  G_jump = float((15/4)* ( (2-.827*alpha)/alpha ) * L)
          
  # Assuming surface roughness to be 1
  G = float(G_jump+1)
  h= (K/G)
  print(h)

cal_gapconductance(p1,t1)
          
