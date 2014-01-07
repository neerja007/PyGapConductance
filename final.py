import cantera as ct
import math
import numpy
print("Enter Helium and Xenon concentration")
He_conc = float(input())
Ar_conc = float(input())


print("Enter bulk Temperature of Coolant - Tb(K)")
Tb = float(input())
print('Enter Pressure in pascals - ')
p1 = float( input() )
X = 45000

print("Enter outer radius of cladding,inner radius of cladding,radius of fuelrod")
Rco = float(input())
Rci = float(input())
Rf = float(input())

arr1 = [0,0,0,0]
arr2 = [0,0,0,0]
i = 0
i1 = 0
j = 1
#do..loop starts here
condition = True
while condition:
	def film_drop(Tb):
		print("Enter hg of coolant-Na")
		hg = float(input())
		Tco = (X/(2*(math.pi)*Rco*hg))+Tb
		return(Tco)
	Tco = film_drop(Tb)
	if(j%2!=0):	
	 arr1[0] = Tco
	else:
	 arr2[0] = Tco
	def clad_drop(Tco):
		if i=0:
		 print("Enter kc of cladding-Stainless steel")
		 kc = float(input())
		 Tci = ((X/(2*math.pi*kc))*2.303*math.log((Rco/Rci)))+Tco
		 i++
		else:
		 kcp = np.polyld(6.5181*10**-6,-5.105*10**-3,1.428,-4.7127)
		 Tci = ((X/(2*math.pi*kcp(Tco)))*2.303*math.log((Rco/Rci)))+Tco	
		return(Tci)
	Tci = clad_drop(Tco)
      	if(j%2!=0):	
	 arr1[1] = Tci
	else:
	 arr2[1] = Tci
 	
	def gap_drop(p1,Tci):
	
  	Ar = ct.Solution('Argon.xml')
  	Ar. TP = Tci, p1
  	d= float(2*(((Ar.volume_mole/(6.023*10**23))*(3/(4*3.14)))**(1/3)))
  	lamda_Ar=float((.0821*Tci)/(math.sqrt(2)*math.pi*d*d*p1*6.023*10**23))

 	 He = ct.Solution('Silane.xml')
  	He.X = 'HE:1'
  	He. TP = Tci, p1
	  d= float(2*(((He.volume_mole/(6.023*10**23))*(3/(4*3.14)))**(1/3)))
	  lamda_He=float((.0821*Tci)/(math.sqrt(2)*math.pi*d*d*p1*6.023*10**23))

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
 	 hg= (K/G)
  	Ts= (X/(hg*2*math.pi*Rf))+Tci
  	return(Ts)
	Ts = gap_drop(p1,Tci)
  	if(j%2!=0):	
	 arr1[2] = Ts
	else:
	 arr2[2] = Ts

	def fuel_drop(Ts):
		print("Enter k of fuelrod")
		k = float(input())
		Xn = X/((math.pi)*Rf**2)
		i1++;
		if i1 > 0:
		 k = (0.042*Ts)+((271*(10**-4)*Ts**2)/2)+((6.9*10**-11)*Ts**3/3)
		Tcl = np.polyld(1.725,0,135.5*10**-4,0.042,(Xn/(4*math.pi))+(0.042+135.5*10**-4+1.725)*Ts)
		print (Tcl.r)
		for num1 in range(0,3):
		 if(Tcl.r(num1)>0):
		  return(Tcl)
		 else:
		 print("Error")
	Tcl = fuel_drop(Ts)
	if(j%2!=0):	
	 arr1[3] = Tcl
	else:
	 arr2[3] = Tcl
	print('The center line temperature is ', Tcl)
	j++;
	for num in range(0,4):
		if((arr1(num)-arr2(num)<2)||(arr2(num)-arr1(num)<2)):
		 count++;
	if count==4:	
	condition = False
