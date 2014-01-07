import math
print("Enter bulk Temperature of Coolant - Tb(K)")
Tb = float(input())
X = 45000
print("Enter outer radius of cladding,inner radius of cladding,radius of fuelrod")
Rco = float(input())
Rci = float(input())
Rf = float(input())

#loop starts here
i = 0
j = 0
def film_drop(Tb):
	print("Enter hg of coolant-Na")
	hg = float(input())
	Tco = (X/(2*(math.pi)*Rco*hg))+Tb
	return(Tco)

def clad_drop(Tco):
	print("Enter kc of cladding-Stainless steel")
	hc = float(input())
	Tci = ((X/(2*math.pi*kc))*2.303*math.log((Rco/Rci)))+Tco
	return(Tci)

def gap_drop():

def fuel_drop(Ts):
	print("Enter k of fuelrod")
	k = float(input())
	X` = X/((math.pi)*Rf**2)
	Tcl = X`/(4*math.pi*k)+Ts
	#when i > 0
	#k = (0.042*Ts)+((271*(10**-4)*Ts**2)/2)+((6.9*10**-11)*Ts**3/3)
	return(Tcl)


