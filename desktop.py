print("Enter Helium and Xenon concetration")
He_conc = input()
Ar_conc = input()

import math

def alphafunc(conc,al,wt):
    alpha =  (conc*al/math.sqrt(wt) )/ (conc/math.sqrt(wt))
    return alpha

def freepath(p,t,d): #calculate mean free path (formula from net)
    l=(.0821*t)/(math.sqrt(2)*math.pi*d*d*p*6.023*10^23)
    return l


l1= freepath(p,t,d) #He
l2= freepath(p,t,d) #Ar
lambda_mix= (1.6*l1)*(1.6*l2)/(2*1.6)

alpha1= alpha(He_conc,al,4)
alpha2= alpha(Ar_conc,al,36)

alpha_mix= alpha1+alpha2

#k of He at 20C - .138
# K of Ar - .016
#Cp of He- 5190 J/Kg K
#Cp of Ar- 520 J/Kg K
#viscocity of He-.0000186
#viscocity of Ar-.000021

k_mix= (He_conc*.138) + ( Ar_conc*.016 )
Cp_mix= He_conc*5190 + Ar_conc*520
vis_mix= He_conc*.0000186 + Ar_conc*.000021
Pr_mix= (vis_mix* Cp_mix)/k_mix

G_jump = 2*((2-alpha_mix)/alpha_mix )*(1.6/(1+1.6))*(k_mix/(Pr_mix*Cp_mix) )*lamda_mix
