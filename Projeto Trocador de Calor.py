# -*- coding: utf-8 -*-
"""Projeto Trocador de Calor
Notas: 1) Considerando que o grupo é da turma 1
       2) Considerando que o material do tubo interno é silicone

"""
#
import numpy as np
from math import *
import pandas as pd

#Dados do Projeto

#Fluido Quente -  
W  = 0.05/2                   # Vazão do Fluido Quente [kg/s] - (0.05 L/s)
T1_in = 333.15               # Temperatura de Entrada do Fluido Quente [K] (60 °C)
T1_out = 323.15              # Temperatura de Saida do Fluido Quente [K] (50 °C)

#Fluido Frio - Tubos
w  = 0.05/2                     # Vazão do Fluido Frio [kg/s] - (0.05 L/s)
T2_in = 298.15                # Temperatura de Entrada do Fluido Frio [K] - (20 °C)
#2_out = T2                   # Temperatura de Saida do Fluido Frio [K] - (30 °C)


# Balanço de Energia
Cph = 3819.61              #Calor Especifico do Fluido Quente na Temperatura de Entrada [J/kg.K] 
Cpc = 4180                 #Calor Especifico do Fluido Frio na Temperatura Média [J/kg.K] 

# Cálculo da Temperatura de Saida do Fluido Quente [K]
T2_out = T2_in+(((W*Cph*(T1_in-T1_out))/(w*Cpc)))
print("Temperatura de saída",T2_out)

#Temperatura Média
T1m = (T1_in + T2_out)/ 2     # 328,15 °K
T2m = (T2_in + T2_out)/ 2     # 298,15 °K

#DT Média - Contra Corrente ( Mais eficiente)
DT_1 = T1_in - T2_out
DT_2 = T1_out - T2_in
print("Delta T1", DT_1)
print("Delta T2", DT_2)


# DTML - Média Logarítmica
DTML = (DT_1-DT_2)/np.log(DT_1/DT_2)
print("DTML",DTML)

# Taxa Transferencia de calor(Q)
qh = W*Cph*(T1_out - T1_in)    #Fluido Quente
qc = w*Cpc*(T2_out - T2_in)    #Fluido Frio
print("taxa de calor do Fluido Frio",qc)
print("taxa de calor do Fluido Quente",qh)

#Dados para o fluido Quente (Temperatura Média: 55°C )
p_h =  1022.789                     # Massa Especifica [kg/m3]
u_h =  1.74E-3                      # Viscosidade [Pa.s                 
k_h =  0.5894                       # Condutividade Termica [W/m.K] 
Pr_h = (u_h*Cph)/k_h

#Dados para o fluido frio (Temperatura Médi40a: 25°C )
p_c = 997                             # Massa Especifica [kg/m3]
u_c = 0.891E-3                        # Viscosidade [Pa.s]
k_c = 0.607                           # Condutividade Termica [W/m.K]
Pr_c = (u_c*Cpc)/k_c



# Dados dos Tubos (A definir)

Di = 0.015                    # Diametro Interno do Tubo [m]   https://www.lojanetlab.com.br/acessorios-para-laboratorios/mangueira-de-silicone/mangueira-de-silicone-pacote-com-5-metros?parceiro=7105&gad=1&gclid=Cj0KCQjw1_SkBhDwARIsANbGpFtcxm2rRcd6U2iFfhwV-H918OLlV7dbC3IEA_d15fMXVrMEJfrfsoQaApr0EALw_wcB
Do = Di                       # Considerando Esperssura do canudo muito pequena
E  = 3.8e-5                   # Rugosidade relativa [m]    https://www.semanticscholar.org/paper/Influence-of-surface-roughness-on-hydrophobic-of-in-Kurimoto-Azman/e9667cad07beffb6c37f54ad2a96905cb4aa3fa9
k  = 0.2                      # condutuvidade Termica do Material [W/m K]  https://www.intertronics.co.uk/wp-content/uploads/2016/11/TB2007-12-Thermally-Conductive-Silicones.pdf
Nt = 20.                     # Número de tubos
Lhp = 0.40                    # Comprimento do tubo [m]
Ds  = 0.050                    # Diametro do Tubo externo [m] 

#Tubo Interno (i)

Si = Nt*(pi*Di**2)/4
print("SI: ",Si)    
Gi = w/Si
print("Gi: ",Gi)               
Rei = (Di*Gi)/u_c       #Reynolds
print("Reynolds",Rei)

# Coeficiente de Transferência Convectivo para o tubo interno
   # -->  Correlação - Sieder-Tate

Nu_i = 1.86*(Rei*Pr_c*Di/Lhp)**(1/3)
print("Nusselt: ",Nu_i)
hi = (k_c*Nu_i)/Di
f_Fi = 0.41*(np.log(0.23*(E/Di)+6.9/Rei))**(-2)
print("Coeficiente de transferência do tubo",hi)


# Tubo Externo

De = Dh = (Ds**2 - Nt*Do**2)/(Ds + Nt*Do)
print("De: ",De)

So = (Ds**2 - Nt*Do**2)*(pi/4)
print("So: ",So)
Go = W/So
print("Go: ",Go)
Reo = (De*Go)/u_h
print("Reynolds do anulo: ",Reo)

# Coeficiente de Transferência Convectivo para o tubo externo
   # --> Correlação - Sieder-Tate

Nu_o = 1.86*(Reo*Pr_h*Do/Lhp)**(1/3)

print("Nusselt tubo externo: ",Nu_o)

ho = (k_h*Nu_o)/Do

print("Coeficiente de transferência de calor ânulo: ",ho)

f_Fo = 0.41*(np.log(0.23*(E/Do)+6.9/Reo))**(-2)


# Coeficiente Global de Transferencia de calor

Rdi = 0       
Rdo = 0

D_oi = Do/Di
UC = ((Do/2*k)*np.log(Do/Di)+1/ho+(Do/(Di*hi)))**(-1)           # Coeficiente Limpo
print("UC",UC)
UD = ((1/UC) + (Rdi*D_oi) + Rdo)**(-1)                          #Coeficiente Sujo
print("UD",UD)

# Área de Troca Térmica

AC = (qc/(UC*DTML))
print("àrea de troca térmica limpa: ", AC)

AD = (qc/(UD*DTML))
print("àrea de troca térmica suja: ", AD)

#Área Linear
Al = pi*Do*Nt
print("Área linear: ",Al)

L = AD/Al
print("Comprimento: ",L)

Nhp = L/(2*Lhp)       # N° de grampos
print ("N° de grampos: ",Nhp)
N_hpf = 1

At = Al*Lhp*Nhp      # Área total de troca térmica
print("Área total de troca térmica: ",At)

U = qc/(Al*DTML)
print("Coeficiente global: ",U)

#Rd = (UC-U)/(UC*U)
#print("Rd de projeto: ",Rd)


#Queda de Pressão - Tubo
k_i = 0.7 + (2000/Rei) + (17.78/Di)

Pi_dist = (2*(Gi**2)*Lhp*f_Fi)/(p_c)

Pi_loc = (Gi**2)*k_i/(2*p_c)

Pt = Pi_dist + Pi_loc
print("Queda de Pressão total do tubo: ",Pt)

#Queda de Pressão - Ânulo
    # --> Considerando que o ânulo é de seção reta
k_o = 0.7 + (2000/Reo) + (17.78/Dh)

Po_dist = (2*(Go**2)*Lhp*f_Fo)/(p_h)

Po_loc = (Go**2)*k_o/(2*p_h)

Po = Po_dist + Po_loc
print("Queda de Pressão total do Ânulo: ",Po)