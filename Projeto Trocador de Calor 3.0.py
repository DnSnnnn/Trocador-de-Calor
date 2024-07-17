"Projeto Trocador de Calor"

"1° Modelagem"


# Bibliotecas
from scipy.optimize import fsolve
from matplotlib.pyplot import*
from numpy import*
#from pandas import DataFrame
import pandas as pd




# Dados

# Balanço de Energia
Cp1 = 3819.61  #J/kg*K
Cp2 = 4180     #J/kg* K

# Dados para o fluido Quente
p_1 = 1022.789   #kg/m3
u_1 = 1.74e-3    #Pa.s
k_1 = 0.5894     #w/m.K
Pr_1 = (u_1 * Cp1) / k_1

# Dados para o fluido Frio
p_2 = 997        #Kg/m3
u_2 = 0.891e-3   #Pa.s
k_2 = 0.607      #W/m.k
Pr_2 = (u_2 * Cp2) / k_2

# Dados dos Tubos
Di = 0.015   #m
Do = Di
E = 3.8e-5
K = 0.2     #W/m.K
N_t = int(20)
Lhp = 0.40  #m
Ds = 0.100  #m

# Fluido Quente
T1_in = 333.15     #K
T1_out =  323.15   #K
W_1 = .05/2        # L/s
Cp_1 = Cp1         #J/kg*K
Pr_1 = Pr_1
k_1 = k_1          #w/m.k
ro_1 = 1022.789    #Kg/m3
M_1 = 1.74e-3      #Pa.S

# Fluido Frio
T2_in = 298.15   #K
T2_out = 0
W_2 = .05/2        # L/s
Cp_2 = 4180      #J/kg.K
Pr_2 = Pr_2
k_2 = .607       #W/m.K
ro_2 = 997       #Kg/m3
M_2 = 0.891e-3   #Pa.s

W_i = W_2
Cp_i = Cp_2
Pr_i = Pr_2
k_i = k_2
ro_i = ro_2
M_i = M_2

W_o = W_1
Cp_o = Cp_1
Pr_o = Pr_1
k_o = k_1
ro_o = ro_1
M_o = M_1

# Balanço de Energia

def balanco_de_energia(T2_out):
    q1 = W_1*Cp_1*(T1_out-T1_in)
    q2 = W_2*Cp_2*(T2_out-T2_in)
    return (q1+q2)

T2_out = float(fsolve(balanco_de_energia,400))

q1 = W_1*Cp_1*(T1_out-T1_in)
q2 = W_2*Cp_2*(T2_out-T2_in)


# Delta Tm

delta_T1 = T1_in-T2_out
delta_T2 = T1_out-T2_in
delta_Tm = (delta_T1-delta_T2)/log(delta_T1/delta_T2 )


# Coeficientes de Transferência de Calor

S_i = (pi*N_t/4)*(Di*Di)
G_i = W_i/S_i
Re_i = (G_i*Di/M_i)
fF_i = .41*((log(.23*((E/Di)**(10/9))+6.9/Re_i))**(-2))
K1_i = 13.6*fF_i+1
K2_i = 11.7+1.8/((Pr_i)**(1/3))
Nu_i = 1.86*(Re_i*Pr_i*Di/Lhp)**(1/3)
h_i = (k_i*Nu_i)/Di

De = Dh = (Ds*Ds-Do*Do*N_t)/(Ds+Do*N_t)

S_o = (pi/4)*(Ds*Ds-Do*Do*N_t)
G_o = W_o/S_o
Re_o = G_o*De/M_o

Nu_o = 1.86*(Re_o*Pr_o*Do/Lhp)**(1/3)
h_o = (k_o*Nu_o)/Do
fF_o = .41*((log(.23*((E/Do)**(10/9))+6.9/Re_o))**(-2))

K1_o = 13.6*fF_o+1
K2_o = 11.7+1.8/((Pr_o)**(1/3))

# Comprimentos de Entrada

LfHd_i = .05*Re_i*Di*Pr_i
LfHd_o = .05*Re_o*Do*Pr_o


# Coeficientes Globais de Transferẽncia de Calor

Uclean = ((Do/(2*K))*(log(Do/Di))+1/h_o+(1/h_i)*(Do/Di))**(-1)
Uc = Uclean
Rd_i = 0.0001
Rd_o = 0.0001
Ud = (1/Uc+(Rd_i*(Do/Di))+Rd_o)**(-1)

# Área Linear

Alinear = pi*Do*N_t

# Área Clean

Aclean = abs(q1)/(Uc*delta_Tm)
Ac = Aclean

# Área Dirty

Adirty = abs(q1)/(Ud*delta_Tm)
Ad = Adirty
L = Ad/Alinear
Nhp = (L/(2*Lhp))
Areal = 2*Alinear*Nhp*Lhp
A = Areal

#Real

q = abs(q1)
U = q/(A*delta_Tm)
Rd = (Uc-U)/(Uc*U)
delta_P_i_dist = ((4*G_i*G_i*Nhp*Lhp*fF_i)/(ro_i*Di))
K_i = .7 +2000/Re_i +17.78/(Di*1000)
delta_P_i_loc = (G_i*G_i*K_i*(2*Nhp-1))/(2*ro_i)
delta_P_i = delta_P_i_dist + delta_P_i_loc
delta_P_o_dist = (4*G_o*G_o*Lhp*Nhp*fF_o)/(ro_o*Dh)
K_o = 2.8 +2000/Re_o +71.12/(Dh*1000)
delta_P_o_loc = (G_o*G_o*K_o*Nhp)/(2*ro_o)
delta_P_o = delta_P_o_loc+delta_P_o_dist

# Dicionário
data = {
    'Variável': ['T2 out', 'q1', 'q2', 'ΔT1', 'ΔT2', 'ΔTm', 'Di', 'Do', 'Ds','μi', 'Si', 'Gi', 'Rei',
                 'fFi', 'K1i', 'K2i', 'hi', 'μo', 'De', 'So', 'Go', 'Reo', 'fFo', 'K1o', 'K2o', 'ho',
                 'LfHdi', 'LfHdo', 'K', 'Uclean', 'Uc', 'Rdi', 'Rdo', 'Ud', 'Alinear', 'Aclean', 'Adirty', 'Ad',
                 'L', 'Lhp', 'Nhp', 'Areal', 'A', 'q', 'U', 'Rd', 'ρi', 'ρo', 'ΔPi dist', 'Ki', 'ΔPi loc',
                 'ΔPi', 'ΔPo dist', 'Ko', 'ΔPo loc', 'ΔPo'],
    '.':['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',
          '.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',
            '.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.' ],
    'Valor': [T2_out, q1, q2, delta_T1, delta_T2, delta_Tm, Di, Do, Ds, M_i,S_i, G_i, Re_i, fF_i, K1_i, K2_i, h_i, M_o, De,
              S_o, G_o, Re_o, fF_o, K1_o, K2_o, h_o, LfHd_i, LfHd_o, K, Uclean, Uc, Rd_i, Rd_o, Ud, Alinear, Aclean,
              Adirty, Ad, L, Lhp, Nhp, Areal, A, q, U, Rd, ro_i, ro_o, delta_P_i_dist, K_i, delta_P_i_loc, delta_P_i,
              delta_P_o_dist, K_o, delta_P_o_loc, delta_P_o],
     '..':['..','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',
          '.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',
            '.','.','.','.','.','.','.','.','.','.','.','.','.','.','.','.' ],

    'unidade': ["K","W","W",'K','K','K','M','M','M','Pa.s','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25',
                '26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51',
                '52','53','54','55','56']

}

# DataFrame
df = pd.DataFrame(data)
print(df)