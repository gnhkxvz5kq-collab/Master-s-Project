import Trovato_Purkinje_Model_corrected as tpf
import numpy as np
V_cell, A_geo, A_cap, V_myo, V_nsr, V_jsr, V_csr, V_ss, V_sl = tpf.Cell_Geometry()
print(V_cell)
'''
Steady state conditions at 1Hz (after 1000 beats) taken from Trovato model supplementary material
V = 86.55
[Na^+]_I = 8.23
[Na^+]_SL]= 8.23
[Na^+]_SS = 8.23
[K^+]_ I = 144
[K^+]_SL = 144
[K^+]_SS = 144
[Ca^2+]_I = 4.36e-05
[Ca^2+]_SL 1.0e-04
[Ca^2+]_SS] = 1e-04
[Ca^2+]_JSR = 1.25
[Ca^2+]_NSR = 1.27
[Ca^2+]_CSR 1.27
m = 0
hf = 0.8
hs = 0.8
j = 0.8
hsp = 0.6
jp = 0.8
mL = 0
hL = 0.5
hLp = 0.2
a = 0 
i = 0.6
i2 = 1
d = 0
ff = 1
fs = 1
fcaf = 1
fcas = 1
jca = 1
nca = 0 
ffp = 1
fcafp = 1
b = 0 
g = 1
xrf = 0 
xs1 = 0.2
xs2 = 0 
xk1 = 1
y = 0.2
CaMKt = 0 
u = 0.5
Jrel_1 = 0 
Jrel_2 = 0
'''

'''will run the function to return derivatives and use python functions for any integration necessary w.r.t time'''
