import Trovato_Purkinje_Model_corrected as tpf
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
V_cell, A_geo, A_cap, V_myo, V_nsr, V_jsr, V_csr, V_ss, V_sl = tpf.Cell_Geometry()

'''First need to arrange the Purkinje fibres into an array which contains details of their positions
Can extract length of cells from the cell geometry defined in tpf to separate centre of cells - L=0.0164
'''
length = 0.0164 * 4
positions = np.arange(0,length,0.0164) #array containing the positions of 4 Purkinje cells starting at 0 separated by Physical length 



