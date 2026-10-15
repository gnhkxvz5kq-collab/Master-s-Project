import numpy as np 
def stimulus_function():
    #I_stim is in micro amps per micro faraday
    #duration is in ms
    I_stim = -40
    duration = 1
    return I_stim

def extracelllular_concentrations_fun():
    #all measure in mM
    Nao = 140
    Cao = 1.8
    Ko = 5.4
    return Nao, Cao, Ko

def Cell_Geometry_func(L = 0.0164, r = 0.00175):
    A_geo = 2 * np.pi * r**2 + 2 * np.pi * r * L
    A_cap = 2 * A_geo 
    v_cell = np.pi * r**2 * L
    v_myo = 0.6 * v_vell
    v_ss = 0.02 * v_cell
    v_sl = 0.15 * v_cell
    v_nsr = 0.04 * v_cell
    v_jsr = 0.002 * v_cell
    v_csr = 0.008 * v_cell
    return A_geo, A_cap v_cell, v_myo, v_ss, v_sl, v_nsr, v_jsr, v_csr

def Votage_differential(I_ion, I_stim, C_m):
    dV_dt = -(I_ion + I_stim)/C_m
    return dV_dt

def fast_INa_func(V, h_fast, h_slow, h_inf, j):

    G_Na = 39.46 #max current conductance in ms / microF
    m_inf = 1/(1 + np.exp(-(V + 48.6264)/9.871))
    tau_m = 1 / (6.775 * np.exp((V+ 1.64)/7.8579) + 8.552 * np.exp(-(V + 77.42)/5.955))
    h_inf = 1 / (1 + exp((V + 78.5)/6.22))
    tau_hfast = 1/ (3.686 * 10e-6 * np.exp((V+ 3.8875)/7.8579) + 16 * np.exp((V-0.4963)/9.1843))
    tau_hslow = 1/ (0.009764 * np.exp(-(V+17.95)/28.05) + 0.3343*np.exp((V+5.730)/56.66))

    A_hfast = 0.99
    A_hslow = 1 - A_hfast
    h = A_hfast * h_fast + A_hslow * h_slow
    j_inf = h_inf
    tau_j = 4.8590 + 1/(0.8628 * np.exp(-(V+166.7)/7.6) + 1.1096 * np.exp((V+6.2719)/9.0358))
    h_CaMK_inf = 1 / (1 + np.exp((V+84.7)/6.22))
    tau_h_CaMK_slow = 3 * tau_hslow
    A_h_CaMK_fast = A_hfast
    A_h_CaMK_slow = A_hslow
    h_CaMK = A_h_CaMK_fast * h_fast + A_h_CaMK_slow * h_slow
    j_CaMk_inf = j_inf
    tau_j_CaMK = 1.46 * tau_j
    K_m_CaMK = 0.15
    phi_Na_CaMK = 1 / (1 + (K_m_CaMK)K_m_CaMK_active)
    I_Na_fast = G_Na * (V - E_Na) * m_inf ** 3 * ((1 - phi_Na_CaMK)* h*j + phi_Na_CaMK * h_CaMK * j_CaMk_inf)
