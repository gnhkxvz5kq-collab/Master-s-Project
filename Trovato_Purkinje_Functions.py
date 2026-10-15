import numpy as np 
#need to convert the cellml file for the Trovato model into a python version that can be used instead 
def initial_values():
    """
    Returns all initial state variables and model parameters
    required to start the Trovato Purkinje model.
    """

    return {

        # --------------------------------------------------
        # Membrane / stimulus
        # --------------------------------------------------
        "v": -86.6814002878592,
        "amp": -40.0,
        "duration": 1.0,

        # --------------------------------------------------
        # CaMK
        # --------------------------------------------------
        "CaMKt": 0.00505983330678751,

        # --------------------------------------------------
        # Sodium current gates (INa)
        # --------------------------------------------------
        "m": 0.00632661703915808,
        "hf": 0.788611739889677,
        "hs": 0.788545979951331,
        "j": 0.790474358603666,
        "hsp": 0.579693514309867,
        "jp": 0.790947058236417,

        # --------------------------------------------------
        # Late sodium current (INaL)
        # --------------------------------------------------
        "mL": 0.000241925773627233,
        "hL": 0.463574582508218,
        "hLp": 0.240216198686475,

        # --------------------------------------------------
        # Ito
        # --------------------------------------------------
        "a": 0.000272851144435704,
        "i1": 0.649604795721571,
        "i2": 0.989965695822495,

        # --------------------------------------------------
        # ICaL
        # --------------------------------------------------
        "d": 6.97735089296892e-9,
        "ff": 0.999999968230738,
        "fs": 0.926692153319136,
        "fcaf": 0.99999996819573,
        "fcas": 0.999999905741936,
        "jca": 0.999978907334662,
        "ffp": 0.999999968365903,
        "fcafp": 0.999999968278239,
        "nca": 0.00547252500964926,

        # --------------------------------------------------
        # ICaT
        # --------------------------------------------------
        "b": 0.000304250912559619,
        "g": 0.994214357917907,

        # --------------------------------------------------
        # IKr
        # --------------------------------------------------
        "xrf": 0.000331691184084272,
        "xrs": 0.568716473334161,

        # --------------------------------------------------
        # IKs
        # --------------------------------------------------
        "xs1": 0.191165248085394,
        "xs2": 0.000222677365291219,

        # --------------------------------------------------
        # If
        # --------------------------------------------------
        "y": 0.233119011214908,

        # --------------------------------------------------
        # IK1
        # --------------------------------------------------
        "xk1": 0.997084813729909,

        # --------------------------------------------------
        # Intracellular ion concentrations
        # --------------------------------------------------
        "nai": 8.23183964616932,
        "nasl": 8.23153516580562,
        "nass": 8.23154325237268,

        "ki": 143.767359809132,
        "ksl": 143.767769906216,
        "kss": 143.767768218104,

        "cai": 4.36004404734282e-5,
        "cass": 1.01777993438818e-4,
        "casl": 1.02004317781147e-4,

        "cansr": 1.26350902016858,
        "cajsr": 1.24811940209535,
        "cacsr": 1.26516959198518,

        # --------------------------------------------------
        # RyR release flux states
        # --------------------------------------------------
        "Jrel1": 1.08240945806962e-4,
        "Jrel2": 1.25045800437317e-69,

        # --------------------------------------------------
        # IP3 receptor state
        # --------------------------------------------------
        "u": 0.466236137183558
    }
#will set up initialisation values later - need to make a list of those that need initialising
#currently all initial values are established within each function component but will change this later
def extracellular():
    # concentrations of [Na^+]_o, [Ca^2+]_o, [K^+]_o in mM for extracellular ionic concentrations
    nao = 140.0
    cao = 1.8
    ko = 5.4

    return nao, cao, ko

def physical_constants():
    R = 8314.0
    T = 310.0
    F = 96485.0

    zna = 1.0
    zca = 2.0
    zk = 1.0

    return R, T, F, zna, zca, zk

def Cell_Geometry(L = 0.0164, rad=0.00175):
    #returns 9 values 
    #this function returns all cell geometry values needed later in the model
    #both L and rad are in cm 
    V_cell = 1000 * np.pi * rad**2 * L #cell volume in cm^3
    A_geo = 2 * np.pi * rad**2 + 2 * np.pi * rad * L #cell surface area in cm^2 
    A_cap = 2 * A_geo #find out what this represents
    V_myo = 0.6 * V_cell
    V_nsr = 0.04 * V_cell
    V_jsr = 0.002 * V_cell
    V_csr = 0.008 * V_cell
    V_ss = 0.02 * V_cell
    V_sl = 0.15 * V_cell
    return V_cell, A_geo, A_cap, V_myo, V_nsr, V_jsr, V_csr, V_ss, V_sl

def Membrane_Potential(R,T,F,C_m, I_Na, I_NaL, I_to, I_sus, I_CaL, I_CaT, I_CaNa, I_CaK, I_Kr, I_Ks, I_f, I_K1, I_NaCa_i, I_NaCa_ss, I_NaK, I_Nab, I_pCa, I_Cab, amp, v, duration):
    if time <= duration:
        I_stim = amp
    else:
        I_stim = 0
    vffrt = v * F**2 / (R*T)
    vfrt = v * F / (R*T)
    dv_dt = -(I_Na + I_NaL + I_to + I_sus + I_CaL + I_CaT + I_CaNa + I_CaK + I_Kr + I_Ks + I_f + I_K1 + I_NaCa_i + I_NaCa_ss + I_NaK + I_Nab + I_pCa + I_Cab + 2 * I_stim) /C_m
    #returns the derivative of the Voltage across the membrane with respect to time
    return dv_dt, vffrt, vfrt

def CaMK_func(cass, time, aCaMK, bCaMK, CaMKo, KmCaM, CaMKt ):
    
    
    CaMKb = CaMKo * (1 - CaMKt)/(1 + (KmCaM/cass))
    CaMKa = CaMKb + CaMKt

    dCaMKa_dt = aCaMK * CaMKb * (CaMKb + CaMKt) - bCaMK * CaMKt

    return CaMKb, CaMKa, dCaMKa_dt

def Reverse_Potentials(R,T,F, nao, ko, cao, nasl, ksl, casl, PKNa):
    
    ENa = (R*T/F) * np.log(nao/nasl)
    EK = (R*T/F) * np.log(ko/ksl)
    ECa = (0.5 * R*T/F) * np.log(cao/casl)
    EKs = (R*T/F) * np.log((ko + PKNa * nao)/(ksl + PKNa * nasl))

    return ENa, EK, ECa, EKs

def INa_func(KmCaMK, v, CaMKa, ENa, m, Ahf, hf, j, jp):
    #initial values
    # INa component: parameters and initial state values

    # Voltage-dependence parameters for activation gate m
    mssV1 = 48.4264
    mssV2 = 7.5653
    mtV1 = 11.64
    mtV2 = 34.77
    mtD1 = 6.765
    mtD2 = 8.552
    mtV3 = 77.42
    mtV4 = 5.955

    # Voltage-dependence parameters for inactivation gates
    hssV1 = 78.5
    hssV2 = 6.22

    # Maximum sodium conductance
    GNa = 39.4572


    mss = 1/ (1 + np.exp((-v-mssV1)/mssV2))
    tm = 1 / (mtD1 * np.exp((v + mtV1)/mtV2 + mtD2 * np.exp(-(v + mtV3)/mtV4)))
    dm_dt = (mss - m0)/tm
    hsss = 1 / (1 + np.exp((v + hssV1)/hssV2))
    thf = 1 / (3.686 * 10 **-6 * np.exp(-(v + 3.8875)/7.8579) + 16 * np.exp((v - 0.4963)/9.1843))
    ths = 1 / (0.-0.009794 * np.exp(-(v + 17.95)/28.05) + 0.3343 * np.exp((v+ 5.73)/56.66))
    Ahs = 1 - Ahf
    dhf_dt = (hss-hf)/thf
    dhs_dt = (hss - hs)/ths
    h = Ahf * hf + Ahs * hs
    jss = hss
    tj = 4.859 + 1/ (0.8628 * np.exp(-(v + 116.728)/7.6005) + 1.1096 * np.exp((v + 6.2719)/9.0358))
    dj_dt = (jss - j)/tj
    hssp = 1 / (1 + np.exp((v + 84.7)/6.22))

    thsp = 3 * ths
    dhsp_dt = (hssp - hsp)/thsp
    hp = Ahf * hf + Ahs * hsp
    tjp = 1.46 * tj 
    djp_dt = (jss - jp)/tjp
    fINaP = 1/(1 + KmCaMK/CaMKa)
    INa = GNa * (v - ENa) * m**3 * ((1 - fINaP) * h * j + fINaP * hp * jp)

    return INa, tm

def INaL_func(tm, v, KmCaMK, CaMKa, ENa, time):
    #initial values
    mL = 0.000241925773627233
    thL = 200
    hL = 0.463574582508218
    hLp = 0.240216198686475
    GNaL = 0.0189

    #calculations for current flow 
    mLss = 1/(1 + np.exp(-(v + 42.85)/5.264))
    dmL_dt = (mLss - mL)/ tm
    hLss = 1 / (1 + np.exp((v + 87.61)/7.488))
    dhL_dt = (hLss - hL)/thL
    hLssp = 1 / (1 + np.exp((v + 93.81)/7.488))
    dhLp_dt = (hLssp - hLp)/thLp
    fINaLP = 1/(1 + KmCaMK/CaMKa)

    INaL = GNaL * (v - ENa) * mL * ((1 - fINaLP) * hL + fINaLP * hLp)
    
    return INaL

def Ito_func(v, EK, a, i1, i2):

    # Fixed parameter
    Gto = 0.192  # mS/microF

    # Steady-state gating values
    ass = 1.0 / (1.0 + np.exp((20.0 - v) / 13.0))
    iss = 1.0 / (1.0 + np.exp((v + 27.0) / 13.0))

    # Voltage-dependent time constants (ms)
    taua = 1.0515 / (
        1.0 / (1.2089 * (1.0 + np.exp(-(v - 18.4099) / 29.3814)))
        + 3.5 / (1.0 + np.exp((v + 100.0) / 29.3814))
    )

    tauis = 43.0 + 1.0 / (
        0.001416 * np.exp(-(v + 96.52) / 59.05)
        + 0.0000000178 * np.exp((v + 114.1) / 8.079)
    )

    tauif = 6.162 + 1.0 / (
        0.3933 * np.exp(-(v + 100.0) / 100.0)
        + 0.08004 * np.exp((v - 8.0) / 8.59)
    )

    # Gating-variable derivatives
    da_dt = (ass - a) / taua
    di1_dt = (iss - i1) / tauis
    di2_dt = (iss - i2) / tauif

    # Transient outward potassium current
    Ito = Gto * a * i1 * i2 * (v - EK)

    return Ito, da_dt, di1_dt, di2_dt

def Isus_func(v, EK):
    # Fixed parameter
    Gsus = 0.0301  # mS/microF

    # Steady-state activation
    asus = 1.0 / (1.0 + np.exp(-(v - 12.0) / 16.0))

    # Sustained outward potassium current
    Isus = Gsus * asus * (v - EK)

    return Isus   

def ICaL_func(
    v, cass, cao, vffrt, vfrt,
    nao, nass, ko, kss,
    KmCaMK, CaMKa,
    d, ff, fs, fcaf, fcas, jca, ffp, fcafp, nca, Kmn, k2n,PCa
):
    """
    L-type calcium channel component from the supplied CellML.

    Voltage: mV
    Concentrations: mM
    Time: ms

    Returns:
        ICaL, ICaNa, ICaK,
        dd_dt, dff_dt, dfs_dt, dfcaf_dt, dfcas_dt,
        djca_dt, dffp_dt, dfcafp_dt, dnca_dt
    """

    Aff = 0.6
    Afs = 1.0 - Aff

    # --------------------------------------------------
    # 2. Steady-state activation and inactivation
    # --------------------------------------------------
    dss = 1.0 / (
        1.0 + np.exp(-(v + 3.94 + 3.3) / 4.23)
    )

    fss = 1.0 / (
        1.0 + np.exp((v + 19.58 + 3.3) / 3.696)
    )

    fcass = fss

    # --------------------------------------------------
    # 3. Voltage-dependent time constants (ms)
    # --------------------------------------------------
    td = 0.6 + 1.0 / (
        np.exp(-0.05 * (v + 6.0))
        + np.exp(0.09 * (v + 14.0))
    )

    tff = 7.0 + 1.0 / (
        0.0045 * np.exp(-(v + 20.0 + 15.19) / 10.0)
        + 0.0045 * np.exp((v + 20.0 + 15.19) / 10.0)
    )

    tfs = 1000.0 + 1.0 / (
        0.000035 * np.exp(-(v + 5.0 + 15.19) / 4.0)
        + 0.000035 * np.exp((v + 5.0 + 15.19) / 6.0)
    )

    tfcaf = 0.72 * (
        7.0 + 1.0 / (
            0.04 * np.exp(-((v + 15.19) - 4.0) / 7.0)
            + 0.04 * np.exp(((v + 15.19) - 4.0) / 7.0)
        )
    )

    tfcas = 0.49 * (
        100.0 + 1.0 / (
            0.00012 * np.exp(-(v + 15.19) / 3.0)
            + 0.00012 * np.exp((v + 15.19) / 7.0)
        )
    )

    tjca = 75.0

    # --------------------------------------------------
    # 4. Fast and slow inactivation weights
    # --------------------------------------------------
    Afcaf = 0.3 + 0.6 / (
        1.0 + np.exp((v - 10.0) / 10.0)
    )

    Afcas = 1.0 - Afcaf

    # --------------------------------------------------
    # 5. Derivatives of gating state variables
    # --------------------------------------------------
    dd_dt = (dss - d) / td

    dff_dt = (fss - ff) / tff
    dfs_dt = (fss - fs) / tfs

    dfcaf_dt = (fcass - fcaf) / tfcaf
    dfcas_dt = (fcass - fcas) / tfcas

    djca_dt = (fcass - jca) / tjca

    # Phosphorylated gates
    tffp = 2.5 * tff
    dffp_dt = (fss - ffp) / tffp

    tfcafp = 2.5 * tfcaf
    dfcafp_dt = (fcass - fcafp) / tfcafp

    # --------------------------------------------------
    # 6. Combined inactivation gates
    # --------------------------------------------------
    f = Aff * ff + Afs * fs
    fp = Aff * ffp + Afs * fs

    fca = Afcaf * fcaf + Afcas * fcas
    fcap = Afcaf * fcafp + Afcas * fcas

    # --------------------------------------------------
    # 7. Calcium-dependent inactivation state
    # --------------------------------------------------
    km2n = jca

    anca = 1.0 / (
        k2n / km2n + (1.0 + Kmn / cass)**4
    )

    dnca_dt = anca * k2n - nca * km2n

    # --------------------------------------------------
    # 8. Calcium, sodium and potassium driving terms
    # --------------------------------------------------
    PhiCaL = (
        4.0 * vffrt
        * (cass * np.exp(2.0 * vfrt) - 0.341 * cao)
        / (np.exp(2.0 * vfrt) - 1.0)
    )

    PhiCaNa = (
        vffrt
        * (0.75 * nass * np.exp(vfrt) - 0.75 * nao)
        / (np.exp(vfrt) - 1.0)
    )

    PhiCaK = (
        vffrt
        * (0.75 * kss * np.exp(vfrt) - 0.75 * ko)
        / (np.exp(vfrt) - 1.0)
    )

    # --------------------------------------------------
    # 9. Relative permeabilities
    # --------------------------------------------------
    PCap = 1.1 * PCa

    PCaNa = 0.00125 * PCa
    PCaK = 3.574e-4 * PCa

    PCaNap = 0.00125 * PCap
    PCaKp = 3.574e-4 * PCap

    # Fraction of channels affected by CaMK phosphorylation
    fICaLp = 1.0 / (1.0 + KmCaMK / CaMKa)

    # --------------------------------------------------
    # 10. Shared gate terms for phosphorylated/unphosphorylated
    # --------------------------------------------------
    gate = f * (1.0 - nca) + jca * fca * nca
    gate_p = fp * (1.0 - nca) + jca * fcap * nca

    # --------------------------------------------------
    # 11. L-type calcium current
    # --------------------------------------------------
    ICaL = (
        (1.0 - fICaLp) * PCa * PhiCaL * d * gate
        + fICaLp * PCap * PhiCaL * d * gate_p
    )

    # --------------------------------------------------
    # 12. Sodium current through the L-type channel
    # --------------------------------------------------
    ICaNa = (
        (1.0 - fICaLp) * PCaNa * PhiCaNa * d * gate
        + fICaLp * PCaNap * PhiCaNa * d * gate_p
    )

    # --------------------------------------------------
    # 13. Potassium current through the L-type channel
    # --------------------------------------------------
    ICaK = (
        (1.0 - fICaLp) * PCaK * PhiCaK * d * gate
        + fICaLp * PCaKp * PhiCaK * d * gate_p
    )

    return (
        ICaL, ICaNa, ICaK,
        dd_dt, dff_dt, dfs_dt,
        dfcaf_dt, dfcas_dt, djca_dt,
        dffp_dt, dfcafp_dt, dnca_dt
    )

def ICaT_func(v, ECa, b, g):
    # Fixed model parameter
    GCaT = 0.0754

    # Steady-state values of the gating variables

    bss = 1 / (1 + np.exp(-(v + 30) / 7))
    gss = 1 / (1 + np.exp((v + 61) / 5))

    # Voltage-dependent time constants (ms)
    taub = 1 / (
        1.068 * np.exp((v + 16.3) / 30)
        + 1.068 * np.exp(-(v + 16.3) / 30)
    )

    taug = 1 / (
        0.015 * np.exp((v + 71.7) / 15.4)
        + 0.015 * np.exp(-(v + 71.7) / 83.3)
    )

    # ODEs for the gating variables
    db_dt = (bss - b) / taub
    dg_dt = (gss - g) / taug

    # T-type calcium current
    ICaT = GCaT * b * g * (v - ECa)

    return ICaT, db_dt, dg_dt

def IKr_func(v, EK, ko, xrf, xrs):
    # Fixed parameter
    GKr = 0.0342

    # Steady-state activation
    xrss = 1 / (1 + np.exp(-(v + 8.337) / 6.789))

    # Fast and slow activation time constants (ms)
    txrf = 12.98 + 1 / (
        0.3652 * np.exp((v + 17.6 - 31.66) / 3.869)
        + 4.123e-5 * np.exp(-((v + 17.6) - 47.78) / 20.38)
    )

    txrs = 1.865 + 1 / (
        0.06629 * np.exp((v + 17.2 - 34.7) / 7.355)
        + 1.128e-5 * np.exp(-((v + 17.2) - 29.74) / 25.94)
    )

    # Weighting factors for fast and slow gates
    Axrf = 1 / (1 + np.exp((v + 54.81) / 38.21))
    Axrs = 1 - Axrf

    # Gating-variable derivatives
    dxrf_dt = (xrss - xrf) / txrf
    dxrs_dt = (xrss - xrs) / txrs

    # Combined activation gate
    xr = Axrf * xrf + Axrs * xrs

    # Rectification factor
    rkr = (
        1 / (1 + np.exp((v + 55) / (0.32 * 75)))
    ) / (
        1 + np.exp((v - 10) / (0.32 * 30))
    )

    # Rapid delayed rectifier potassium current
    IKr = GKr * np.sqrt(ko / 5.4) * xr * rkr * (v - EK)

    return IKr, dxrf_dt, dxrs_dt

def IKs_func(v, EKs, casl, xs1, xs2):
    # Fixed parameter
    GKs = 0.0029

    # Steady-state activation
    xs1ss = 1 / (1 + np.exp(-(v + 11.6) / 8.932))
    xs2ss = xs1ss

    # Time constants (ms)
    txs1 = 817.3 + 1 / (
        2.326e-4 * np.exp((v + 48.28) / 17.8)
        + 0.001292 * np.exp(-(v + 210) / 230)
    )

    txs2 = 1 / (
        0.01 * np.exp((v - 50) / 20)
        + 0.0193 * np.exp(-(v + 66.54) / 31)
    )

    # Calcium-dependent scaling factor
    KsCa = 1 + 0.6 / (
        1 + (3.8e-5 / casl)**1.4
    )

    # Gating-variable derivatives
    dxs1_dt = (xs1ss - xs1) / txs1
    dxs2_dt = (xs2ss - xs2) / txs2

    # Slow delayed rectifier potassium current
    IKs = GKs * KsCa * xs1 * xs2 * (v - EKs)

    return IKs, dxs1_dt, dxs2_dt

def If_func(v, ENa, EK, y):
    # Fixed conductances (mS/µF)
    GfNa = 0.0116
    GfK = 0.0232

    # Steady-state activation and time constant
    yss = 1 / (1 + np.exp((v + 87) / 9.5))

    tauy = 2000 / (
        np.exp((v + 57) / 60)
        + np.exp(-(v + 132) / 10)
    )

    # Gating-variable derivative
    dy_dt = (yss - y) / tauy

    # Sodium and potassium components
    IfNa = GfNa * y**2 * (v - ENa)
    IfK = GfK * y**2 * (v - EK)

    # Total funny current
    If = IfNa + IfK

    return If, IfNa, IfK, dy_dt

def IK1_func(v, EK, ko, xk1):
    # Fixed conductance (mS/µF)
    GK1 = 0.0455

    # Steady-state gate
    xk1ss = 1 / (
        1 + np.exp(
            -(v + 2.5538 * ko + 144.59)
            / (1.5692 * ko + 3.8115)
        )
    )

    # Time constant (ms)
    txk1 = 122.2 / (
        np.exp(-(v + 127.2) / 20.36)
        + np.exp((v + 236.8) / 69.33)
    )

    # Rectification factor
    rk1 = 1 / (
        1 + np.exp((v + 116 - 5.5 * ko) / 11)
    )

    # Gate derivative
    dxk1_dt = (xk1ss - xk1) / txk1

    # Inward rectifier potassium current
    IK1 = (
        GK1
        * 2.3238
        * np.sqrt(ko / 5.4)
        * rk1
        * xk1
        * (v - EK)
    )

    return IK1, dxk1_dt

def _ncx_fluxes(na_i, ca_i, nao, cao, v, R, T, F):
    """
    Shared algebraic calculations for the Na+/Ca2+ exchanger.
    na_i and ca_i are the local intracellular concentrations.
    """

    # Fixed kinetic parameters
    kna1 = 15.0
    kna2 = 5.0
    kna3 = 88.12
    kasymm = 12.5

    wna = 60000.0
    wca = 60000.0
    wnaca = 5000.0

    kcaon = 1500000.0
    kcaoff = 5000.0

    qna = 0.5224
    qca = 0.167
    KmCaAct = 0.00015

    # Voltage-dependent factors
    hna = np.exp(qna * v * F / (R * T))
    hca = np.exp(qca * v * F / (R * T))

    # Sodium-binding intermediates, intracellular side
    h1 = 1 + (na_i / kna3) * (1 + hna)
    h2 = na_i * hna / (kna3 * h1)
    h3 = 1 / h1

    h4 = 1 + (na_i / kna1) * (1 + na_i / kna2)
    h5 = na_i**2 / (h4 * kna1 * kna2)
    h6 = 1 / h4

    # Sodium-binding intermediates, extracellular side
    h7 = 1 + (nao / kna3) * (1 + 1 / hna)
    h8 = nao / (kna3 * hna * h7)
    h9 = 1 / h7

    h10 = (
        kasymm + 1
        + (nao / kna1) * (1 + nao / kna2)
    )
    h11 = nao**2 / (h10 * kna1 * kna2)
    h12 = 1 / h10

    # Transition rates
    k1 = h12 * cao * kcaon
    k2 = kcaoff

    k3p = h9 * wca
    k3pp = h8 * wnaca
    k3 = k3p + k3pp

    k4p = h3 * wca / hca
    k4pp = h2 * wnaca
    k4 = k4p + k4pp

    k5 = kcaoff
    k6 = h6 * ca_i * kcaon
    k7 = h5 * h2 * wna
    k8 = h8 * h11 * wna

    # Occupancy intermediates
    x1 = (
        k2 * k4 * (k7 + k6)
        + k5 * k7 * (k2 + k3)
    )

    x2 = (
        k1 * k7 * (k4 + k5)
        + k4 * k6 * (k1 + k8)
    )

    x3 = (
        k1 * k3 * (k7 + k6)
        + k8 * k6 * (k2 + k3)
    )

    x4 = (
        k2 * k8 * (k4 + k5)
        + k3 * k5 * (k1 + k8)
    )

    denom = x1 + x2 + x3 + x4

    E1 = x1 / denom
    E2 = x2 / denom
    E3 = x3 / denom
    E4 = x4 / denom

    # Allosteric calcium activation
    allo = 1 / (1 + (KmCaAct / ca_i)**2)

    # Sodium and calcium exchanger fluxes
    JncxNa = (
        3 * (E4 * k7 - E1 * k8)
        + E3 * k4pp
        - E2 * k3pp
    )

    JncxCa = E2 * k2 - E1 * k1

    return JncxNa, JncxCa, allo

def INaCa_i_func(
    v, R, T, F,
    nass, cass, nao, cao,
    nasl, casl, zna, zca
):
    # Fixed exchanger conductance
    Gncx = 0.00095709

    # Intracellular/local-myoplasmic exchanger contribution
    Jna_i, Jca_i, allo_i = _ncx_fluxes(
        na_i=nasl,
        ca_i=casl,
        nao=nao,
        cao=cao,
        v=v,
        R=R,
        T=T,
        F=F
    )

    INaCa_i = (
        0.8 * Gncx * allo_i
        * (zna * Jna_i + zca * Jca_i)
    )

    # Subspace exchanger contribution
    Jna_ss, Jca_ss, allo_ss = _ncx_fluxes(
        na_i=nass,
        ca_i=cass,
        nao=nao,
        cao=cao,
        v=v,
        R=R,
        T=T,
        F=F
    )

    INaCa_ss = (
        0.2 * Gncx * allo_ss
        * (zna * Jna_ss + zca * Jca_ss)
    )

    return INaCa_i, INaCa_ss

def INaK_func(v, R, T, F, nasl, ksl, nao, ko, zna, zk):
    """
    Na+/K+ ATPase pump current.

    Parameters
    ----------
    v : float
        Membrane voltage (mV).
    R, T, F : float
        Gas constant, absolute temperature, and Faraday constant,
        using the same conventions as the CellML model.
    nasl, ksl : float
        Subsarcolemmal Na+ and K+ concentrations (mM).
    nao, ko : float
        Extracellular Na+ and K+ concentrations (mM).
    zna, zk : float
        Na+ and K+ valences.

    Returns
    -------
    INaK : float
        Na+/K+ pump current (microA/microF).
    """

    # Fixed parameters
    k1p = 949.5
    k1m = 182.4
    k2p = 687.2
    k2m = 39.4
    k3p = 1899.0
    k3m = 79300.0
    k4p = 639.0
    k4m = 40.0

    Knai0 = 9.073
    Knao0 = 27.78
    delta = -0.155

    Kki = 0.5
    Kko = 0.3582

    MgADP = 0.05
    MgATP = 9.8
    Kmgatp = 1.698e-7
    H = 1e-7
    eP = 4.2
    Khp = 1.698e-7
    Knap = 224.0
    Kxkur = 292.0

    Pnak = 32.4872

    # Voltage-dependent Na+ affinities
    Knai = Knai0 * np.exp(delta * v * F / (3.0 * R * T))
    Knao = Knao0 * np.exp((1.0 - delta) * v * F / (3.0 * R * T))

    # Phosphorylation factor
    P = eP / (
        1.0
        + H / Khp
        + nasl / Knap
        + ksl / Kxkur
    )

    # Transition rates a1...a4 and b1...b4
    a1 = (
        k1p * (nasl / Knai)**3
        / (
            (1.0 + nasl / Knai)**3
            + (1.0 + ksl / Kki)**2
            - 1.0
        )
    )

    b1 = k1m * MgADP

    a2 = k2p

    b2 = (
        k2m * (nao / Knao)**3
        / (
            (1.0 + nao / Knao)**3
            + (1.0 + ko / Kko)**2
            - 1.0
        )
    )

    a3 = (
        k3p * (ko / Kko)**2
        / (
            (1.0 + nao / Knao)**3
            + (1.0 + ko / Kko)**2
            - 1.0
        )
    )

    b3 = k3m * P * H / (1.0 + MgATP / Kmgatp)

    a4 = (
        k4p * (MgATP / Kmgatp)
        / (1.0 + MgATP / Kmgatp)
    )

    b4 = (
        k4m * (ksl / Kki)**2
        / (
            (1.0 + nasl / Knai)**3
            + (1.0 + ksl / Kki)**2
            - 1.0
        )
    )

    # Occupancy-state weighting factors
    x1 = (
        a4 * a1 * a2
        + b2 * b4 * b3
        + a2 * b4 * b3
        + b3 * a1 * a2
    )

    x2 = (
        b2 * b1 * b4
        + a1 * a2 * a3
        + a3 * b1 * b4
        + a2 * a3 * b4
    )

    x3 = (
        a2 * a3 * a4
        + b3 * b2 * b1
        + b2 * b1 * a4
        + a3 * a4 * b1
    )

    x4 = (
        b4 * b3 * b2
        + a3 * a4 * a1
        + b2 * a4 * a1
        + b3 * b2 * a1
    )

    x_sum = x1 + x2 + x3 + x4

    E1 = x1 / x_sum
    E2 = x2 / x_sum
    E3 = x3 / x_sum
    E4 = x4 / x_sum

    # Na+ and K+ transport fluxes
    JnakNa = 3.0 * (E1 * a3 - E2 * b3)
    JnakK = 2.0 * (E4 * b1 - E3 * a1)

    # Total pump current
    INaK = Pnak * (zna * JnakNa + zk * JnakK)

    return INaK

def INab_func(nasl, nao, vffrt, vfrt):
    """
    Background Na+ current.

    Returns
    -------
    INab : float
        Current (microA/microF).
    """

    PNab = 9.375e-10

    exp_vfrt = np.exp(vfrt)

    INab = (
        PNab * vffrt
        * (nasl * exp_vfrt - nao)
        / np.expm1(vfrt)
    )

    return INab

def ICab_func(casl, cao, vffrt, vfrt):
    """
    Background Ca2+ current.

    Returns
    -------
    ICab : float
        Current (microA/microF).
    """

    PCab = 2.5e-8

    exp_2vfrt = np.exp(2.0 * vfrt)

    ICab = (
        PCab * 4.0 * vffrt
        * (casl * exp_2vfrt - 0.341 * cao)
        / np.expm1(2.0 * vfrt)
    )

    return ICab

def IpCa_func(casl):
    """
    Sarcolemmal Ca2+ pump current.

    Returns
    -------
    IpCa : float
        Current (microA/microF).
    """

    GpCa = 0.0005
    KmCap = 0.0005

    IpCa = GpCa * casl / (KmCap + casl)

    return IpCa

def intracellular_ion_derivatives(
    X,
    currents,
    fluxes,
    volumes,
    Acap,
    F,
    cm=1.0
):
    """
    Calculate derivatives of intracellular ion concentrations.

    State order in X:
        nai, nasl, nass, ki, ksl, kss,
        cai, cass, casl, cansr, cajsr, cacsr

    currents: dictionary containing the currents used below.
    fluxes:   dictionary containing the ion and calcium fluxes.
    volumes:  dictionary containing compartment volumes.

    Returns
    -------
    dXdt : numpy array
        Derivatives in the same order as X, in mM/ms.
    """

    # Unpack state variables
    (
        nai, nasl, nass,
        ki, ksl, kss,
        cai, cass, casl,
        cansr, cajsr, cacsr
    ) = X

    # Fixed calcium-buffer parameters
    BSRmax = 0.019975
    KmBSR = 0.00087

    BSLmax = 0.4777
    KmBSL = 0.0087

    csqnmax = 2.88
    kmcsqn = 0.8
    csqnmaxsl = 1.2

    cmdnmax = 0.1125
    kmcmdn = 0.00238
    cmdnmaxsl = 0.0125

    trpnmax = 0.0315
    kmtrpn = 0.0005
    trpnmaxsl = 0.0035

    # Compartment volumes
    vmyo = volumes["vmyo"]
    vsl = volumes["vsl"]
    vss = volumes["vss"]
    vnsr = volumes["vnsr"]
    vjsr = volumes["vjsr"]
    vcsr = volumes["vcsr"]

    # Ionic currents
    INa = currents["INa"]
    INaL = currents["INaL"]
    INaCa_i = currents["INaCa_i"]
    INaCa_ss = currents["INaCa_ss"]
    INaK = currents["INaK"]
    INab = currents["INab"]
    IfNa = currents["IfNa"]
    IfK = currents["IfK"]
    IK1 = currents["IK1"]
    IKr = currents["IKr"]
    IKs = currents["IKs"]
    Istim = currents["Istim"]
    Ito = currents["Ito"]
    Isus = currents["Isus"]
    ICab = currents["ICab"]
    IpCa = currents["IpCa"]
    ICaL = currents["ICaL"]
    ICaT = currents["ICaT"]
    ICaNa = currents["ICaNa"]
    ICaK = currents["ICaK"]

    # Ionic and calcium fluxes
    JdiffNa = fluxes["JdiffNa"]
    JgapNa = fluxes["JgapNa"]

    JdiffK = fluxes["JdiffK"]
    JgapK = fluxes["JgapK"]

    Jdiff = fluxes["Jdiff"]
    Jgap = fluxes["Jgap"]

    Jup1 = fluxes["Jup1"]
    Jup2 = fluxes["Jup2"]

    Jrel1 = fluxes["Jrel1"]
    Jrel2 = fluxes["Jrel2"]
    Jip3 = fluxes["Jip3"]

    Jtr1 = fluxes["Jtr1"]
    Jtr2 = fluxes["Jtr2"]

    # --------------------------------------------------
    # Sodium concentration derivatives
    # --------------------------------------------------

    dnai_dt = JgapNa * vsl / vmyo

    dnasl_dt = (
        -(INa + INaL + 3 * INaCa_i + 3 * INaK + IfNa + INab)
        * Acap * cm / (F * vsl)
        + JdiffNa * vss / vsl
        - JgapNa
    )

    dnass_dt = (
        -(ICaNa + 3 * INaCa_ss)
        * cm * Acap / (F * vss)
        - JdiffNa
    )

    # --------------------------------------------------
    # Potassium concentration derivatives
    # --------------------------------------------------

    dki_dt = JgapK * vsl / vmyo

    dksl_dt = (
        -(Ito + Isus + IKr + IKs + IfK + IK1 + Istim - 2 * INaK)
        * cm * Acap / (F * vsl)
        + JdiffK * vss / vsl
        - JgapK
    )

    dkss_dt = (
        -ICaK * cm * Acap / (F * vss)
        - JdiffK
    )

    # --------------------------------------------------
    # Calcium buffering factors
    # --------------------------------------------------

    Bcai = 1 / (
        1
        + cmdnmax * kmcmdn / (kmcmdn + cai)**2
        + trpnmax * kmtrpn / (kmtrpn + cai)**2
    )

    Bcass = 1 / (
        1
        + BSRmax * KmBSR / (KmBSR + cass)**2
        + BSLmax * KmBSL / (KmBSL + cass)**2
    )

    Bcasl = 1 / (
        1
        + cmdnmaxsl * kmcmdn / (kmcmdn + casl)**2
        + trpnmaxsl * kmtrpn / (kmtrpn + casl)**2
    )

    Bcajsr = 1 / (
        1
        + csqnmaxsl * kmcsqn / (kmcsqn + cajsr)**2
    )

    Bcacsr = 1 / (
        1
        + csqnmax * kmcsqn / (kmcsqn + cacsr)**2
    )

    # --------------------------------------------------
    # Calcium concentration derivatives
    # --------------------------------------------------

    dcai_dt = Bcai * (
        Jgap * vsl / vmyo
        + Jrel2 * vcsr / vmyo
        - Jup2 * vnsr / vmyo
    )

    dcass_dt = Bcass * (
        -(ICaL - 2 * INaCa_ss) * cm * Acap / (2 * F * vss)
        + (Jrel1 + Jip3) * vjsr / vss
        - Jdiff
    )

    # Literal algebraic structure of the supplied CellML equation
    # for casl: its nested minus signs simplify to a positive
    # current term, as written here.
    dcasl_dt = Bcasl * (
        (IpCa + ICab + ICaT - 2 * INaCa_i)
        * cm * Acap / (2 * F * vsl)
        - Jup1 * vnsr / vsl
        - Jgap
        + Jdiff * vss / vsl
    )

    dcansr_dt = (
        Jup1 + Jup2
        - Jtr1 * vjsr / vnsr
        - Jtr2 * vcsr / vnsr
    )

    dcajsr_dt = Bcajsr * (
        Jtr1 - Jrel1 - Jip3
    )

    dcacsr_dt = Bcacsr * (
        Jtr2 - Jrel2
    )

    # Return derivatives in exactly the same order as X
    dXdt = np.array([
        dnai_dt,
        dnasl_dt,
        dnass_dt,
        dki_dt,
        dksl_dt,
        dkss_dt,
        dcai_dt,
        dcass_dt,
        dcasl_dt,
        dcansr_dt,
        dcajsr_dt,
        dcacsr_dt
    ])

    return dXdt


# ============================================================
# 1. DIFFUSION BETWEEN INTRACELLULAR COMPARTMENTS
# ============================================================

# Time constants (ms)
SSTAU = 0.2
GAPTAU = 12.0


def diffusion_fluxes(
    cai, casl, cass,
    nai, nasl, nass,
    ki, ksl, kss
):
    """
    Calculate sodium, potassium and calcium diffusion fluxes.

    Concentrations:
        cai, casl, cass : calcium concentrations (mM)
        nai, nasl, nass : sodium concentrations (mM)
        ki, ksl, kss    : potassium concentrations (mM)

    Compartment labels:
        i   : bulk intracellular/myoplasmic compartment
        sl  : subsarcolemmal compartment
        ss  : subspace compartment

    Returns:
        Dictionary of six fluxes, each in mM/ms:
            JdiffNa, JgapNa, JdiffK, JgapK, Jdiff, Jgap

    Jdiff* describes exchange between the subspace and
    subsarcolemmal compartments.
    Jgap* describes exchange between the subsarcolemmal
    and bulk intracellular compartments.
    """
    SSTAU = 0.2
    GAPTAU = 12.0

    # Sodium
    JdiffNa = (nass - nasl) / SSTAU
    JgapNa = (nasl - nai) / GAPTAU

    # Potassium
    JdiffK = (kss - ksl) / SSTAU
    JgapK = (ksl - ki) / GAPTAU

    # Calcium
    Jdiff = (cass - casl) / SSTAU
    Jgap = (casl - cai) / GAPTAU

    return {
        "JdiffNa": JdiffNa,
        "JgapNa": JgapNa,
        "JdiffK": JdiffK,
        "JgapK": JgapK,
        "Jdiff": Jdiff,
        "Jgap": Jgap
    }


# ============================================================
# 2. RYANODINE RECEPTOR (RyR) CALCIUM RELEASE
# ============================================================

# Initial values of the dynamic SR release fluxes
JREL1_0 = 0.000108240945806962
JREL2_0 = 1.25045800437317e-69


def ryr_derivatives(
    Jrel1, Jrel2,
    ICaL, Acap, F,
    vss, vsl, vcsr, vnsr, vjsr, vmyo,
    cajsr, cacsr,
    Jdiff, Jgap, Jip3, Jup2,
    CaMKa
):
    """
    Calculate RyR-related release fluxes and their derivatives.

    Dynamic state inputs:
        Jrel1 : calcium release flux from the junctional SR (mM/ms)
        Jrel2 : calcium release flux from the corbular SR (mM/ms)

    Other inputs:
        ICaL  : L-type calcium current (microA/microF)
        Acap  : membrane capacitance area factor (cm^2)
        F     : Faraday constant (C/mol)
        vss, vsl, vcsr, vnsr, vjsr, vmyo : compartment volumes
        cajsr, cacsr : calcium concentrations in junctional
                       and corbular SR (mM)
        Jdiff, Jgap  : calcium diffusion fluxes (mM/ms)
        Jip3         : IP3-mediated calcium flux (mM/ms)
        Jup2         : calcium uptake flux into SR (mM/ms)
        CaMKa        : activated CaMK concentration, as defined
                       in the source model (mM)

    Returns:
        REL, REL2     : driving terms used to calculate target
                        release fluxes (source-defined units)
        irelss, irelss2 : target release fluxes (mM/ms)
        ireltau, ireltau2 : release time constants as defined
                            in the source
        dJrel1_dt, dJrel2_dt : time derivatives of the two
                               release-flux state variables

    The derivatives are intended to be included in the ODE
    solver's state vector.
    """

    # Driving term for junctional SR release
    REL = -(
        (ICaL * Acap) / (vss * F * 2)
        - ((Jrel1 + Jip3) * vjsr / vss)
        + Jdiff
    )

    # CaMK-dependent factor used in both release pathways
    camk_factor = 1 + 1 / (1 + (0.28 / CaMKa)**8)

    # Release time constant for junctional SR
    ireltau = (
        2 * camk_factor
        / (1 + 0.0123 / cajsr)
    )

    # Target release flux for junctional SR.
    # The source sets the target to zero when REL <= 0.
    if REL > 0:
        irelss = (
            15 * camk_factor * REL
            / (1 + (1 / cajsr)**8)
        )
    else:
        irelss = 0.0

    # Time evolution of junctional SR release flux
    dJrel1_dt = (irelss - Jrel1) / ireltau

    # Driving term for corbular SR release
    REL2 = (
        Jgap * vsl / vmyo
        - Jup2 * vnsr / vmyo
        + Jrel2 * vcsr / vmyo
    )

    # Release time constant for corbular SR
    ireltau2 = (
        6 * camk_factor
        / (1 + 0.0123 / cacsr)
    )

    # Target release flux for corbular SR
    if REL2 > 0:
        irelss2 = (
            91 * camk_factor * REL2
            / (1 + (1 / cacsr)**8)
        )
    else:
        irelss2 = 0.0

    # Time evolution of corbular SR release flux
    dJrel2_dt = (irelss2 - Jrel2) / ireltau2

    return {
        "REL": REL,
        "REL2": REL2,
        "irelss": irelss,
        "irelss2": irelss2,
        "ireltau": ireltau,
        "ireltau2": ireltau2,
        "dJrel1_dt": dJrel1_dt,
        "dJrel2_dt": dJrel2_dt
    }


# ============================================================
# 3. SERCA CALCIUM UPTAKE
# ============================================================

# Parameters from the CellML source
DQUPCAMKBAR = 0.75
DKMPLBBAR = 0.00017
KMUP = 0.00028
NSRBAR = 15.0


def serca_fluxes(
    KmCaMK, CaMKa,
    cai, casl, cansr
):
    """
    Calculate SERCA-mediated calcium uptake into the SR.

    Inputs:
        KmCaMK : CaMK activation constant (mM)
        CaMKa  : activated CaMK concentration (mM)
        cai    : bulk intracellular calcium concentration (mM)
        casl   : subsarcolemmal calcium concentration (mM)
        cansr  : network SR calcium concentration (mM)

    Returns:
        dkmplb, dqupcamk : CaMK-dependent regulatory factors
        Jup1 : uptake flux associated with the subsarcolemmal
               calcium compartment (mM/ms)
        Jup2 : uptake flux associated with bulk intracellular
               calcium (mM/ms)

    The flux equations are transcribed from the supplied
    CellML source, including their original signs.
    """
    DQUPCAMKBAR = 0.75
    DKMPLBBAR = 0.00017
    KMUP = 0.00028
    NSRBAR = 15.0

    # CaMK-dependent regulation of phospholamban
    dkmplb = (
        DKMPLBBAR * CaMKa
        / (CaMKa + KmCaMK)
    )

    # CaMK-dependent regulation of SERCA uptake
    dqupcamk = (
        DQUPCAMKBAR * CaMKa
        / (CaMKa + KmCaMK)
    )

    # Calcium uptake associated with the subsarcolemmal pool
    Jup1 = (
        0.0002 * (dqupcamk + 1)
        / (1 + (KMUP - dkmplb) / casl)
        - 0.00105 * cansr / NSRBAR
    )

    # Calcium uptake associated with the bulk intracellular pool
    Jup2 = (
        0.0026 * (dqupcamk + 1)
        / (1 + (KMUP - dkmplb) / cai)
        - 0.0042 * cansr / NSRBAR
    )

    return {
        "dkmplb": dkmplb,
        "dqupcamk": dqupcamk,
        "Jup1": Jup1,
        "Jup2": Jup2
    }


# ============================================================
# 4. IP3-MEDIATED CALCIUM RELEASE
# ============================================================

# Fixed parameters from the CellML source
IP3_CONC = 0.0001
IP3_K1 = 150000.0
IP3_K1A = 16.5
IP3_K0 = 96000.0
IP3_K0A = 9.6
IP3_K2 = 1800.0
IP3_K2A = 0.21
TAUIP3 = 3.7

# Initial value of the dynamic IP3 activation variable
U_0 = 0.466236137183558


def ip3_derivatives(cass, cajsr, u):
    """
    Calculate IP3 activation, channel open probability, and
    IP3-mediated calcium release.

    Inputs:
        cass : subspace calcium concentration (mM)
        cajsr: junctional SR calcium concentration (mM)
        u    : dynamic IP3 activation state variable

    Returns:
        du_dt : time derivative of u
        POip3 : IP3-channel open probability as defined in
                the source equations
        Jip3  : IP3-mediated calcium flux (mM/ms)

    The parameter IP3_CONC is fixed at the source's initial
    value; it is not a dynamic state in this component.
    """
    IP3_CONC = 0.0001
    IP3_K1 = 150000.0
    IP3_K1A = 16.5
    IP3_K0 = 96000.0
    IP3_K0A = 9.6
    IP3_K2 = 1800.0
    IP3_K2A = 0.21
    TAUIP3 = 3.7

    # Evolution of the IP3 activation state
    du_dt = (
        cass * IP3_K2 * (1 - u)
        - IP3_K2A * u
    )

    # IP3-channel open probability
    POip3 = (
        TAUIP3 * IP3_CONC * cass * (1 - u)
        / (
            (1 + IP3_CONC * IP3_K0 / IP3_K0A)
            * (1 + cass * IP3_K1 / IP3_K1A)
        )
    )

    # Calcium flux through IP3-mediated release
    Jip3 = 10.92 * POip3 * (cajsr - cass)

    return {
        "du_dt": du_dt,
        "POip3": POip3,
        "Jip3": Jip3
    }


# ============================================================
# 5. TRANSFER BETWEEN SR COMPARTMENTS
# ============================================================

def trans_fluxes(cansr, cajsr, cacsr):
    """
    Calculate calcium transfer from network SR to the
    junctional and corbular SR compartments.

    Inputs:
        cansr : network SR calcium concentration (mM)
        cajsr : junctional SR calcium concentration (mM)
        cacsr : corbular SR calcium concentration (mM)

    Returns:
        Jtr1 : transfer flux from network SR to junctional SR
               (mM/ms)
        Jtr2 : transfer flux from network SR to corbular SR
               (mM/ms)

    Both transfer fluxes use the 120 ms time constant from
    the CellML source.
    """

    Jtr1 = (cansr - cajsr) / 120.0
    Jtr2 = (cansr - cacsr) / 120.0

    return {
        "Jtr1": Jtr1,
        "Jtr2": Jtr2
    }

'''this returns the derivatives of all components'''
def Trovato_Purkinje_Model(L=0.0164, radius=0.00175, v, time, ENa, EK, GfNa = 0.0116, GfK = 0.232, y =0.233119011214908, cass, time):
    #returns the implementation of running the Trovato Purkinje model, including ODEs and current calculations.
    #can be used in other files to produce results for the model but will only output data not graphical representations.
    #want this to return net current for the model
    #establish initial conditions and physical constants
    R, T, F, zna, zca, zk = physical_constants()
    nao, cao, ko = extracellular()
    V_cell, A_geo, A_cap, V_myo, V_nsr, V_jsr, V_csr, V_ss, V_sl = Cell_Geometry(L, radius)
    #calculate shared quantities needed for further calculations
    ENa, EK, ECa, EKs = Reverse_Potentials(R,T,F, nao, ko, cao, nasl, ksl, casl)
    CaMKb, CaMKa, dCaMKa_dt = CaMK_func(cass, time)
    vffrt = v * F**2 / (R*T)
    vfrt = v * F / (R*T)

    '''establish the current components'''

    If, IfNa, IfK, dy_dt = If_func(v, ENa, EK, y)
    INa, tm = INa_func(KmCaMK, v, CaMKa, ENa, time)
    INaL = INaL_func(tm, v, KmCaMK, CaMKa, ENa, time)
    Ito, da_dt, di1_dt, di2_dt = Ito_func(v, EK, a, i1, i2)
    Isus = Isus_func(v, EK)
    (ICaL, ICaNa, ICaK,
            dd_dt, dff_dt, dfs_dt,
            dfcaf_dt, dfcas_dt, djca_dt,
            dffp_dt, dfcafp_dt, dnca_dt)= ICaL_func(v, cass, cao, vffrt, vfrt,nao, nass, ko, kss, KmCaMK, CaMKa,d, ff, fs, fcaf, fcas, jca, ffp, fcafp, nca)
    ICaT, db_dt, dg_dt = ICaT_func(v, ECa, b, g)
    IKr, dxrf_dt, dxrs_dt = IKr_func(v, EK, ko, xrf, xrs)
    IKs, dxs1_dt, dxs2_dt = IKs_func(v, EKs, casl, xs1, xs2)
    IK1, dxk1_dt = IK1_func(v, EK, ko, xk1)
    INaCa_i, INaCa_ss = INaCa_i_func(v, R, T, F, nass, cass, nao, cao, nasl, casl, zna, zca)
    INaK =  INaK_func(v, R, T, F, nasl, ksl, nao, ko, zna, zk)
    INab = INab_func(nasl, nao, vffrt, vfrt)
    ICa_b = ICab_func(casl, cao, vffrt, vfrt)
    IpCa = IpCa_func(casl)

    '''calculate the fluxes'''
    fluxes = diffusion_fluxes(cai,casl,cass ,nai,nasl,nass, ki,ksl,kss)
    serca = serca_fluxes(KmCaMK, CaMKa,cai, casl, cansr)
    ip3 = ip3_derivatives(cass, cajsr, u)
    ryr = ryr_derivatives(Jrel1, Jrel2, ICaL, Acap, F, vss, vsl, vcsr, vnsr, vjsr, vmyo, cajsr, cacsr, Jdiff, Jgap, Jip3, Jup2, CaMKa)
    '''intracellular ion derivates'''
    ion_derivs = intracellular_ion_derivatives(X, currents, fluxes, volumes, Acap, F, cm=1.0)

    '''Membrane voltage derivateives'''
    dv_dt, vffrt, vfrt = Membrane_Potential(R,T,F,C_m, I_Na, I_NaL, I_to, I_sus, I_CaL, I_CaT, I_CaNa, I_CaK, I_Kr, I_Ks, I_f, I_K1, I_NaCa_i, I_NaCa_ss, I_NaK, I_Nab, I_pCa, I_Cab, amp, v, duration)

    return np.array([
    dv_dt,
    dCaMKt_dt,

    dm_dt,
    dhf_dt,
    dhs_dt,
    dj_dt,
    dhsp_dt,
    djp_dt,

    dmL_dt,
    dhL_dt,
    dhLp_dt,

    da_dt,
    di1_dt,
    di2_dt,

    dd_dt,
    dff_dt,
    dfs_dt,
    dfcaf_dt,
    dfcas_dt,
    djca_dt,
    dffp_dt,
    dfcafp_dt,
    dnca_dt,

    db_dt,
    dg_dt,

    dxrf_dt,
    dxrs_dt,

    dxs1_dt,
    dxs2_dt,

    dy_dt,

    dxk1_dt,

    *ion_derivs,

    dJrel1_dt,
    dJrel2_dt,

    du_dt
])
