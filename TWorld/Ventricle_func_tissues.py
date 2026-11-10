"""Generated Python translation of TWorld.cellml.
Algebraic equations and state derivatives translated from the supplied CellML.
Units follow the model: time ms, voltage mV, concentrations mM."""
import numpy as np

from scipy.integrate import solve_ivp
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

import Bidomain_functions_1D as bd

STATE_NAMES = ['v', 'camk_trap', 'camk_f_ICaL', 'camk_f_RyR', 'camk_f_PLB', 'casig_serca_trap', 'buffers_NaBj', 'buffers_NaBsl', 'buffers_TnClow', 'buffers_TnCHc', 'buffers_TnCHm', 'buffers_CaM', 'buffers_Myosin_ca', 'buffers_Myosin_mg', 'buffers_SRB', 'buffers_SLLj', 'buffers_SLLsl', 'buffers_SLHj', 'buffers_SLHsl', 'buffers_Csqn', 'casr', 'naj', 'nasl', 'nai', 'ki', 'cli', 'caj', 'casl', 'cai', 'm', 'h', 'j', 'hp', 'jp', 'm_P', 'h_P', 'j_P', 'hp_P', 'jp_P', 'mL', 'hL', 'hLp', 'd', 'ff', 'fs', 'fcaf', 'fcas', 'jca', 'ffp', 'fcafp', 'nca', 'nca_i', 'ical_pureCDI_junc', 'ical_pureCDI_sl', 'd_P', 'ff_P', 'fs_P', 'fcaf_P', 'fcas_P', 'fBPf', 'fcaBPf', 'xtos', 'ytos', 'xtof', 'ytof', 'xtos_p', 'xtof_p', 'ytos_p', 'ytof_p', 'C0', 'C1', 'C2', 'O', 'I', 'xs_junc', 'xs_sl', 'jrel_icaldep_act', 'jrel_icaldep_f1', 'jrel_icaldep_f2', 'ryr_R', 'ryr_O', 'ryr_I', 'ryr_CaRI', 'ryr_R_p', 'ryr_O_p', 'ryr_I_p', 'ryr_CaRI_p', 'contraction_XS', 'contraction_XW', 'contraction_Ca_TRPN', 'contraction_TmBlocked', 'contraction_ZETAS', 'contraction_ZETAW']
STATE_INITIAL = np.array([-87.0655895407835, 0.0414775695414329, 0.149310462149601, 0.0583425856770706, 0.0102519516890396, 0.0414775695414004, 3.5588839934477, 0.776467094025442, 0.0104491410305814, 0.124415836833344, 0.00475139089328611, 0.000326761159683871, 0.00435725142843641, 0.134587917639655, 0.00247759889601234, 0.00415486833433121, 0.00810233315864154, 0.0497605868295565, 0.0960185900338078, 1.32358809148372, 0.709622552333155, 8.88424477947199, 8.88341217642103, 8.88365804124162, 134.213097948672, 15.5745509675962, 9.20580586542242e-05, 9.62428218317655e-05, 9.66105402693763e-05, 0.00115997496446638, 0.791660456106286, 0.791587530878068, 0.61240980788246, 0.791052149435449, 0.0033402868342183, 0.647375310489836, 0.647279109941897, 0.419158183843564, 0.418669666153646, 0.000224902633151066, 0.478572465113524, 0.258793079229038, 6.00992839633468e-26, 0.999999988243592, 0.911094950010905, 0.999999988251085, 0.991477881027689, 0.999993204834445, 0.999999987868131, 0.999999987980701, 0.00528788208273223, 0.00621547539490541, 0.99834227586449, 0.99826318591148, 0.0, 0.999999987992091, 0.930411822666439, 0.999999987985882, 0.999145398038995, 0.999999987992091, 0.999999987985882, 0.000286101402392894, 0.768184915918179, 0.00028609937449475, 0.99999864687518, 0.000132590983159417, 0.000132590043205027, 0.777117031645024, 0.99999864687975, 0.997558791218494, 0.000960122333286174, 0.000804204051120587, 0.00064129227809794, 3.55901681542091e-05, 0.0311875010014239, 0.031204810574373, -9.23949112019061e-15, 0.999994847409611, 0.998725931920675, 0.821525956842405, 3.33929041065479e-08, 4.57362306340243e-09, 0.0653313479138348, 0.8216639452872, 5.00978004918837e-08, 6.8611892649539e-09, 0.0651828099424919, 0.00191507884532543, 0.00285370717430549, 0.0319778784167468, 0.992378492940848, 0.0, 0.0], dtype=float)
INITIAL_VALUES = {'celltype': 0.0, 'sex': 0.0, 'INa_multiplier': 1.0, 'INaL_multiplier': 1.0, 'INab_multiplier': 1.0, 'INaK_multiplier': 1.0, 'IClCa_multiplier': 1.0, 'IClb_multiplier': 1.0, 'Itos_multiplier': 1.0, 'Itof_multiplier': 1.0, 'IKb_multiplier': 1.0, 'ICaLPCa_multiplier': 1.0, 'ICab_multiplier': 1.0, 'Jrel_multiplier': 1.0, 'Jup_multiplier': 1.0, 'INaCa_multiplier': 1.0, 'IKr_multiplier': 1.0, 'IKs_multiplier': 1.0, 'IK1_multiplier': 1.0, 'IpCa_multiplier': 1.0, 'cellLength': 100.0, 'cellRadius': 10.25, 'nao': 140.0, 'cao': 1.8, 'ko': 5.0, 'clo': 150.0, 'R': 8314.0, 'T': 310.0, 'F': 96485.0, 'zna': 1.0, 'zca': 2.0, 'zk': 1.0, 'zcl': -1.0, 'PNaK': 0.01833, 'Cmem': 1.381e-10, 'v': -87.0655895407835, 'Fjunc': 0.11, 'i_Stim_Start': 0.0, 'i_Stim_Amplitude': -53.0, 'i_Stim_Period': 1000.0, 'i_Stim_PulseDuration': 1.0, 'Whole_cell_PP1': 0.13698, 'camk_trap': 0.0414775695414329, 'CaMK0': 0.1, 'Km_CaMK_Ca': 0.0075, 'CaMK__alpha': 0.05, 'CaMK__beta': 0.00068, 'tau_plb': 100000.0, 'tau_ryr': 10000.0, 'K_Phos_CaMK': 0.35, 'camk_f_ICaL': 0.149310462149601, 'camk_f_RyR': 0.0583425856770706, 'camk_f_PLB': 0.0102519516890396, 'alpha_serca': 0.05, 'casig_serca_trap': 0.0414775695414004, 'kon_na': 0.0001, 'koff_na': 0.001, 'kon_tnchca': 2.37, 'koff_tnchca': 3.2e-05, 'kon_tnchmg': 0.003, 'koff_tnchmg': 0.00333, 'kon_cam': 34.0, 'koff_cam': 0.238, 'kon_myoca': 13.8, 'koff_myoca': 0.00046, 'kon_myomg': 0.0157, 'koff_myomg': 5.7e-05, 'kon_sr': 100.0, 'koff_sr': 0.06, 'kon_sll': 100.0, 'koff_sll': 1.3, 'kon_slh': 100.0, 'koff_slh': 0.03, 'kon_csqn': 100.0, 'koff_csqn': 65.0, 'Bmax_Naj': 7.561, 'Bmax_Nasl': 1.65, 'Bmax_TnClow': 0.07, 'Bmax_TnChigh': 0.14, 'Bmax_myosin': 0.14, 'Bmax_SR': 0.01785854, 'buffers_NaBj': 3.5588839934477, 'buffers_NaBsl': 0.776467094025442, 'buffers_TnClow': 0.0104491410305814, 'buffers_TnCHc': 0.124415836833344, 'buffers_TnCHm': 0.00475139089328611, 'buffers_CaM': 0.000326761159683871, 'buffers_Myosin_ca': 0.00435725142843641, 'buffers_Myosin_mg': 0.134587917639655, 'buffers_SRB': 0.00247759889601234, 'buffers_SLLj': 0.00415486833433121, 'buffers_SLLsl': 0.00810233315864154, 'buffers_SLHj': 0.0497605868295565, 'buffers_SLHsl': 0.0960185900338078, 'buffers_Csqn': 1.32358809148372, 'J_na_juncsl': 1.83127823220608e-14, 'J_na_slmyo': 1.638627922219794e-12, 'J_ca_juncsl': 3.260674633581688e-13, 'J_ca_slmyo': 1.341273673480337e-11, 'naj': 8.88424477947199, 'nasl': 8.88341217642103, 'nai': 8.88365804124162, 'mgi': 0.5, 'ki': 134.213097948672, 'cli': 15.5745509675962, 'casr': 0.709622552333155, 'caj': 9.20580586542242e-05, 'casl': 9.62428218317655e-05, 'cai': 9.66105402693763e-05, 'GNa_b': 22.08788, 'm': 0.00115997496446638, 'h': 0.791660456106286, 'j': 0.791587530878068, 'hp': 0.61240980788246, 'jp': 0.791052149435449, 'm_P': 0.0033402868342183, 'h_P': 0.647375310489836, 'j_P': 0.647279109941897, 'hp_P': 0.419158183843564, 'jp_P': 0.418669666153646, 'GNaL_b': 0.04229, 'mL': 0.000224902633151066, 'thL': 145.0, 'hL': 0.478572465113524, 'hLp': 0.258793079229038, 'd': 6.00992839633468e-26, 'Aff': 0.52477, 'ff': 0.999999988243592, 'fs': 0.911094950010905, 'fcaf': 0.999999988251085, 'fcas': 0.991477881027689, 'tjca': 66.0, 'jca': 0.999993204834445, 'ffp': 0.999999987868131, 'fcafp': 0.999999987980701, 'Kmn': 0.00222, 'k2n': 957.85903, 'nca': 0.00528788208273223, 'nca_i': 0.00621547539490541, 'dielConstant': 74.0, 'PCa_b': 0.00015768, 'd_P': 0.0, 'ff_P': 0.999999987992091, 'fs_P': 0.930411822666439, 'fcaf_P': 0.999999987985882, 'fcas_P': 0.999145398038995, 'fBPf': 0.999999987992091, 'fcaBPf': 0.999999987985882, 'ICaL_fractionSS': 0.8, 'rateRecovery': 0.02313, 'ical_pureCDI_junc': 0.99834227586449, 'ical_pureCDI_sl': 0.99826318591148, 'xtos': 0.000286101402392894, 'ytos': 0.768184915918179, 'xtof': 0.00028609937449475, 'ytof': 0.99999864687518, 'xtos_p': 0.000132590983159417, 'xtof_p': 0.000132590043205027, 'ytos_p': 0.777117031645024, 'ytof_p': 0.99999864687975, 'GKr_b': 0.043, 'C0': 0.997558791218494, 'C1': 0.000960122333286174, 'C2': 0.000804204051120587, 'I': 3.55901681542091e-05, 'O': 0.00064129227809794, 'alpha_1': 0.154375, 'IKr__beta_1': 0.1911, 'GKs_b': 2.97002, 'gKs_factor': 0.01, 'xs_junc': 0.0311875010014239, 'xs_sl': 0.031204810574373, 'GK1_b': 0.6992, 'INaCa_fractionSS': 0.31, 'kna1': 11.9712, 'kna2': 2.76, 'kna3': 88.767, 'kasymm': 19.4258, 'wna': 32978.0, 'wca': 51756.0, 'wnaca': 2776.3, 'kcaon': 3416400.0, 'kcaoff': 3853.2, 'qna': 0.6718, 'qca': 0.0955, 'KmCaAct': 0.00015, 'Gncx_b': 0.00179, 'IbarNaK_b': 2.10774, 'KmNaip': 11.0, 'KmNaip_PKA': 8.4615, 'KmKo': 1.5, 'GNab_b': 0.000594, 'GCab_b': 0.000515575, 'GKb_b': 0.010879, 'GClCa_b': 0.01615, 'GClb_b': 0.00241, 'KdClCa': 0.1, 'IbarSLCaP_b': 0.02064, 'KmPCa': 0.0005, 'Q10SLCaP': 2.35, 'koCa': 23.87221, 'kom': 0.16219, 'kiCa': 0.39871, 'kim': 0.04311, 'ec50SR': 0.75385, 'steepnessCaSR': 5.09473, 'caExpFactor': 0.68655, 'caTransFactor': 0.94428, 'caExpFactor2': 2.06273, 'caTransFactor2': 0.52967, 'directRelMidpoint': 0.95271, 'bt': 12.4767, 'jrel_icaldep_act': -9.23949112019061e-15, 'tauInact': 64.11202, 'jrel_icaldep_f1': 0.999994847409611, 'tauInact2': 119.48978, 'jrel_icaldep_f2': 0.998725931920675, 'MaxSR': 15.0, 'MinSR': 1.0, 'ecCaI': 0.001, 'steepnessCaI': 5.93447, 'minCaI': 0.93249, 'maxCaI': 30.13294, 'baseRateCaI': 0.00030232, 'CI_to_RI': 0.00248, 'ryr_R': 0.821525956842405, 'ryr_O': 3.33929041065479e-08, 'ryr_I': 4.57362306340243e-09, 'ryr_CaRI': 0.0653313479138348, 'ryr_R_p': 0.8216639452872, 'ryr_O_p': 5.00978004918837e-08, 'ryr_I_p': 6.8611892649539e-09, 'ryr_CaRI_p': 0.0651828099424919, 'Q10SRCaP': 2.6, 'Vmax_SRCaP_b': 0.00543, 'Kmr': 2.31442, 'hillSRCaP': 1.02809, 'Kmf': 0.00030672, 'Km_SERCA_Ca': 0.4, 'Max_Vmax_SERCA_Ca': 1.11142, 'fracTnIpo': 0.0031, 'perm50': 0.35, 'TRPN_n': 1.65, 'koff': 0.07854, 'dr': 0.25, 'wfrac': 0.5, 'TOT_A': 25.0, 'ktm_unblock': 0.02626, 'beta_0': 2.3, 'contraction__beta_1': -2.4, 'gamma': 0.0085, 'gamma_wu': 0.615, 'phi': 2.23, 'nperm': 2.036, 'Tref': 80.0, 'nu': 10.15996, 'mu': 3.94046, 'lambda_': 1.0, 'lambda_min': 0.87, 'lambda_max': 1.2, 'lambda_rate': 0.0, 'contraction_TmBlocked': 0.992378492940848, 'contraction_XW': 0.00285370717430549, 'contraction_XS': 0.00191507884532543, 'contraction_ZETAS': 0.0, 'contraction_ZETAW': 0.0, 'contraction_Ca_TRPN': 0.0319778784167468, 'fINa_PKA': 0.0, 'fICaL_PKA': 0.0, 'fINaK_PKA': 0.0, 'fIKs_PKA': 0.0, 'fPLB_PKA': 0.0, 'fTnI_PKA': 0.0, 'fMyBPC_PKA': 0.0}

def initial_state():
    """Return a fresh copy of the 93 CellML state initial values."""
    return STATE_INITIAL.copy()

def evaluate(time, state, parameters=None):
    """Return (state derivatives, all evaluated variables) at one instant."""
    parameters = {} if parameters is None else dict(parameters)
    values = dict(INITIAL_VALUES)
    values.update(parameters)
    values["time"] = float(time)
    celltype = values['celltype']
    sex = values['sex']
    INa_multiplier = values['INa_multiplier']
    INaL_multiplier = values['INaL_multiplier']
    INab_multiplier = values['INab_multiplier']
    INaK_multiplier = values['INaK_multiplier']
    IClCa_multiplier = values['IClCa_multiplier']
    IClb_multiplier = values['IClb_multiplier']
    Itos_multiplier = values['Itos_multiplier']
    Itof_multiplier = values['Itof_multiplier']
    IKb_multiplier = values['IKb_multiplier']
    ICaLPCa_multiplier = values['ICaLPCa_multiplier']
    ICab_multiplier = values['ICab_multiplier']
    Jrel_multiplier = values['Jrel_multiplier']
    Jup_multiplier = values['Jup_multiplier']
    INaCa_multiplier = values['INaCa_multiplier']
    IKr_multiplier = values['IKr_multiplier']
    IKs_multiplier = values['IKs_multiplier']
    IK1_multiplier = values['IK1_multiplier']
    IpCa_multiplier = values['IpCa_multiplier']
    cellLength = values['cellLength']
    cellRadius = values['cellRadius']
    nao = values['nao']
    cao = values['cao']
    ko = values['ko']
    clo = values['clo']
    R = values['R']
    T = values['T']
    F = values['F']
    zna = values['zna']
    zca = values['zca']
    zk = values['zk']
    zcl = values['zcl']
    PNaK = values['PNaK']
    Cmem = values['Cmem']
    Fjunc = values['Fjunc']
    i_Stim_Start = values['i_Stim_Start']
    i_Stim_Amplitude = values['i_Stim_Amplitude']
    i_Stim_Period = values['i_Stim_Period']
    i_Stim_PulseDuration = values['i_Stim_PulseDuration']
    Whole_cell_PP1 = values['Whole_cell_PP1']
    CaMK0 = values['CaMK0']
    Km_CaMK_Ca = values['Km_CaMK_Ca']
    CaMK__alpha = values['CaMK__alpha']
    CaMK__beta = values['CaMK__beta']
    tau_plb = values['tau_plb']
    tau_ryr = values['tau_ryr']
    K_Phos_CaMK = values['K_Phos_CaMK']
    alpha_serca = values['alpha_serca']
    kon_na = values['kon_na']
    koff_na = values['koff_na']
    kon_tnchca = values['kon_tnchca']
    koff_tnchca = values['koff_tnchca']
    kon_tnchmg = values['kon_tnchmg']
    koff_tnchmg = values['koff_tnchmg']
    kon_cam = values['kon_cam']
    koff_cam = values['koff_cam']
    kon_myoca = values['kon_myoca']
    koff_myoca = values['koff_myoca']
    kon_myomg = values['kon_myomg']
    koff_myomg = values['koff_myomg']
    kon_sr = values['kon_sr']
    koff_sr = values['koff_sr']
    kon_sll = values['kon_sll']
    koff_sll = values['koff_sll']
    kon_slh = values['kon_slh']
    koff_slh = values['koff_slh']
    kon_csqn = values['kon_csqn']
    koff_csqn = values['koff_csqn']
    Bmax_Naj = values['Bmax_Naj']
    Bmax_Nasl = values['Bmax_Nasl']
    Bmax_TnClow = values['Bmax_TnClow']
    Bmax_TnChigh = values['Bmax_TnChigh']
    Bmax_myosin = values['Bmax_myosin']
    Bmax_SR = values['Bmax_SR']
    J_na_juncsl = values['J_na_juncsl']
    J_na_slmyo = values['J_na_slmyo']
    J_ca_juncsl = values['J_ca_juncsl']
    J_ca_slmyo = values['J_ca_slmyo']
    mgi = values['mgi']
    GNa_b = values['GNa_b']
    GNaL_b = values['GNaL_b']
    thL = values['thL']
    Aff = values['Aff']
    tjca = values['tjca']
    Kmn = values['Kmn']
    k2n = values['k2n']
    dielConstant = values['dielConstant']
    PCa_b = values['PCa_b']
    ICaL_fractionSS = values['ICaL_fractionSS']
    rateRecovery = values['rateRecovery']
    GKr_b = values['GKr_b']
    alpha_1 = values['alpha_1']
    IKr__beta_1 = values['IKr__beta_1']
    GKs_b = values['GKs_b']
    gKs_factor = values['gKs_factor']
    GK1_b = values['GK1_b']
    INaCa_fractionSS = values['INaCa_fractionSS']
    kna1 = values['kna1']
    kna2 = values['kna2']
    kna3 = values['kna3']
    kasymm = values['kasymm']
    wna = values['wna']
    wca = values['wca']
    wnaca = values['wnaca']
    kcaon = values['kcaon']
    kcaoff = values['kcaoff']
    qna = values['qna']
    qca = values['qca']
    KmCaAct = values['KmCaAct']
    Gncx_b = values['Gncx_b']
    IbarNaK_b = values['IbarNaK_b']
    KmNaip = values['KmNaip']
    KmNaip_PKA = values['KmNaip_PKA']
    KmKo = values['KmKo']
    GNab_b = values['GNab_b']
    GCab_b = values['GCab_b']
    GKb_b = values['GKb_b']
    GClCa_b = values['GClCa_b']
    GClb_b = values['GClb_b']
    KdClCa = values['KdClCa']
    IbarSLCaP_b = values['IbarSLCaP_b']
    KmPCa = values['KmPCa']
    Q10SLCaP = values['Q10SLCaP']
    koCa = values['koCa']
    kom = values['kom']
    kiCa = values['kiCa']
    kim = values['kim']
    ec50SR = values['ec50SR']
    steepnessCaSR = values['steepnessCaSR']
    caExpFactor = values['caExpFactor']
    caTransFactor = values['caTransFactor']
    caExpFactor2 = values['caExpFactor2']
    caTransFactor2 = values['caTransFactor2']
    directRelMidpoint = values['directRelMidpoint']
    bt = values['bt']
    tauInact = values['tauInact']
    tauInact2 = values['tauInact2']
    MaxSR = values['MaxSR']
    MinSR = values['MinSR']
    ecCaI = values['ecCaI']
    steepnessCaI = values['steepnessCaI']
    minCaI = values['minCaI']
    maxCaI = values['maxCaI']
    baseRateCaI = values['baseRateCaI']
    CI_to_RI = values['CI_to_RI']
    Q10SRCaP = values['Q10SRCaP']
    Vmax_SRCaP_b = values['Vmax_SRCaP_b']
    Kmr = values['Kmr']
    hillSRCaP = values['hillSRCaP']
    Kmf = values['Kmf']
    Km_SERCA_Ca = values['Km_SERCA_Ca']
    Max_Vmax_SERCA_Ca = values['Max_Vmax_SERCA_Ca']
    fracTnIpo = values['fracTnIpo']
    perm50 = values['perm50']
    TRPN_n = values['TRPN_n']
    koff = values['koff']
    dr = values['dr']
    wfrac = values['wfrac']
    TOT_A = values['TOT_A']
    ktm_unblock = values['ktm_unblock']
    beta_0 = values['beta_0']
    contraction__beta_1 = values['contraction__beta_1']
    gamma = values['gamma']
    gamma_wu = values['gamma_wu']
    phi = values['phi']
    nperm = values['nperm']
    Tref = values['Tref']
    nu = values['nu']
    mu = values['mu']
    lambda_ = values['lambda_']
    lambda_min = values['lambda_min']
    lambda_max = values['lambda_max']
    lambda_rate = values['lambda_rate']
    fINa_PKA = values['fINa_PKA']
    fICaL_PKA = values['fICaL_PKA']
    fINaK_PKA = values['fINaK_PKA']
    fIKs_PKA = values['fIKs_PKA']
    fPLB_PKA = values['fPLB_PKA']
    fTnI_PKA = values['fTnI_PKA']
    fMyBPC_PKA = values['fMyBPC_PKA']
    state = np.asarray(state, dtype=float)
    if state.shape != STATE_INITIAL.shape:
        raise ValueError(f"state must have shape {STATE_INITIAL.shape}, got {state.shape}")
    v = state[0]
    values['v'] = v
    camk_trap = state[1]
    values['camk_trap'] = camk_trap
    camk_f_ICaL = state[2]
    values['camk_f_ICaL'] = camk_f_ICaL
    camk_f_RyR = state[3]
    values['camk_f_RyR'] = camk_f_RyR
    camk_f_PLB = state[4]
    values['camk_f_PLB'] = camk_f_PLB
    casig_serca_trap = state[5]
    values['casig_serca_trap'] = casig_serca_trap
    buffers_NaBj = state[6]
    values['buffers_NaBj'] = buffers_NaBj
    buffers_NaBsl = state[7]
    values['buffers_NaBsl'] = buffers_NaBsl
    buffers_TnClow = state[8]
    values['buffers_TnClow'] = buffers_TnClow
    buffers_TnCHc = state[9]
    values['buffers_TnCHc'] = buffers_TnCHc
    buffers_TnCHm = state[10]
    values['buffers_TnCHm'] = buffers_TnCHm
    buffers_CaM = state[11]
    values['buffers_CaM'] = buffers_CaM
    buffers_Myosin_ca = state[12]
    values['buffers_Myosin_ca'] = buffers_Myosin_ca
    buffers_Myosin_mg = state[13]
    values['buffers_Myosin_mg'] = buffers_Myosin_mg
    buffers_SRB = state[14]
    values['buffers_SRB'] = buffers_SRB
    buffers_SLLj = state[15]
    values['buffers_SLLj'] = buffers_SLLj
    buffers_SLLsl = state[16]
    values['buffers_SLLsl'] = buffers_SLLsl
    buffers_SLHj = state[17]
    values['buffers_SLHj'] = buffers_SLHj
    buffers_SLHsl = state[18]
    values['buffers_SLHsl'] = buffers_SLHsl
    buffers_Csqn = state[19]
    values['buffers_Csqn'] = buffers_Csqn
    casr = state[20]
    values['casr'] = casr
    naj = state[21]
    values['naj'] = naj
    nasl = state[22]
    values['nasl'] = nasl
    nai = state[23]
    values['nai'] = nai
    ki = state[24]
    values['ki'] = ki
    cli = state[25]
    values['cli'] = cli
    caj = state[26]
    values['caj'] = caj
    casl = state[27]
    values['casl'] = casl
    cai = state[28]
    values['cai'] = cai
    m = state[29]
    values['m'] = m
    h = state[30]
    values['h'] = h
    j = state[31]
    values['j'] = j
    hp = state[32]
    values['hp'] = hp
    jp = state[33]
    values['jp'] = jp
    m_P = state[34]
    values['m_P'] = m_P
    h_P = state[35]
    values['h_P'] = h_P
    j_P = state[36]
    values['j_P'] = j_P
    hp_P = state[37]
    values['hp_P'] = hp_P
    jp_P = state[38]
    values['jp_P'] = jp_P
    mL = state[39]
    values['mL'] = mL
    hL = state[40]
    values['hL'] = hL
    hLp = state[41]
    values['hLp'] = hLp
    d = state[42]
    values['d'] = d
    ff = state[43]
    values['ff'] = ff
    fs = state[44]
    values['fs'] = fs
    fcaf = state[45]
    values['fcaf'] = fcaf
    fcas = state[46]
    values['fcas'] = fcas
    jca = state[47]
    values['jca'] = jca
    ffp = state[48]
    values['ffp'] = ffp
    fcafp = state[49]
    values['fcafp'] = fcafp
    nca = state[50]
    values['nca'] = nca
    nca_i = state[51]
    values['nca_i'] = nca_i
    ical_pureCDI_junc = state[52]
    values['ical_pureCDI_junc'] = ical_pureCDI_junc
    ical_pureCDI_sl = state[53]
    values['ical_pureCDI_sl'] = ical_pureCDI_sl
    d_P = state[54]
    values['d_P'] = d_P
    ff_P = state[55]
    values['ff_P'] = ff_P
    fs_P = state[56]
    values['fs_P'] = fs_P
    fcaf_P = state[57]
    values['fcaf_P'] = fcaf_P
    fcas_P = state[58]
    values['fcas_P'] = fcas_P
    fBPf = state[59]
    values['fBPf'] = fBPf
    fcaBPf = state[60]
    values['fcaBPf'] = fcaBPf
    xtos = state[61]
    values['xtos'] = xtos
    ytos = state[62]
    values['ytos'] = ytos
    xtof = state[63]
    values['xtof'] = xtof
    ytof = state[64]
    values['ytof'] = ytof
    xtos_p = state[65]
    values['xtos_p'] = xtos_p
    xtof_p = state[66]
    values['xtof_p'] = xtof_p
    ytos_p = state[67]
    values['ytos_p'] = ytos_p
    ytof_p = state[68]
    values['ytof_p'] = ytof_p
    C0 = state[69]
    values['C0'] = C0
    C1 = state[70]
    values['C1'] = C1
    C2 = state[71]
    values['C2'] = C2
    O = state[72]
    values['O'] = O
    I = state[73]
    values['I'] = I
    xs_junc = state[74]
    values['xs_junc'] = xs_junc
    xs_sl = state[75]
    values['xs_sl'] = xs_sl
    jrel_icaldep_act = state[76]
    values['jrel_icaldep_act'] = jrel_icaldep_act
    jrel_icaldep_f1 = state[77]
    values['jrel_icaldep_f1'] = jrel_icaldep_f1
    jrel_icaldep_f2 = state[78]
    values['jrel_icaldep_f2'] = jrel_icaldep_f2
    ryr_R = state[79]
    values['ryr_R'] = ryr_R
    ryr_O = state[80]
    values['ryr_O'] = ryr_O
    ryr_I = state[81]
    values['ryr_I'] = ryr_I
    ryr_CaRI = state[82]
    values['ryr_CaRI'] = ryr_CaRI
    ryr_R_p = state[83]
    values['ryr_R_p'] = ryr_R_p
    ryr_O_p = state[84]
    values['ryr_O_p'] = ryr_O_p
    ryr_I_p = state[85]
    values['ryr_I_p'] = ryr_I_p
    ryr_CaRI_p = state[86]
    values['ryr_CaRI_p'] = ryr_CaRI_p
    contraction_XS = state[87]
    values['contraction_XS'] = contraction_XS
    contraction_XW = state[88]
    values['contraction_XW'] = contraction_XW
    contraction_Ca_TRPN = state[89]
    values['contraction_Ca_TRPN'] = contraction_Ca_TRPN
    contraction_TmBlocked = state[90]
    values['contraction_TmBlocked'] = contraction_TmBlocked
    contraction_ZETAS = state[91]
    values['contraction_ZETAS'] = contraction_ZETAS
    contraction_ZETAW = state[92]
    values['contraction_ZETAW'] = contraction_ZETAW
    A = ((TOT_A * dr) / (((1 - dr) * wfrac) + dr))
    values['A'] = A
    Afcaf = (0.3 + (0.6 / (1 + np.exp(((v - 9.24247) / 27.96201)))))
    values['Afcaf'] = Afcaf
    Afs = (1 - Aff)
    values['Afs'] = Afs
    Bmax_CaM = ((0.024 / 1.1) if (sex == 1) else ((0.024 * 1.1) if (sex == 2) else 0.024))
    values['Bmax_CaM'] = Bmax_CaM
    CaMKIILeakMultiplier = (1 + (2 * camk_f_RyR))
    values['CaMKIILeakMultiplier'] = CaMKIILeakMultiplier
    ECaj = (((R * T) / (zca * F)) * np.log((cao / caj)))
    values['ECaj'] = ECaj
    ECasl = (((R * T) / (zca * F)) * np.log((cao / casl)))
    values['ECasl'] = ECasl
    ECl = (((R * T) / (zcl * F)) * np.log((clo / cli)))
    values['ECl'] = ECl
    EK = (((R * T) / (zk * F)) * np.log((ko / ki)))
    values['EK'] = EK
    EKs = (((R * T) / (zk * F)) * np.log(((ko + (PNaK * nao)) / (ki + (PNaK * nai)))))
    values['EKs'] = EKs
    ENaj = (((R * T) / (zna * F)) * np.log((nao / naj)))
    values['ENaj'] = ENaj
    ENasl = (((R * T) / (zna * F)) * np.log((nao / nasl)))
    values['ENasl'] = ENasl
    Fsl = (1.0 - Fjunc)
    values['Fsl'] = Fsl
    GCab = (ICab_multiplier * GCab_b)
    values['GCab'] = GCab
    GClCa = (IClCa_multiplier * GClCa_b)
    values['GClCa'] = GClCa
    GClb = (IClb_multiplier * GClb_b)
    values['GClb'] = GClb
    GKb = (GKb_b * IKb_multiplier)
    values['GKb'] = GKb
    GNa = (GNa_b * INa_multiplier)
    values['GNa'] = GNa
    GNab = (INab_multiplier * GNab_b)
    values['GNab'] = GNab
    Gto_fast = ((0.29856 * Itof_multiplier) if (celltype == 1) else ((0.14928 * Itof_multiplier) if (celltype == 2) else (0.01276 * Itof_multiplier)))
    values['Gto_fast'] = Gto_fast
    Gto_slow = ((0.02036 * Itos_multiplier) if (celltype == 1) else ((0.04632 * Itos_multiplier) if (celltype == 2) else (0.07210 * Itos_multiplier)))
    values['Gto_slow'] = Gto_slow
    IK1__celltype_factor = (1.1 if (celltype == 1) else (1.3 if (celltype == 2) else 1.0))
    values['IK1__celltype_factor'] = IK1__celltype_factor
    IK1__sex_factor = (1.07 if (sex == 1) else (0.92 if (sex == 2) else 1.0))
    values['IK1__sex_factor'] = IK1__sex_factor
    IKr__celltype_factor = (1.25 if (celltype == 1) else (0.7 if (celltype == 2) else 1.0))
    values['IKr__celltype_factor'] = IKr__celltype_factor
    IKr__sex_factor = (1.1 if (sex == 1) else (0.87 if (sex == 2) else 1.0))
    values['IKr__sex_factor'] = IKr__sex_factor
    IKs__celltype_factor = (1.4 if (celltype == 1) else (0.5 if (celltype == 2) else 1.0))
    values['IKs__celltype_factor'] = IKs__celltype_factor
    IKs__sex_factor = (1.09 if (sex == 1) else (0.905 if (sex == 2) else 1.0))
    values['IKs__sex_factor'] = IKs__sex_factor
    INaCa__celltype_factor = (1.1 if (celltype == 1) else (1.4 if (celltype == 2) else 1.0))
    values['INaCa__celltype_factor'] = INaCa__celltype_factor
    INaCa__sex_factor = ((1.0 / 1.0724) if (sex == 1) else (1.0724 if (sex == 2) else 1.0))
    values['INaCa__sex_factor'] = INaCa__sex_factor
    IbarNaK = (INaK_multiplier * IbarNaK_b)
    values['IbarNaK'] = IbarNaK
    IbarSLCaP = ((IbarSLCaP_b * 0.8 * IpCa_multiplier) if (sex == 1) else ((IbarSLCaP_b * 1.28 * IpCa_multiplier) if (sex == 2) else (IbarSLCaP_b * IpCa_multiplier)))
    values['IbarSLCaP'] = IbarSLCaP
    Ii = ((0.5 * (nai + ki + cli + (4 * cai))) / 1000)
    values['Ii'] = Ii
    Io = ((0.5 * (nao + ko + clo + (4 * cao))) / 1000)
    values['Io'] = Io
    Iss = ((0.5 * (nasl + ki + cli + (4 * casl))) / 1000)
    values['Iss'] = Iss
    Istim = 0
    values['Istim'] = Istim
    Jrel_ICaLdep = (0.00174 * jrel_icaldep_act * jrel_icaldep_f1 * jrel_icaldep_f2)
    values['Jrel_ICaLdep'] = Jrel_ICaLdep
    Kmf_p = (Kmf * 0.5)
    values['Kmf_p'] = Kmf_p
    PCa = ((PCa_b * ICaLPCa_multiplier * 1.025) if (celltype == 1) else ((PCa_b * ICaLPCa_multiplier * 1.1) if (celltype == 2) else (PCa_b * ICaLPCa_multiplier)))
    values['PCa'] = PCa
    PKAForceMultiplier = (1.0 + (0.26 * fMyBPC_PKA))
    values['PKAForceMultiplier'] = PKAForceMultiplier
    PP1_tot = Whole_cell_PP1
    values['PP1_tot'] = PP1_tot
    Qpow = ((T - 310) / 10)
    values['Qpow'] = Qpow
    RIcleft = ((((1 - ryr_R) - ryr_O) - ryr_I) - ryr_CaRI)
    values['RIcleft'] = RIcleft
    RIcleft_p = ((((1 - ryr_R_p) - ryr_O_p) - ryr_I_p) - ryr_CaRI_p)
    values['RIcleft_p'] = RIcleft_p
    Vmax_SRCaP = ((1.2 * Vmax_SRCaP_b * Jup_multiplier) if (celltype == 1) else (Vmax_SRCaP_b * Jup_multiplier))
    values['Vmax_SRCaP'] = Vmax_SRCaP
    XSSS = (dr * 0.5)
    values['XSSS'] = XSSS
    XU = (((1 - contraction_TmBlocked) - contraction_XW) - contraction_XS)
    values['XU'] = XU
    XWSS = ((1 - dr) * wfrac * 0.5)
    values['XWSS'] = XWSS
    a_rel = ((1.25 * bt) / 1)
    values['a_rel'] = a_rel
    ah = (0 if (v >= (-40)) else (0.057 * np.exp(((-(v + 80)) / 6.8))))
    values['ah'] = ah
    aj = (0 if (v >= (-40)) else (((((-(2.5428e4)) * np.exp((0.2444 * v))) - ((6.948e-6) * np.exp(((-0.04391) * v)))) * (v + 37.78)) / (1 + np.exp((0.311 * (v + 79.23))))))
    values['aj'] = aj
    allo_i = (1 / (1 + ((KmCaAct / casl) ** 2)))
    values['allo_i'] = allo_i
    allo_ss = (1 / (1 + ((KmCaAct / caj) ** 2)))
    values['allo_ss'] = allo_ss
    bh = ((0.77 / (0.13 * (1 + np.exp(((-(v + 10.66)) / 11.1))))) if (v >= (-40)) else ((2.7 * np.exp((0.079 * v))) + ((3.1e5) * np.exp((0.3485 * v)))))
    values['bh'] = bh
    bj = (((0.6 * np.exp((0.057 * v))) / (1 + np.exp(((-0.1) * (v + 32))))) if (v >= (-40)) else ((0.02424 * np.exp(((-0.01052) * v))) / (1 + np.exp(((-0.1378) * (v + 40.14))))))
    values['bj'] = bj
    bound = ((CaMK0 * (1 - camk_trap)) / (1 + (Km_CaMK_Ca / caj)))
    values['bound'] = bound
    bound_serca = ((CaMK0 * (1 - casig_serca_trap)) / (1 + (Km_CaMK_Ca / caj)))
    values['bound_serca'] = bound_serca
    caTransFactor2p = (caTransFactor2 * 1.5)
    values['caTransFactor2p'] = caTransFactor2p
    constA = ((1.82e6) * ((dielConstant * T) ** (-1.5)))
    values['constA'] = constA
    dPss = min((1.0323 * np.exp(((-1.0553) * np.exp(((-0.0810) * (v + 12.62483)))))), 1)
    values['dPss'] = dPss
    d_buffers_Myosin_ca = ((kon_myoca * cai * ((Bmax_myosin - buffers_Myosin_ca) - buffers_Myosin_mg)) - (koff_myoca * buffers_Myosin_ca))
    values['d_buffers_Myosin_ca'] = d_buffers_Myosin_ca
    d_buffers_NaBj = ((kon_na * naj * (Bmax_Naj - buffers_NaBj)) - (koff_na * buffers_NaBj))
    values['d_buffers_NaBj'] = d_buffers_NaBj
    d_buffers_NaBsl = ((kon_na * nasl * (Bmax_Nasl - buffers_NaBsl)) - (koff_na * buffers_NaBsl))
    values['d_buffers_NaBsl'] = d_buffers_NaBsl
    d_buffers_SRB = ((kon_sr * cai * (Bmax_SR - buffers_SRB)) - (koff_sr * buffers_SRB))
    values['d_buffers_SRB'] = d_buffers_SRB
    d_buffers_TnCHc = ((kon_tnchca * cai * ((Bmax_TnChigh - buffers_TnCHc) - buffers_TnCHm)) - (koff_tnchca * buffers_TnCHc))
    values['d_buffers_TnCHc'] = d_buffers_TnCHc
    dss = min((1.0763 * np.exp(((-1.0070) * np.exp(((-0.0829) * (v + 3.62483)))))), 1)
    values['dss'] = dss
    dti_develop = (1.354 + ((1.0e-4) / (np.exp(((v - 167.4) / 15.89)) + np.exp(((-(v - 12.23)) / 0.2154)))))
    values['dti_develop'] = dti_develop
    dti_recover = (1.0 - (0.5 / (1.0 + np.exp(((v + 70.0) / 20.0)))))
    values['dti_recover'] = dti_recover
    fICaLP = fICaL_PKA
    values['fICaLP'] = fICaLP
    fICaLp = camk_f_ICaL
    values['fICaLp'] = fICaLp
    fINaLp = camk_f_RyR
    values['fINaLp'] = fINaLp
    fINa_P = fINa_PKA
    values['fINa_P'] = fINa_P
    fINap = camk_f_RyR
    values['fINap'] = fINap
    fItop = camk_f_RyR
    values['fItop'] = fItop
    fPKA_TnI = (1.45 - ((0.45 * (1 - fTnI_PKA)) / (1 - fracTnIpo)))
    values['fPKA_TnI'] = fPKA_TnI
    fnak = (0.75 + ((0.00375 - (((140 - nao) / 50) * 0.001)) * v))
    values['fnak'] = fnak
    fss = (1 / (1 + np.exp(((v + 19.58) / 3.696))))
    values['fss'] = fss
    fss_P = (1.0 / (1.0 + np.exp(((v + 25.58) / 3.696))))
    values['fss_P'] = fss_P
    gamma_rate = ((gamma * contraction_ZETAS) if (contraction_ZETAS > 0) else ((gamma * ((-contraction_ZETAS) - 1)) if (contraction_ZETAS < (-1)) else 0.0))
    values['gamma_rate'] = gamma_rate
    gamma_rate_w = (gamma_wu * np.abs(contraction_ZETAW))
    values['gamma_rate_w'] = gamma_rate_w
    h10_i = (kasymm + 1 + ((nao / kna1) * (1 + (nao / kna2))))
    values['h10_i'] = h10_i
    h10_ss = (kasymm + 1 + ((nao / kna1) * (1 + (nao / kna2))))
    values['h10_ss'] = h10_ss
    h4_i = (1 + ((nai / kna1) * (1 + (nai / kna2))))
    values['h4_i'] = h4_i
    h4_ss = (1 + ((nasl / kna1) * (1 + (nasl / kna2))))
    values['h4_ss'] = h4_ss
    hLss = (1 / (1 + np.exp(((v + 87.61) / 7.488))))
    values['hLss'] = hLss
    hLssp = (1 / (1 + np.exp(((v + 93.81) / 7.488))))
    values['hLssp'] = hLssp
    hca = np.exp(((qca * (v - 8.3117) * F) / (R * T)))
    values['hca'] = hca
    hna = np.exp(((qna * (v - 8.3117) * F) / (R * T)))
    values['hna'] = hna
    hss = (1 / ((1 + np.exp(((v + 71.55) / 7.43))) ** 2))
    values['hss'] = hss
    hss_P = (1 / ((1 + np.exp(((v + 76.55) / 7.43))) ** 2))
    values['hss_P'] = hss_P
    hssp = (1 / ((1 + np.exp(((v + 77.55) / 7.43))) ** 2))
    values['hssp'] = hssp
    hssp_P = (1 / ((1 + np.exp(((v + 82.55) / 7.43))) ** 2))
    values['hssp_P'] = hssp_P
    jcass = (1.0 / (1.0 + np.exp(((v + 17.66945) / 3.21501))))
    values['jcass'] = jcass
    k2_i = kcaoff
    values['k2_i'] = k2_i
    k2_ss = kcaoff
    values['k2_ss'] = k2_ss
    k5_i = kcaoff
    values['k5_i'] = k5_i
    k5_ss = kcaoff
    values['k5_ss'] = k5_ss
    kCaSR = (MaxSR - ((MaxSR - MinSR) / (1 + ((ec50SR / casr) ** steepnessCaSR))))
    values['kCaSR'] = kCaSR
    kPKA_IKs = fIKs_PKA
    values['kPKA_IKs'] = kPKA_IKs
    k_uw = (0.026 * nu)
    values['k_uw'] = k_uw
    k_ws = (0.004 * (1 + (fMyBPC_PKA / 2)) * mu)
    values['k_ws'] = k_ws
    km2n = (jca * 0.84191)
    values['km2n'] = km2n
    ks = (Jrel_multiplier * 26.6)
    values['ks'] = ks
    lambda0 = min(lambda_max, lambda_)
    values['lambda0'] = lambda0
    mLss = (1 / (1 + np.exp(((-(v + 42.85)) / 5.264))))
    values['mLss'] = mLss
    mss = (1 / ((1 + np.exp(((-(v + 56.86)) / 9.03))) ** 2))
    values['mss'] = mss
    mss_P = (1 / ((1 + np.exp(((-(v + 61.86)) / 9.03))) ** 2))
    values['mss_P'] = mss_P
    nonlinearModifier = (0.2144 * np.exp((1.83 * casr)))
    values['nonlinearModifier'] = nonlinearModifier
    phosphorylationTotal = ((camk_f_PLB + fPLB_PKA) - (camk_f_PLB * fPLB_PKA))
    values['phosphorylationTotal'] = phosphorylationTotal
    sigmoidBaseCaI = (minCaI + ((maxCaI - minCaI) / (1 + ((ecCaI / caj) ** steepnessCaI))))
    values['sigmoidBaseCaI'] = sigmoidBaseCaI
    sigmoidTransition = (1 - (1 / (1 + ((1.86532 * caj) / 0.032))))
    values['sigmoidTransition'] = sigmoidTransition
    sigmoidTransition2 = (1 - (1 / (1 + ((1.86532 * casl) / 0.032))))
    values['sigmoidTransition2'] = sigmoidTransition2
    tau_cal = tau_ryr
    values['tau_cal'] = tau_cal
    tau_rel = max((bt / (1.0 + (0.0123 / casr))), 0.001)
    values['tau_rel'] = tau_rel
    tauxtof = ((8.5 * np.exp((-(((v + 45) / 50) ** 2)))) + 0.5)
    values['tauxtof'] = tauxtof
    tauxtos = ((9 / (1 + np.exp(((v + 3.0) / 15)))) + 0.5)
    values['tauxtos'] = tauxtos
    tauytof = ((85 * np.exp(((-((v + 40) ** 2)) / 220))) + 7)
    values['tauytof'] = tauytof
    tauytos = ((800 / (1 + np.exp(((v + 60.0) / 10)))) + 30)
    values['tauytos'] = tauytos
    td = (1.5 + (1 / (np.exp(((-0.05) * (v + 6))) + np.exp((0.09 * (v + 14))))))
    values['td'] = td
    tfcaf = (13.50673 + (1 / ((0.1542 * np.exp(((-(v - 1.31611)) / 11.3396))) + (0.1542 * np.exp(((v - 1.31611) / 11.3396))))))
    values['tfcaf'] = tfcaf
    tfcas = (177.95813 + (1 / (((4.73955e-4) * np.exp(((-(v - 0.79049)) / 0.81777))) + ((4.73955e-4) * np.exp(((v + 2.40474) / 1.90812))))))
    values['tfcas'] = tfcas
    tff = (6.17111 + (1 / ((0.00126 * np.exp(((-(v + 26.63596)) / 9.69961))) + (0.00126 * np.exp(((v + 26.63596) / 9.69961))))))
    values['tff'] = tff
    tfs = (2719.22489 + (1 / (((7.19411e-5) * np.exp(((-(v + 5.74631)) / 10.8769))) + ((7.19411e-5) * np.exp(((v + 5.74631) / 16.31535))))))
    values['tfs'] = tfs
    thLp = (3 * thL)
    values['thLp'] = thLp
    tm = ((0.1292 * np.exp((-(((v + 45.79) / 15.54) ** 2)))) + (0.06487 * np.exp((-(((v - 4.823) / 51.12) ** 2)))))
    values['tm'] = tm
    tmL = ((0.1292 * np.exp((-(((v + 45.79) / 15.54) ** 2)))) + (0.06487 * np.exp((-(((v - 4.823) / 51.12) ** 2)))))
    values['tmL'] = tmL
    vcell = ((1e-15) * np.pi * (cellRadius ** 2) * cellLength)
    values['vcell'] = vcell
    vffrt = ((v * F * F) / (R * T))
    values['vffrt'] = vffrt
    vfrt = ((v * F) / (R * T))
    values['vfrt'] = vfrt
    xkb = (1 / (1 + np.exp(((-(v - 10.8968)) / 23.9871))))
    values['xkb'] = xkb
    xtoss = (1 / (1 + np.exp(((-(v - 19.0)) / 13))))
    values['xtoss'] = xtoss
    xtoss_p = (1 / (1 + np.exp(((-(v - 29)) / 13))))
    values['xtoss_p'] = xtoss_p
    ytoss = (1 / (1 + np.exp(((v + 19.5) / 5))))
    values['ytoss'] = ytoss
    Afcas = (1 - Afcaf)
    values['Afcas'] = Afcas
    CaMK_active = (bound + camk_trap)
    values['CaMK_active'] = CaMK_active
    GK1 = (GK1_b * IK1__celltype_factor * IK1__sex_factor * IK1_multiplier)
    values['GK1'] = GK1
    GKr = (GKr_b * IKr__celltype_factor * IKr__sex_factor * IKr_multiplier)
    values['GKr'] = GKr
    GKs = (GKs_b * IKs__celltype_factor * IKs__sex_factor * IKs_multiplier)
    values['GKs'] = GKs
    GNaL = ((GNaL_b * INaL_multiplier * (1 + fINaLp) * 0.7) if (celltype == 1) else (GNaL_b * INaL_multiplier * (1 + fINaLp)))
    values['GNaL'] = GNaL
    GNa_P = (1.25 * GNa)
    values['GNa_P'] = GNa_P
    Gncx = (Gncx_b * INaCa__celltype_factor * INaCa__sex_factor * INaCa_multiplier)
    values['Gncx'] = Gncx
    ICab_junc = (Fjunc * GCab * (v - ECaj))
    values['ICab_junc'] = ICab_junc
    ICab_sl = (Fsl * GCab * (v - ECasl))
    values['ICab_sl'] = ICab_sl
    IClCa_junc = (((0.5 * GClCa) / (1 + (KdClCa / caj))) * (v - ECl))
    values['IClCa_junc'] = IClCa_junc
    IClCa_sl = (((0.5 * GClCa) / (1 + (KdClCa / casl))) * (v - ECl))
    values['IClCa_sl'] = IClCa_sl
    IClb = (GClb * (v - ECl))
    values['IClb'] = IClb
    IKb = (GKb * xkb * (v - EK))
    values['IKb'] = IKb
    IKr__alpha = (0.1161 * np.exp((0.2990 * vfrt)))
    values['IKr__alpha'] = IKr__alpha
    IKr__beta = (0.2442 * np.exp(((-1.604) * vfrt)))
    values['IKr__beta'] = IKr__beta
    INaBase_CaMK = (GNa * (m ** 3.0) * hp * jp)
    values['INaBase_CaMK'] = INaBase_CaMK
    INaBase_NP = (GNa * (m ** 3.0) * h * j)
    values['INaBase_NP'] = INaBase_NP
    INaK_junc_PKA = ((((Fjunc * IbarNaK * fnak) / (1 + ((KmNaip_PKA / naj) ** 4))) * ko) / (ko + KmKo))
    values['INaK_junc_PKA'] = INaK_junc_PKA
    INaK_junc_noPKA = ((((Fjunc * IbarNaK * fnak) / (1 + ((KmNaip / naj) ** 4))) * ko) / (ko + KmKo))
    values['INaK_junc_noPKA'] = INaK_junc_noPKA
    INaK_sl_PKA = ((((Fsl * IbarNaK * fnak) / (1 + ((KmNaip_PKA / nasl) ** 4))) * ko) / (ko + KmKo))
    values['INaK_sl_PKA'] = INaK_sl_PKA
    INaK_sl_noPKA = ((((Fsl * IbarNaK * fnak) / (1 + ((KmNaip / nasl) ** 4))) * ko) / (ko + KmKo))
    values['INaK_sl_noPKA'] = INaK_sl_noPKA
    INabj = (Fjunc * GNab * (v - ENaj))
    values['INabj'] = INabj
    INabsl = (Fsl * GNab * (v - ENasl))
    values['INabsl'] = INabsl
    IpCa_junc = ((Fjunc * (Q10SLCaP ** Qpow) * IbarSLCaP * (caj ** 1.6)) / ((KmPCa ** 1.6) + (caj ** 1.6)))
    values['IpCa_junc'] = IpCa_junc
    IpCa_sl = ((Fsl * (Q10SLCaP ** Qpow) * IbarSLCaP * (casl ** 1.6)) / ((KmPCa ** 1.6) + (casl ** 1.6)))
    values['IpCa_sl'] = IpCa_sl
    Itof = (Gto_fast * (v - EK) * (((1 - fItop) * xtof * ytof) + (fItop * xtof_p * ytof_p)))
    values['Itof'] = Itof
    Itos = (Gto_slow * (v - EK) * (((1 - fItop) * xtos * ytos) + (fItop * xtos_p * ytos_p)))
    values['Itos'] = Itos
    J_SRCarel_np = ((ks * ryr_O * (casr - caj)) + Jrel_ICaLdep)
    values['J_SRCarel_np'] = J_SRCarel_np
    J_SRCarel_p = ((ks * ryr_O_p * (casr - caj)) + Jrel_ICaLdep)
    values['J_SRCarel_p'] = J_SRCarel_p
    J_SRleak = ((1.59306e-6) * (casr - caj) * nonlinearModifier * CaMKIILeakMultiplier)
    values['J_SRleak'] = J_SRleak
    Lfac = max(0, (1 + (beta_0 * ((lambda0 + min(lambda_min, lambda0)) - (1 + lambda_min)))))
    values['Lfac'] = Lfac
    PCaK = ((2.211399546628710e-04) * PCa)
    values['PCaK'] = PCaK
    PCaNa = ((7.734329695819495e-04) * PCa)
    values['PCaNa'] = PCaNa
    PCa_Pb = (1.9 * PCa)
    values['PCa_Pb'] = PCa_Pb
    PCap = (1.1 * PCa)
    values['PCap'] = PCap
    P_g_0 = (gKs_factor * (0.2 + (0.2 * kPKA_IKs)))
    values['P_g_0'] = P_g_0
    P_g_max = (gKs_factor * (0.8 + (7 * kPKA_IKs)))
    values['P_g_max'] = P_g_max
    P_tau_0 = (26 + (9 * kPKA_IKs))
    values['P_tau_0'] = P_tau_0
    P_tau_max = (40 + (4 * kPKA_IKs))
    values['P_tau_max'] = P_tau_max
    P_vh_0 = ((-1) - (10 * kPKA_IKs))
    values['P_vh_0'] = P_vh_0
    P_vh_max = ((-12) - (9 * kPKA_IKs))
    values['P_vh_max'] = P_vh_max
    RI_to_CI = (baseRateCaI * sigmoidBaseCaI)
    values['RI_to_CI'] = RI_to_CI
    aK1 = (4.094 / (1 + np.exp((0.1217 * ((v - EK) - 49.934)))))
    values['aK1'] = aK1
    alpha_2 = (0.0578 * np.exp((0.9710 * vfrt)))
    values['alpha_2'] = alpha_2
    alpha_C2ToI = ((0.52e-4) * np.exp((1.525 * vfrt)))
    values['alpha_C2ToI'] = alpha_C2ToI
    alpha_i = (0.2533 * np.exp((0.5953 * vfrt)))
    values['alpha_i'] = alpha_i
    anca = (1 / ((k2n / km2n) + ((1 + (Kmn / caj)) ** 3.80763)))
    values['anca'] = anca
    anca_i = (1 / ((k2n / km2n) + ((1 + (Kmn / casl)) ** 3.80763)))
    values['anca_i'] = anca_i
    bK1 = (((15.72 * np.exp((0.0674 * ((v - EK) - 3.257)))) + np.exp((0.0618 * ((v - EK) - 594.31)))) / (1 + np.exp(((-0.1629) * ((v - EK) + 14.207)))))
    values['bK1'] = bK1
    beta_2 = ((0.349e-3) * np.exp(((-1.062) * vfrt)))
    values['beta_2'] = beta_2
    beta_i = (0.04568 * np.exp(((-0.8209) * vfrt)))
    values['beta_i'] = beta_i
    ca50 = ((fPKA_TnI * 0.7645) + (contraction__beta_1 * min(0.2, (lambda_ - 1))))
    values['ca50'] = ca50
    casig_SERCA_act = (bound_serca + casig_serca_trap)
    values['casig_SERCA_act'] = casig_SERCA_act
    cds = ((phi * k_ws * wfrac * (1 - dr)) / dr)
    values['cds'] = cds
    cdw = ((phi * k_uw * (1 - wfrac)) / wfrac)
    values['cdw'] = cdw
    d_buffers_CaM = ((kon_cam * cai * (Bmax_CaM - buffers_CaM)) - (koff_cam * buffers_CaM))
    values['d_buffers_CaM'] = d_buffers_CaM
    f = ((Aff * ff) + (Afs * fs))
    values['f'] = f
    fBP = ((Aff * fBPf) + (Afs * fs_P))
    values['fBP'] = fBP
    fBPss = fss_P
    values['fBPss'] = fBPss
    fICaL_P = fICaLP
    values['fICaL_P'] = fICaL_P
    fINa_BP = (fINap * fINa_P)
    values['fINa_BP'] = fINa_BP
    f_P = ((Aff * ff_P) + (Afs * fs_P))
    values['f_P'] = f_P
    fcass = fss
    values['fcass'] = fcass
    fcass_P = fss_P
    values['fcass_P'] = fcass_P
    fp = ((Aff * ffp) + (Afs * fs))
    values['fp'] = fp
    gamma_cai = (10 ** ((-constA) * 4 * ((np.sqrt(Ii) / (1 + np.sqrt(Ii))) - (0.3 * Ii))))
    values['gamma_cai'] = gamma_cai
    gamma_cao = (10 ** ((-constA) * 4 * ((np.sqrt(Io) / (1 + np.sqrt(Io))) - (0.3 * Io))))
    values['gamma_cao'] = gamma_cao
    gamma_cass = (10 ** ((-constA) * 4 * ((np.sqrt(Iss) / (1 + np.sqrt(Iss))) - (0.3 * Iss))))
    values['gamma_cass'] = gamma_cass
    gamma_ki = (10 ** ((-constA) * 1 * ((np.sqrt(Ii) / (1 + np.sqrt(Ii))) - (0.3 * Ii))))
    values['gamma_ki'] = gamma_ki
    gamma_ko = (10 ** ((-constA) * 1 * ((np.sqrt(Io) / (1 + np.sqrt(Io))) - (0.3 * Io))))
    values['gamma_ko'] = gamma_ko
    gamma_kss = (10 ** ((-constA) * 1 * ((np.sqrt(Iss) / (1 + np.sqrt(Iss))) - (0.3 * Iss))))
    values['gamma_kss'] = gamma_kss
    gamma_nai = (10 ** ((-constA) * 1 * ((np.sqrt(Ii) / (1 + np.sqrt(Ii))) - (0.3 * Ii))))
    values['gamma_nai'] = gamma_nai
    gamma_nao = (10 ** ((-constA) * 1 * ((np.sqrt(Io) / (1 + np.sqrt(Io))) - (0.3 * Io))))
    values['gamma_nao'] = gamma_nao
    gamma_nass = (10 ** ((-constA) * 1 * ((np.sqrt(Iss) / (1 + np.sqrt(Iss))) - (0.3 * Iss))))
    values['gamma_nass'] = gamma_nass
    h11_i = ((nao * nao) / (h10_i * kna1 * kna2))
    values['h11_i'] = h11_i
    h11_ss = ((nao * nao) / (h10_ss * kna1 * kna2))
    values['h11_ss'] = h11_ss
    h12_i = (1 / h10_i)
    values['h12_i'] = h12_i
    h12_ss = (1 / h10_ss)
    values['h12_ss'] = h12_ss
    h1_i = (1 + ((nai / kna3) * (1 + hna)))
    values['h1_i'] = h1_i
    h1_ss = (1 + ((nasl / kna3) * (1 + hna)))
    values['h1_ss'] = h1_ss
    h5_i = ((nai * nai) / (h4_i * kna1 * kna2))
    values['h5_i'] = h5_i
    h5_ss = ((nasl * nasl) / (h4_ss * kna1 * kna2))
    values['h5_ss'] = h5_ss
    h6_i = (1 / h4_i)
    values['h6_i'] = h6_i
    h6_ss = (1 / h4_ss)
    values['h6_ss'] = h6_ss
    h7_i = (1 + ((nao / kna3) * (1 + (1 / hna))))
    values['h7_i'] = h7_i
    h7_ss = (1 + ((nao / kna3) * (1 + (1 / hna))))
    values['h7_ss'] = h7_ss
    jss = hss
    values['jss'] = jss
    jss_P = hss_P
    values['jss_P'] = jss_P
    jssp_P = hssp_P
    values['jssp_P'] = jssp_P
    k_su = (k_ws * ((1 / dr) - 1) * wfrac)
    values['k_su'] = k_su
    k_wu = ((k_uw * ((1 / wfrac) - 1)) - k_ws)
    values['k_wu'] = k_wu
    kiSRCa = (kiCa * kCaSR)
    values['kiSRCa'] = kiSRCa
    koSRCa = (koCa / kCaSR)
    values['koSRCa'] = koSRCa
    ktm_block = ((ktm_unblock * (perm50 ** nperm) * 0.5) / ((0.5 - XSSS) - XWSS))
    values['ktm_block'] = ktm_block
    tauTransition = (1.0967 + ((1 - sigmoidTransition) * 141.4299))
    values['tauTransition'] = tauTransition
    tauTransition2 = (1.0967 + ((1 - sigmoidTransition2) * 141.4299))
    values['tauTransition2'] = tauTransition2
    tauytof_p = (tauytof * dti_develop * dti_recover)
    values['tauytof_p'] = tauytof_p
    tauytos_p = (tauytos * dti_develop * dti_recover)
    values['tauytos_p'] = tauytos_p
    tfcafp = (2.5 * tfcaf)
    values['tfcafp'] = tfcafp
    tffp = (2.5 * tff)
    values['tffp'] = tffp
    th = (1 / (ah + bh))
    values['th'] = th
    tj = (1 / (aj + bj))
    values['tj'] = tj
    vjunc = ((5.39e-4) * vcell)
    values['vjunc'] = vjunc
    vmyo = (0.65 * vcell)
    values['vmyo'] = vmyo
    vsl = (0.02 * vcell)
    values['vsl'] = vsl
    vsr = (0.035 * vcell)
    values['vsr'] = vsr
    xb_su_gamma = (gamma_rate * contraction_XS)
    values['xb_su_gamma'] = xb_su_gamma
    xb_uw = (k_uw * XU * (1.0 + (fMyBPC_PKA / 2)))
    values['xb_uw'] = xb_uw
    xb_ws = (k_ws * contraction_XW)
    values['xb_ws'] = xb_ws
    xb_wu_gamma = (gamma_rate_w * contraction_XW)
    values['xb_wu_gamma'] = xb_wu_gamma
    Bmax_Csqn = (((1.3655214e-1) * vmyo) / vsr)
    values['Bmax_Csqn'] = Bmax_Csqn
    Bmax_SLhighj = (((1.75755e-4) * vmyo) / vjunc)
    values['Bmax_SLhighj'] = Bmax_SLhighj
    Bmax_SLhighsl = (((1.215423e-2) * vmyo) / vsl)
    values['Bmax_SLhighsl'] = Bmax_SLhighsl
    Bmax_SLlowj = (((4.89983e-4) * vmyo) / vjunc)
    values['Bmax_SLlowj'] = Bmax_SLlowj
    Bmax_SLlowsl = (((3.3923e-2) * vmyo) / vsl)
    values['Bmax_SLlowsl'] = Bmax_SLlowsl
    CaMK_Phos_ss_ICaL = (CaMK_active / (CaMK_active + K_Phos_CaMK))
    values['CaMK_Phos_ss_ICaL'] = CaMK_Phos_ss_ICaL
    CaMK_Phos_ss_PLB = (CaMK_active / (CaMK_active + 10))
    values['CaMK_Phos_ss_PLB'] = CaMK_Phos_ss_PLB
    CaMK_Phos_ss_RyR = (CaMK_active / (CaMK_active + 1))
    values['CaMK_Phos_ss_RyR'] = CaMK_Phos_ss_RyR
    ICab = (ICab_junc + ICab_sl)
    values['ICab'] = ICab
    IClCa = (IClCa_junc + IClCa_sl)
    values['IClCa'] = IClCa
    IKr = (GKr * np.sqrt((ko / 5)) * O * (v - EK))
    values['IKr'] = IKr
    INaBase_BP = (GNa_P * (m_P ** 3.0) * hp_P * jp_P)
    values['INaBase_BP'] = INaBase_BP
    INaBase_PKA = (GNa_P * (m_P ** 3.0) * h_P * j_P)
    values['INaBase_PKA'] = INaBase_PKA
    INaKj = (((1 - fINaK_PKA) * INaK_junc_noPKA) + (fINaK_PKA * INaK_junc_PKA))
    values['INaKj'] = INaKj
    INaKsl = (((1 - fINaK_PKA) * INaK_sl_noPKA) + (fINaK_PKA * INaK_sl_PKA))
    values['INaKsl'] = INaKsl
    INaLj = (Fjunc * GNaL * (v - ENaj) * mL * (((1 - fINaLp) * hL) + (fINaLp * hLp)))
    values['INaLj'] = INaLj
    INaLsl = (Fsl * GNaL * (v - ENasl) * mL * (((1 - fINaLp) * hL) + (fINaLp * hLp)))
    values['INaLsl'] = INaLsl
    INab = (INabj + INabsl)
    values['INab'] = INab
    IpCa = (IpCa_junc + IpCa_sl)
    values['IpCa'] = IpCa
    Ito = (Itos + Itof)
    values['Ito'] = Ito
    J_SRCarel = ((J_SRCarel_p * camk_f_RyR) + (J_SRCarel_np * (1 - camk_f_RyR)))
    values['J_SRCarel'] = J_SRCarel
    K1ss = (aK1 / (aK1 + bK1))
    values['K1ss'] = K1ss
    PCaKp = ((2.211399546628710e-04) * PCap)
    values['PCaKp'] = PCaKp
    PCaNap = ((7.734329695819495e-04) * PCap)
    values['PCaNap'] = PCaNap
    PCa_P = ((PCa_Pb * 1.025) if (celltype == 1) else ((PCa_Pb * 1.1) if (celltype == 2) else PCa_Pb))
    values['PCa_P'] = PCa_P
    PhiCaK_i = ((1 * vffrt * ((gamma_ki * ki * np.exp((1 * vfrt))) - (gamma_ko * ko))) / (np.exp((1 * vfrt)) - 1))
    values['PhiCaK_i'] = PhiCaK_i
    PhiCaK_ss = ((1 * vffrt * ((gamma_kss * ki * np.exp((1 * vfrt))) - (gamma_ko * ko))) / (np.exp((1 * vfrt)) - 1))
    values['PhiCaK_ss'] = PhiCaK_ss
    PhiCaL_i = ((4 * vffrt * ((gamma_cai * cai * np.exp((2 * vfrt))) - (gamma_cao * cao))) / (np.exp((2 * vfrt)) - 1))
    values['PhiCaL_i'] = PhiCaL_i
    PhiCaL_ss = ((4 * vffrt * ((gamma_cass * casl * np.exp((2 * vfrt))) - (gamma_cao * cao))) / (np.exp((2 * vfrt)) - 1))
    values['PhiCaL_ss'] = PhiCaL_ss
    PhiCaNa_i = ((1 * vffrt * ((gamma_nai * nai * np.exp((1 * vfrt))) - (gamma_nao * nao))) / (np.exp((1 * vfrt)) - 1))
    values['PhiCaNa_i'] = PhiCaNa_i
    PhiCaNa_ss = ((1 * vffrt * ((gamma_nass * nasl * np.exp((1 * vfrt))) - (gamma_nao * nao))) / (np.exp((1 * vfrt)) - 1))
    values['PhiCaNa_ss'] = PhiCaNa_ss
    Ta = (((Lfac * PKAForceMultiplier * Tref) / dr) * (((contraction_ZETAS + 1) * contraction_XS) + (contraction_ZETAW * contraction_XW)))
    values['Ta'] = Ta
    Vmax_mult = (1 + (Max_Vmax_SERCA_Ca / (1 + ((Km_SERCA_Ca / casig_SERCA_act) ** 2))))
    values['Vmax_mult'] = Vmax_mult
    VsTs_Ca_junc = (P_tau_0 + ((P_tau_max - P_tau_0) / (1 + (((150e-6) / caj) ** 3))))
    values['VsTs_Ca_junc'] = VsTs_Ca_junc
    VsTs_Ca_sl = (P_tau_0 + ((P_tau_max - P_tau_0) / (1 + (((150e-6) / casl) ** 3))))
    values['VsTs_Ca_sl'] = VsTs_Ca_sl
    VsXs_Ca_junc = (P_vh_0 + ((P_vh_max - P_vh_0) / (1 + (((350e-6) / caj) ** 4))))
    values['VsXs_Ca_junc'] = VsXs_Ca_junc
    VsXs_Ca_sl = (P_vh_0 + ((P_vh_max - P_vh_0) / (1 + (((350e-6) / casl) ** 4))))
    values['VsXs_Ca_sl'] = VsXs_Ca_sl
    beta_ItoC2 = ((beta_2 * beta_i * alpha_C2ToI) / (alpha_2 * alpha_i))
    values['beta_ItoC2'] = beta_ItoC2
    d_contraction_Ca_TRPN = (koff * (((((cai * 1000) / ca50) ** TRPN_n) * (1 - contraction_Ca_TRPN)) - contraction_Ca_TRPN))
    values['d_contraction_Ca_TRPN'] = d_contraction_Ca_TRPN
    fICaL_BP = (fICaLp * fICaL_P)
    values['fICaL_BP'] = fICaL_BP
    fINa_CaMKonly = (fINap - fINa_BP)
    values['fINa_CaMKonly'] = fINa_CaMKonly
    fINa_PKAonly = (fINa_P - fINa_BP)
    values['fINa_PKAonly'] = fINa_PKAonly
    fca = ((Afcaf * fcaf) + (Afcas * fcas))
    values['fca'] = fca
    fcaBP = ((Afcaf * fcaBPf) + (Afcas * fcas_P))
    values['fcaBP'] = fcaBP
    fcaBPss = fcass_P
    values['fcaBPss'] = fcaBPss
    fcap = ((Afcaf * fcafp) + (Afcas * fcas))
    values['fcap'] = fcap
    fcap_P = ((Afcaf * fcaf_P) + (Afcas * fcas_P))
    values['fcap_P'] = fcap_P
    gks_junc = (P_g_0 + ((P_g_max - P_g_0) / (1 + (((150e-6) / caj) ** 1.3))))
    values['gks_junc'] = gks_junc
    gks_sl = (P_g_0 + ((P_g_max - P_g_0) / (1 + (((150e-6) / casl) ** 1.3))))
    values['gks_sl'] = gks_sl
    h2_i = ((nai * hna) / (kna3 * h1_i))
    values['h2_i'] = h2_i
    h2_ss = ((nasl * hna) / (kna3 * h1_ss))
    values['h2_ss'] = h2_ss
    h3_i = (1 / h1_i)
    values['h3_i'] = h3_i
    h3_ss = (1 / h1_ss)
    values['h3_ss'] = h3_ss
    h8_i = (nao / (kna3 * hna * h7_i))
    values['h8_i'] = h8_i
    h8_ss = (nao / (kna3 * hna * h7_ss))
    values['h8_ss'] = h8_ss
    h9_i = (1 / h7_i)
    values['h9_i'] = h9_i
    h9_ss = (1 / h7_ss)
    values['h9_ss'] = h9_ss
    k1_i = (h12_i * cao * kcaon)
    values['k1_i'] = k1_i
    k1_ss = (h12_ss * cao * kcaon)
    values['k1_ss'] = k1_ss
    k6_i = (h6_i * cai * kcaon)
    values['k6_i'] = k6_i
    k6_ss = (h6_ss * casl * kcaon)
    values['k6_ss'] = k6_ss
    tjp = (1.46 * tj)
    values['tjp'] = tjp
    xb_su = (k_su * contraction_XS)
    values['xb_su'] = xb_su
    xb_wu = (k_wu * contraction_XW)
    values['xb_wu'] = xb_wu
    ICaK_i_CaMK = (PCaKp * PhiCaK_i * d * ((fp * (1.0 - nca_i)) + (jca * fcap * nca_i)))
    values['ICaK_i_CaMK'] = ICaK_i_CaMK
    ICaK_i_NP = (PCaK * PhiCaK_i * d * ((f * (1.0 - nca_i)) + (jca * fca * nca_i)))
    values['ICaK_i_NP'] = ICaK_i_NP
    ICaK_ss_CaMK = (PCaKp * PhiCaK_ss * d * ((fp * (1.0 - nca)) + (jca * fcap * nca)))
    values['ICaK_ss_CaMK'] = ICaK_ss_CaMK
    ICaK_ss_NP = (PCaK * PhiCaK_ss * d * ((f * (1.0 - nca)) + (jca * fca * nca)))
    values['ICaK_ss_NP'] = ICaK_ss_NP
    ICaL_i_BP = (PCa_P * PhiCaL_i * d_P * ((fBP * (1.0 - nca_i)) + (jca * fcaBP * nca_i)))
    values['ICaL_i_BP'] = ICaL_i_BP
    ICaL_i_CaMK = (PCap * PhiCaL_i * d * ((fp * (1.0 - nca_i)) + (jca * fcap * nca_i)))
    values['ICaL_i_CaMK'] = ICaL_i_CaMK
    ICaL_i_NP = (PCa * PhiCaL_i * d * ((f * (1.0 - nca_i)) + (jca * fca * nca_i)))
    values['ICaL_i_NP'] = ICaL_i_NP
    ICaL_i_PKA = (PCa_P * PhiCaL_i * d_P * ((f_P * (1.0 - nca_i)) + (jca * fcap_P * nca_i)))
    values['ICaL_i_PKA'] = ICaL_i_PKA
    ICaL_ss_BP = (PCa_P * PhiCaL_ss * d_P * ((fBP * (1.0 - nca)) + (jca * fcaBP * nca)))
    values['ICaL_ss_BP'] = ICaL_ss_BP
    ICaL_ss_CaMK = (PCap * PhiCaL_ss * d * ((fp * (1.0 - nca)) + (jca * fcap * nca)))
    values['ICaL_ss_CaMK'] = ICaL_ss_CaMK
    ICaL_ss_NP = (PCa * PhiCaL_ss * d * ((f * (1.0 - nca)) + (jca * fca * nca)))
    values['ICaL_ss_NP'] = ICaL_ss_NP
    ICaL_ss_PKA = (PCa_P * PhiCaL_ss * d_P * ((f_P * (1.0 - nca)) + (jca * fcap_P * nca)))
    values['ICaL_ss_PKA'] = ICaL_ss_PKA
    ICaNa_i_CaMK = (PCaNap * PhiCaNa_i * d * ((fp * (1.0 - nca_i)) + (jca * fcap * nca_i)))
    values['ICaNa_i_CaMK'] = ICaNa_i_CaMK
    ICaNa_i_NP = (PCaNa * PhiCaNa_i * d * ((f * (1.0 - nca_i)) + (jca * fca * nca_i)))
    values['ICaNa_i_NP'] = ICaNa_i_NP
    ICaNa_ss_CaMK = (PCaNap * PhiCaNa_ss * d * ((fp * (1.0 - nca)) + (jca * fcap * nca)))
    values['ICaNa_ss_CaMK'] = ICaNa_ss_CaMK
    ICaNa_ss_NP = (PCaNa * PhiCaNa_ss * d * ((f * (1.0 - nca)) + (jca * fca * nca)))
    values['ICaNa_ss_NP'] = ICaNa_ss_NP
    ICl_tot = (IClCa + IClb)
    values['ICl_tot'] = ICl_tot
    IK1 = (GK1 * np.sqrt((ko / 5)) * K1ss * (v - EK))
    values['IK1'] = IK1
    IKs_junc = (Fjunc * GKs * gks_junc * (xs_junc ** 2) * (v - EKs))
    values['IKs_junc'] = IKs_junc
    IKs_sl = (Fsl * GKs * gks_sl * (xs_sl ** 2) * (v - EKs))
    values['IKs_sl'] = IKs_sl
    INaBase = (((((1 - fINa_CaMKonly) - fINa_PKAonly) - fINa_BP) * INaBase_NP) + (fINa_CaMKonly * INaBase_CaMK) + (fINa_PKAonly * INaBase_PKA) + (fINa_BP * INaBase_BP))
    values['INaBase'] = INaBase
    INaK = (INaKj + INaKsl)
    values['INaK'] = INaK
    INaL = (INaLj + INaLsl)
    values['INaL'] = INaL
    Jserca_np = (((Q10SRCaP ** Qpow) * Vmax_SRCaP * Vmax_mult * (((cai / Kmf) ** hillSRCaP) - ((casr / Kmr) ** hillSRCaP))) / (1 + ((cai / Kmf) ** hillSRCaP) + ((casr / Kmr) ** hillSRCaP)))
    values['Jserca_np'] = Jserca_np
    Jserca_p = (((Q10SRCaP ** Qpow) * Vmax_SRCaP * Vmax_mult * (((cai / Kmf_p) ** hillSRCaP) - ((casr / Kmr) ** hillSRCaP))) / (1 + ((cai / Kmf_p) ** hillSRCaP) + ((casr / Kmr) ** hillSRCaP)))
    values['Jserca_p'] = Jserca_p
    PCaK_P = ((2.211399546628710e-04) * PCa_P)
    values['PCaK_P'] = PCaK_P
    PCaNa_P = ((7.734329695819495e-04) * PCa_P)
    values['PCaNa_P'] = PCaNa_P
    d_buffers_Csqn = ((kon_csqn * casr * (Bmax_Csqn - buffers_Csqn)) - (koff_csqn * buffers_Csqn))
    values['d_buffers_Csqn'] = d_buffers_Csqn
    d_buffers_SLHj = ((kon_slh * caj * (Bmax_SLhighj - buffers_SLHj)) - (koff_slh * buffers_SLHj))
    values['d_buffers_SLHj'] = d_buffers_SLHj
    d_buffers_SLHsl = ((kon_slh * casl * (Bmax_SLhighsl - buffers_SLHsl)) - (koff_slh * buffers_SLHsl))
    values['d_buffers_SLHsl'] = d_buffers_SLHsl
    d_buffers_SLLj = ((kon_sll * caj * (Bmax_SLlowj - buffers_SLLj)) - (koff_sll * buffers_SLLj))
    values['d_buffers_SLLj'] = d_buffers_SLLj
    d_buffers_SLLsl = ((kon_sll * casl * (Bmax_SLlowsl - buffers_SLLsl)) - (koff_sll * buffers_SLLsl))
    values['d_buffers_SLLsl'] = d_buffers_SLLsl
    d_buffers_TnClow = (Bmax_TnClow * d_contraction_Ca_TRPN)
    values['d_buffers_TnClow'] = d_buffers_TnClow
    fICaL_CaMKonly = (fICaLp - fICaL_BP)
    values['fICaL_CaMKonly'] = fICaL_CaMKonly
    fICaL_PKAonly = (fICaL_P - fICaL_BP)
    values['fICaL_PKAonly'] = fICaL_PKAonly
    k3p_i = (h9_i * wca)
    values['k3p_i'] = k3p_i
    k3p_ss = (h9_ss * wca)
    values['k3p_ss'] = k3p_ss
    k3pp_i = (h8_i * wnaca)
    values['k3pp_i'] = k3pp_i
    k3pp_ss = (h8_ss * wnaca)
    values['k3pp_ss'] = k3pp_ss
    k4p_i = ((h3_i * wca) / hca)
    values['k4p_i'] = k4p_i
    k4p_ss = ((h3_ss * wca) / hca)
    values['k4p_ss'] = k4p_ss
    k4pp_i = (h2_i * wnaca)
    values['k4pp_i'] = k4pp_i
    k4pp_ss = (h2_ss * wnaca)
    values['k4pp_ss'] = k4pp_ss
    k7_i = (h5_i * h2_i * wna)
    values['k7_i'] = k7_i
    k7_ss = (h5_ss * h2_ss * wna)
    values['k7_ss'] = k7_ss
    k8_i = (h8_i * h11_i * wna)
    values['k8_i'] = k8_i
    k8_ss = (h8_ss * h11_ss * wna)
    values['k8_ss'] = k8_ss
    tauxs_junc = (2 * (50 + ((50 + (350 * np.exp(((-((v + 30) ** 2)) / 4000)))) / (1 + np.exp(((-(v + VsTs_Ca_junc)) / 10))))))
    values['tauxs_junc'] = tauxs_junc
    tauxs_sl = (2 * (50 + ((50 + (350 * np.exp(((-((v + 30) ** 2)) / 4000)))) / (1 + np.exp(((-(v + VsTs_Ca_sl)) / 10))))))
    values['tauxs_sl'] = tauxs_sl
    xsss_junc = (1 / (1 + np.exp(((-(v - VsXs_Ca_junc)) / 25))))
    values['xsss_junc'] = xsss_junc
    xsss_sl = (1 / (1 + np.exp(((-(v - VsXs_Ca_sl)) / 25))))
    values['xsss_sl'] = xsss_sl
    ICaK_i_BP = (PCaK_P * PhiCaK_i * d_P * ((fBP * (1.0 - nca_i)) + (jca * fcaBP * nca_i)))
    values['ICaK_i_BP'] = ICaK_i_BP
    ICaK_i_PKA = (PCaK_P * PhiCaK_i * d_P * ((f_P * (1.0 - nca_i)) + (jca * fcap_P * nca_i)))
    values['ICaK_i_PKA'] = ICaK_i_PKA
    ICaK_ss_BP = (PCaK_P * PhiCaK_ss * d_P * ((fBP * (1.0 - nca)) + (jca * fcaBP * nca)))
    values['ICaK_ss_BP'] = ICaK_ss_BP
    ICaK_ss_PKA = (PCaK_P * PhiCaK_ss * d_P * ((f_P * (1.0 - nca)) + (jca * fcap_P * nca)))
    values['ICaK_ss_PKA'] = ICaK_ss_PKA
    ICaL_junc = ((((((1 - fICaL_CaMKonly) - fICaL_PKAonly) - fICaL_BP) * ICaL_ss_NP) + (fICaL_CaMKonly * ICaL_ss_CaMK) + (fICaL_PKAonly * ICaL_ss_PKA) + (fICaL_BP * ICaL_ss_BP)) * ical_pureCDI_junc * ICaL_fractionSS)
    values['ICaL_junc'] = ICaL_junc
    ICaL_sl = ((((((1 - fICaL_CaMKonly) - fICaL_PKAonly) - fICaL_BP) * ICaL_i_NP) + (fICaL_CaMKonly * ICaL_i_CaMK) + (fICaL_PKAonly * ICaL_i_PKA) + (fICaL_BP * ICaL_i_BP)) * ical_pureCDI_sl * (1 - ICaL_fractionSS))
    values['ICaL_sl'] = ICaL_sl
    ICaNa_i_BP = (PCaNa_P * PhiCaNa_i * d_P * ((fBP * (1.0 - nca_i)) + (jca * fcaBP * nca_i)))
    values['ICaNa_i_BP'] = ICaNa_i_BP
    ICaNa_i_PKA = (PCaNa_P * PhiCaNa_i * d_P * ((f_P * (1.0 - nca_i)) + (jca * fcap_P * nca_i)))
    values['ICaNa_i_PKA'] = ICaNa_i_PKA
    ICaNa_ss_BP = (PCaNa_P * PhiCaNa_ss * d_P * ((fBP * (1.0 - nca)) + (jca * fcaBP * nca)))
    values['ICaNa_ss_BP'] = ICaNa_ss_BP
    ICaNa_ss_PKA = (PCaNa_P * PhiCaNa_ss * d_P * ((f_P * (1.0 - nca)) + (jca * fcap_P * nca)))
    values['ICaNa_ss_PKA'] = ICaNa_ss_PKA
    IKs = (IKs_junc + IKs_sl)
    values['IKs'] = IKs
    INaj = (Fjunc * INaBase * (v - ENaj))
    values['INaj'] = INaj
    INasl = (Fsl * INaBase * (v - ENasl))
    values['INasl'] = INasl
    J_CaB_cytosol = (d_buffers_TnClow + d_buffers_TnCHc + d_buffers_CaM + d_buffers_Myosin_ca + d_buffers_SRB)
    values['J_CaB_cytosol'] = J_CaB_cytosol
    J_CaB_junction = (d_buffers_SLLj + d_buffers_SLHj)
    values['J_CaB_junction'] = J_CaB_junction
    J_CaB_sl = (d_buffers_SLLsl + d_buffers_SLHsl)
    values['J_CaB_sl'] = J_CaB_sl
    Jserca = ((Jserca_np * (1 - phosphorylationTotal)) + (Jserca_p * phosphorylationTotal))
    values['Jserca'] = Jserca
    k3_i = (k3p_i + k3pp_i)
    values['k3_i'] = k3_i
    k3_ss = (k3p_ss + k3pp_ss)
    values['k3_ss'] = k3_ss
    k4_i = (k4p_i + k4pp_i)
    values['k4_i'] = k4_i
    k4_ss = (k4p_ss + k4pp_ss)
    values['k4_ss'] = k4_ss
    ICaK_junc = ((((((1 - fICaL_CaMKonly) - fICaL_PKAonly) - fICaL_BP) * ICaK_ss_NP) + (fICaL_CaMKonly * ICaK_ss_CaMK) + (fICaL_PKAonly * ICaK_ss_PKA) + (fICaL_BP * ICaK_ss_BP)) * ical_pureCDI_junc * ICaL_fractionSS)
    values['ICaK_junc'] = ICaK_junc
    ICaK_sl = ((((((1 - fICaL_CaMKonly) - fICaL_PKAonly) - fICaL_BP) * ICaK_i_NP) + (fICaL_CaMKonly * ICaK_i_CaMK) + (fICaL_PKAonly * ICaK_i_PKA) + (fICaL_BP * ICaK_i_BP)) * ical_pureCDI_sl * (1 - ICaL_fractionSS))
    values['ICaK_sl'] = ICaK_sl
    ICaL = (ICaL_junc + ICaL_sl)
    values['ICaL'] = ICaL
    ICaL_junc_positive = np.abs(ICaL_junc)
    values['ICaL_junc_positive'] = ICaL_junc_positive
    ICaNa_junc = ((((((1 - fICaL_CaMKonly) - fICaL_PKAonly) - fICaL_BP) * ICaNa_ss_NP) + (fICaL_CaMKonly * ICaNa_ss_CaMK) + (fICaL_PKAonly * ICaNa_ss_PKA) + (fICaL_BP * ICaNa_ss_BP)) * ical_pureCDI_junc * ICaL_fractionSS)
    values['ICaNa_junc'] = ICaNa_junc
    ICaNa_sl = ((((((1 - fICaL_CaMKonly) - fICaL_PKAonly) - fICaL_BP) * ICaNa_i_NP) + (fICaL_CaMKonly * ICaNa_i_CaMK) + (fICaL_PKAonly * ICaNa_i_PKA) + (fICaL_BP * ICaNa_i_BP)) * ical_pureCDI_sl * (1 - ICaL_fractionSS))
    values['ICaNa_sl'] = ICaNa_sl
    INa = (INaj + INasl)
    values['INa'] = INa
    x1_i = ((k2_i * k4_i * (k7_i + k6_i)) + (k5_i * k7_i * (k2_i + k3_i)))
    values['x1_i'] = x1_i
    x1_ss = ((k2_ss * k4_ss * (k7_ss + k6_ss)) + (k5_ss * k7_ss * (k2_ss + k3_ss)))
    values['x1_ss'] = x1_ss
    x2_i = ((k1_i * k7_i * (k4_i + k5_i)) + (k4_i * k6_i * (k1_i + k8_i)))
    values['x2_i'] = x2_i
    x2_ss = ((k1_ss * k7_ss * (k4_ss + k5_ss)) + (k4_ss * k6_ss * (k1_ss + k8_ss)))
    values['x2_ss'] = x2_ss
    x3_i = ((k1_i * k3_i * (k7_i + k6_i)) + (k8_i * k6_i * (k2_i + k3_i)))
    values['x3_i'] = x3_i
    x3_ss = ((k1_ss * k3_ss * (k7_ss + k6_ss)) + (k8_ss * k6_ss * (k2_ss + k3_ss)))
    values['x3_ss'] = x3_ss
    x4_i = ((k2_i * k8_i * (k4_i + k5_i)) + (k3_i * k5_i * (k1_i + k8_i)))
    values['x4_i'] = x4_i
    x4_ss = ((k2_ss * k8_ss * (k4_ss + k5_ss)) + (k3_ss * k5_ss * (k1_ss + k8_ss)))
    values['x4_ss'] = x4_ss
    E1_i = (x1_i / (x1_i + x2_i + x3_i + x4_i))
    values['E1_i'] = E1_i
    E1_ss = (x1_ss / (x1_ss + x2_ss + x3_ss + x4_ss))
    values['E1_ss'] = E1_ss
    E2_i = (x2_i / (x1_i + x2_i + x3_i + x4_i))
    values['E2_i'] = E2_i
    E2_ss = (x2_ss / (x1_ss + x2_ss + x3_ss + x4_ss))
    values['E2_ss'] = E2_ss
    E3_i = (x3_i / (x1_i + x2_i + x3_i + x4_i))
    values['E3_i'] = E3_i
    E3_ss = (x3_ss / (x1_ss + x2_ss + x3_ss + x4_ss))
    values['E3_ss'] = E3_ss
    E4_i = (x4_i / (x1_i + x2_i + x3_i + x4_i))
    values['E4_i'] = E4_i
    E4_ss = (x4_ss / (x1_ss + x2_ss + x3_ss + x4_ss))
    values['E4_ss'] = E4_ss
    ICaK = (ICaK_junc + ICaK_sl)
    values['ICaK'] = ICaK
    ICaL_junc_sigmoided = (1 - (1 / (1 + ((ICaL_junc_positive / 0.45) ** 4.5))))
    values['ICaL_junc_sigmoided'] = ICaL_junc_sigmoided
    ICaNa = (ICaNa_junc + ICaNa_sl)
    values['ICaNa'] = ICaNa
    ICa_tot = (ICaL + ICaK + ICaNa)
    values['ICa_tot'] = ICa_tot
    IK_tot = (((Ito + IKr + IKs + IK1) - (2 * INaK)) + ICaK + IKb + Istim)
    values['IK_tot'] = IK_tot
    JncxCa_i = ((E2_i * k2_i) - (E1_i * k1_i))
    values['JncxCa_i'] = JncxCa_i
    JncxCa_ss = ((E2_ss * k2_ss) - (E1_ss * k1_ss))
    values['JncxCa_ss'] = JncxCa_ss
    JncxNa_i = (((3 * ((E4_i * k7_i) - (E1_i * k8_i))) + (E3_i * k4pp_i)) - (E2_i * k3pp_i))
    values['JncxNa_i'] = JncxNa_i
    JncxNa_ss = (((3 * ((E4_ss * k7_ss) - (E1_ss * k8_ss))) + (E3_ss * k4pp_ss)) - (E2_ss * k3pp_ss))
    values['JncxNa_ss'] = JncxNa_ss
    Jrel_inact_inf = (1 / (1 + (ICaL_junc_sigmoided / (1e-3))))
    values['Jrel_inact_inf'] = Jrel_inact_inf
    Jrel_inact_inf2 = (1 / (1 + (ICaL_junc_sigmoided / (6e-4))))
    values['Jrel_inact_inf2'] = Jrel_inact_inf2
    Jrel_inf = ((a_rel * ICaL_junc_sigmoided) / (1.0 + ((directRelMidpoint / casr) ** 7.72672)))
    values['Jrel_inf'] = Jrel_inf
    INaCa_i = ((1 - INaCa_fractionSS) * Gncx * allo_i * ((zna * JncxNa_i) + (zca * JncxCa_i)))
    values['INaCa_i'] = INaCa_i
    INaCa_ss = (INaCa_fractionSS * Gncx * allo_ss * ((zna * JncxNa_ss) + (zca * JncxCa_ss)))
    values['INaCa_ss'] = INaCa_ss
    ICa_tot_junc = ((ICaL_junc + ICab_junc + IpCa_junc) - (2 * INaCa_ss))
    values['ICa_tot_junc'] = ICa_tot_junc
    ICa_tot_sl = ((ICaL_sl + ICab_sl + IpCa_sl) - (2 * INaCa_i))
    values['ICa_tot_sl'] = ICa_tot_sl
    INaCa = (INaCa_i + INaCa_ss)
    values['INaCa'] = INaCa
    INa_tot_junc = (INaj + INaLj + INabj + (3 * INaCa_ss) + (3 * INaKj) + ICaNa_junc)
    values['INa_tot_junc'] = INa_tot_junc
    INa_tot_sl = (INasl + INaLsl + INabsl + (3 * INaCa_i) + (3 * INaKsl) + ICaNa_sl)
    values['INa_tot_sl'] = INa_tot_sl
    I_tot = (INa_tot_junc + INa_tot_sl + ICa_tot_junc + ICa_tot_sl + IK_tot + ICl_tot)
    values['I_tot'] = I_tot
    derivatives = np.empty_like(state)
    derivatives[0] = (-I_tot)
    derivatives[1] = ((CaMK__alpha * bound * CaMK_active) - (CaMK__beta * camk_trap * (0.1 + ((0.9 * PP1_tot) / 0.1371))))
    derivatives[2] = ((CaMK_Phos_ss_ICaL - camk_f_ICaL) / tau_cal)
    derivatives[3] = ((CaMK_Phos_ss_RyR - camk_f_RyR) / tau_ryr)
    derivatives[4] = ((CaMK_Phos_ss_PLB - camk_f_PLB) / tau_plb)
    derivatives[5] = ((alpha_serca * bound_serca * casig_SERCA_act) - (CaMK__beta * casig_serca_trap * (0.1 + ((0.9 * PP1_tot) / 0.1371))))
    derivatives[6] = d_buffers_NaBj
    derivatives[7] = d_buffers_NaBsl
    derivatives[8] = d_buffers_TnClow
    derivatives[9] = d_buffers_TnCHc
    derivatives[10] = ((kon_tnchmg * mgi * ((Bmax_TnChigh - buffers_TnCHc) - buffers_TnCHm)) - (koff_tnchmg * buffers_TnCHm))
    derivatives[11] = d_buffers_CaM
    derivatives[12] = d_buffers_Myosin_ca
    derivatives[13] = ((kon_myomg * mgi * ((Bmax_myosin - buffers_Myosin_ca) - buffers_Myosin_mg)) - (koff_myomg * buffers_Myosin_mg))
    derivatives[14] = d_buffers_SRB
    derivatives[15] = d_buffers_SLLj
    derivatives[16] = d_buffers_SLLsl
    derivatives[17] = d_buffers_SLHj
    derivatives[18] = d_buffers_SLHsl
    derivatives[19] = d_buffers_Csqn
    derivatives[20] = ((Jserca - (((J_SRleak * vmyo) / vsr) + J_SRCarel)) - d_buffers_Csqn)
    derivatives[21] = (((((-INa_tot_junc) * Cmem) / (vjunc * F)) + ((J_na_juncsl / vjunc) * (nasl - naj))) - d_buffers_NaBj)
    derivatives[22] = (((((-INa_tot_sl) * Cmem) / (vsl * F)) + ((J_na_juncsl / vsl) * (naj - nasl)) + ((J_na_slmyo / vsl) * (nai - nasl))) - d_buffers_NaBsl)
    derivatives[23] = ((J_na_slmyo / vmyo) * (nasl - nai))
    derivatives[24] = (((-IK_tot) * Cmem) / (vmyo * F))
    derivatives[25] = ((ICl_tot * Cmem) / (vmyo * F))
    derivatives[26] = ((((((-ICa_tot_junc) * Cmem) / (vjunc * 2 * F)) + ((J_ca_juncsl / vjunc) * (casl - caj))) - J_CaB_junction) + ((J_SRCarel * vsr) / vjunc) + ((J_SRleak * vmyo) / vjunc))
    derivatives[27] = (((((-ICa_tot_sl) * Cmem) / (vsl * 2 * F)) + ((J_ca_juncsl / vsl) * (caj - casl)) + ((J_ca_slmyo / vsl) * (cai - casl))) - J_CaB_sl)
    derivatives[28] = (((((-Jserca) * vsr) / vmyo) - J_CaB_cytosol) + ((J_ca_slmyo / vmyo) * (casl - cai)))
    derivatives[29] = ((mss - m) / tm)
    derivatives[30] = ((hss - h) / th)
    derivatives[31] = ((jss - j) / tj)
    derivatives[32] = ((hssp - hp) / th)
    derivatives[33] = ((jss - jp) / tjp)
    derivatives[34] = ((mss_P - m_P) / tm)
    derivatives[35] = ((hss_P - h_P) / th)
    derivatives[36] = ((jss_P - j_P) / tj)
    derivatives[37] = ((hssp_P - hp_P) / th)
    derivatives[38] = ((jssp_P - jp_P) / tjp)
    derivatives[39] = ((mLss - mL) / tmL)
    derivatives[40] = ((hLss - hL) / thL)
    derivatives[41] = ((hLssp - hLp) / thLp)
    derivatives[42] = ((dss - d) / td)
    derivatives[43] = ((fss - ff) / tff)
    derivatives[44] = ((fss - fs) / tfs)
    derivatives[45] = ((fcass - fcaf) / tfcaf)
    derivatives[46] = ((fcass - fcas) / tfcas)
    derivatives[47] = ((jcass - jca) / tjca)
    derivatives[48] = ((fss - ffp) / tffp)
    derivatives[49] = ((fcass - fcafp) / tfcafp)
    derivatives[50] = ((anca * k2n) - (nca * km2n))
    derivatives[51] = ((anca_i * k2n) - (nca_i * km2n))
    derivatives[52] = ((((-ical_pureCDI_junc) * sigmoidTransition) / tauTransition) + ((1 - ical_pureCDI_junc) * rateRecovery))
    derivatives[53] = ((((-ical_pureCDI_sl) * sigmoidTransition2) / tauTransition2) + ((1 - ical_pureCDI_sl) * rateRecovery))
    derivatives[54] = ((dPss - d_P) / td)
    derivatives[55] = ((fss_P - ff_P) / tff)
    derivatives[56] = ((fss_P - fs_P) / tfs)
    derivatives[57] = ((fcass_P - fcaf_P) / tfcaf)
    derivatives[58] = ((fcass_P - fcas_P) / tfcas)
    derivatives[59] = ((fBPss - fBPf) / tffp)
    derivatives[60] = ((fcaBPss - fcaBPf) / tfcafp)
    derivatives[61] = ((xtoss - xtos) / tauxtos)
    derivatives[62] = ((ytoss - ytos) / tauytos)
    derivatives[63] = ((xtoss - xtof) / tauxtof)
    derivatives[64] = ((ytoss - ytof) / tauytof)
    derivatives[65] = ((xtoss_p - xtos_p) / tauxtos)
    derivatives[66] = ((xtoss_p - xtof_p) / tauxtof)
    derivatives[67] = ((ytoss - ytos_p) / tauytos_p)
    derivatives[68] = ((ytoss - ytof_p) / tauytof_p)
    derivatives[69] = ((IKr__beta * C1) - (IKr__alpha * C0))
    derivatives[70] = (((IKr__alpha * C0) + (IKr__beta_1 * C2)) - ((IKr__beta + alpha_1) * C1))
    derivatives[71] = (((alpha_1 * C1) + (beta_2 * O) + (beta_ItoC2 * I)) - ((IKr__beta_1 + alpha_2 + alpha_C2ToI) * C2))
    derivatives[72] = (((alpha_2 * C2) + (beta_i * I)) - ((beta_2 + alpha_i) * O))
    derivatives[73] = (((alpha_C2ToI * C2) + (alpha_i * O)) - ((beta_ItoC2 + beta_i) * I))
    derivatives[74] = ((xsss_junc - xs_junc) / tauxs_junc)
    derivatives[75] = ((xsss_sl - xs_sl) / tauxs_sl)
    derivatives[76] = ((Jrel_inf - jrel_icaldep_act) / tau_rel)
    derivatives[77] = ((Jrel_inact_inf - jrel_icaldep_f1) / tauInact)
    derivatives[78] = ((Jrel_inact_inf2 - jrel_icaldep_f2) / tauInact2)
    derivatives[79] = (((kim * RIcleft) - (kiSRCa * caTransFactor * ((caj / 1) ** caExpFactor) * ryr_R)) - ((caTransFactor2 * koSRCa * ((caj / 1) ** caExpFactor2) * ryr_R) - (kom * ryr_O)))
    derivatives[80] = (((caTransFactor2 * koSRCa * ((caj / 1) ** caExpFactor2) * ryr_R) - (kom * ryr_O)) - ((kiSRCa * caTransFactor * ((caj / 1) ** caExpFactor) * ryr_O) - (kim * ryr_I)))
    derivatives[81] = (((kiSRCa * caTransFactor * ((caj / 1) ** caExpFactor) * ryr_O) - (kim * ryr_I)) - ((kom * ryr_I) - (caTransFactor2 * koSRCa * ((caj / 1) ** caExpFactor2) * RIcleft)))
    derivatives[82] = ((RI_to_CI * RIcleft) - (CI_to_RI * ryr_CaRI))
    derivatives[83] = (((kim * RIcleft_p) - (kiSRCa * caTransFactor * ((caj / 1) ** caExpFactor) * ryr_R_p)) - ((caTransFactor2p * koSRCa * ((caj / 1) ** caExpFactor2) * ryr_R_p) - (kom * ryr_O_p)))
    derivatives[84] = (((caTransFactor2p * koSRCa * ((caj / 1) ** caExpFactor2) * ryr_R_p) - (kom * ryr_O_p)) - ((kiSRCa * caTransFactor * ((caj / 1) ** caExpFactor) * ryr_O_p) - (kim * ryr_I_p)))
    derivatives[85] = (((kiSRCa * caTransFactor * ((caj / 1) ** caExpFactor) * ryr_O_p) - (kim * ryr_I_p)) - ((kom * ryr_I_p) - (caTransFactor2p * koSRCa * ((caj / 1) ** caExpFactor2) * RIcleft_p)))
    derivatives[86] = ((RI_to_CI * RIcleft_p) - (CI_to_RI * ryr_CaRI_p))
    derivatives[87] = ((xb_ws - xb_su) - xb_su_gamma)
    derivatives[88] = (((xb_uw - xb_wu) - xb_ws) - xb_wu_gamma)
    derivatives[89] = d_contraction_Ca_TRPN
    derivatives[90] = ((ktm_block * min(100, (contraction_Ca_TRPN ** ((-nperm) / 2))) * XU) - (ktm_unblock * (contraction_Ca_TRPN ** (nperm / 2)) * contraction_TmBlocked))
    derivatives[91] = ((A * lambda_rate) - (cds * contraction_ZETAS))
    derivatives[92] = ((A * lambda_rate) - (cdw * contraction_ZETAW))
    values['v'] = v
    values['camk_trap'] = camk_trap
    values['camk_f_ICaL'] = camk_f_ICaL
    values['camk_f_RyR'] = camk_f_RyR
    values['camk_f_PLB'] = camk_f_PLB
    values['casig_serca_trap'] = casig_serca_trap
    values['buffers_NaBj'] = buffers_NaBj
    values['buffers_NaBsl'] = buffers_NaBsl
    values['buffers_TnClow'] = buffers_TnClow
    values['buffers_TnCHc'] = buffers_TnCHc
    values['buffers_TnCHm'] = buffers_TnCHm
    values['buffers_CaM'] = buffers_CaM
    values['buffers_Myosin_ca'] = buffers_Myosin_ca
    values['buffers_Myosin_mg'] = buffers_Myosin_mg
    values['buffers_SRB'] = buffers_SRB
    values['buffers_SLLj'] = buffers_SLLj
    values['buffers_SLLsl'] = buffers_SLLsl
    values['buffers_SLHj'] = buffers_SLHj
    values['buffers_SLHsl'] = buffers_SLHsl
    values['buffers_Csqn'] = buffers_Csqn
    values['casr'] = casr
    values['naj'] = naj
    values['nasl'] = nasl
    values['nai'] = nai
    values['ki'] = ki
    values['cli'] = cli
    values['caj'] = caj
    values['casl'] = casl
    values['cai'] = cai
    values['m'] = m
    values['h'] = h
    values['j'] = j
    values['hp'] = hp
    values['jp'] = jp
    values['m_P'] = m_P
    values['h_P'] = h_P
    values['j_P'] = j_P
    values['hp_P'] = hp_P
    values['jp_P'] = jp_P
    values['mL'] = mL
    values['hL'] = hL
    values['hLp'] = hLp
    values['d'] = d
    values['ff'] = ff
    values['fs'] = fs
    values['fcaf'] = fcaf
    values['fcas'] = fcas
    values['jca'] = jca
    values['ffp'] = ffp
    values['fcafp'] = fcafp
    values['nca'] = nca
    values['nca_i'] = nca_i
    values['ical_pureCDI_junc'] = ical_pureCDI_junc
    values['ical_pureCDI_sl'] = ical_pureCDI_sl
    values['d_P'] = d_P
    values['ff_P'] = ff_P
    values['fs_P'] = fs_P
    values['fcaf_P'] = fcaf_P
    values['fcas_P'] = fcas_P
    values['fBPf'] = fBPf
    values['fcaBPf'] = fcaBPf
    values['xtos'] = xtos
    values['ytos'] = ytos
    values['xtof'] = xtof
    values['ytof'] = ytof
    values['xtos_p'] = xtos_p
    values['xtof_p'] = xtof_p
    values['ytos_p'] = ytos_p
    values['ytof_p'] = ytof_p
    values['C0'] = C0
    values['C1'] = C1
    values['C2'] = C2
    values['O'] = O
    values['I'] = I
    values['xs_junc'] = xs_junc
    values['xs_sl'] = xs_sl
    values['jrel_icaldep_act'] = jrel_icaldep_act
    values['jrel_icaldep_f1'] = jrel_icaldep_f1
    values['jrel_icaldep_f2'] = jrel_icaldep_f2
    values['ryr_R'] = ryr_R
    values['ryr_O'] = ryr_O
    values['ryr_I'] = ryr_I
    values['ryr_CaRI'] = ryr_CaRI
    values['ryr_R_p'] = ryr_R_p
    values['ryr_O_p'] = ryr_O_p
    values['ryr_I_p'] = ryr_I_p
    values['ryr_CaRI_p'] = ryr_CaRI_p
    values['contraction_XS'] = contraction_XS
    values['contraction_XW'] = contraction_XW
    values['contraction_Ca_TRPN'] = contraction_Ca_TRPN
    values['contraction_TmBlocked'] = contraction_TmBlocked
    values['contraction_ZETAS'] = contraction_ZETAS
    values['contraction_ZETAW'] = contraction_ZETAW
    return  derivatives, I_tot

def rhs(time, state, parameters=None):
    """solve_ivp-compatible derivative function."""
    return evaluate(time, state, parameters)[0]




def ventricle_bidomain_1D(
    total_time,
    cell_number,
    states,
    dt,
    sigma_i,
    sigma_e,
    beta,
    Cm,
    amp=-53.0,
    duration=1.0,
    L=100.0,
    dx=None,
    parameters=None,
    Cm_TW=0.02146,
    max_step=0.01,
    rtol=1e-6,
    atol=1e-8
):
    """
    Integrate a 1D ventricular bidomain tissue using solve_ivp.

    dt is the output sampling interval, NOT a fixed integration step.
    The solver chooses its own adaptive time steps.

    Returns
    -------
    times : ndarray
        Output times (ms).

    voltage_history : ndarray
        Membrane voltage, shape (number of output times, cell_number).
    """

    states = np.asarray(states, dtype=float).copy()

    n_states = len(STATE_NAMES)

    # ---------------------------------------------------------
    # Validate inputs
    # ---------------------------------------------------------

    if cell_number < 2:
        raise ValueError("cell_number must be at least 2.")

    if states.shape != (cell_number, n_states):
        raise ValueError(
            f"states must have shape ({cell_number}, {n_states})."
        )

    if total_time <= 0 or dt <= 0:
        raise ValueError("total_time and dt must be positive.")

    if beta <= 0 or Cm <= 0 or Cm_TW <= 0:
        raise ValueError("beta, Cm and Cm_TW must be positive.")

    if max_step <= 0:
        raise ValueError("max_step must be positive.")

    if not (0 <= duration < 1000):
        raise ValueError("duration must be between 0 and 1000 ms.")

    if dx is None:
        dx = L * 1e-6

    if dx <= 0:
        raise ValueError("dx must be positive.")

    if parameters is None:
        parameters = {}

    # ---------------------------------------------------------
    # Construct bidomain matrices ONCE
    # ---------------------------------------------------------

    length = cell_number * dx

    D2, x, dx_actual = bd.spatial_derivative_1D(
        cell_number,
        length
    )

    Ki, Ke = bd.conductivity_matrices(
        D2,
        sigma_i,
        sigma_e
    )

    # ---------------------------------------------------------
    # Define the solve_ivp right-hand side
    # ---------------------------------------------------------

    def tissue_rhs(t, y, stimulus_on):

        # Reshape flattened solver vector into individual cells
        cell_states = y.reshape(cell_number, n_states)

        # Allocate derivatives and ionic currents
        dstate = np.zeros_like(cell_states)
        I_ion = np.zeros(cell_number)

        # Evaluate T-WORLD in every cell
        for i in range(cell_number):

            derivatives, total_current = evaluate(
                t,
                cell_states[i],
                parameters
            )

            dstate[i] = derivatives

            # Convert T-WORLD current into physical current density.
            # I_tot is outward-positive, as required by bidomain.
            I_ion[i] = Cm_TW * total_current

        # -----------------------------------------------------
        # External stimulus
        # -----------------------------------------------------

        I_stim = np.zeros(cell_number)

        if stimulus_on:
            # T-WORLD: negative amp is inward/depolarising.
            # Bidomain: positive I_stim is inward/depolarising.
            I_stim[0] = -Cm_TW * amp

        # -----------------------------------------------------
        # Bidomain voltage derivative
        # -----------------------------------------------------

        Vm = cell_states[:, 0]

        phi_e = bd.solve_extracellular_potential(
            Ki,
            Ke,
            Vm
        )

        dVm_dt = bd.calculate_dVm(
            Ki,
            Vm,
            phi_e,
            beta,
            Cm,
            I_ion,
            I_stim
        )

        # Replace T-WORLD's single-cell voltage derivatives
        # with the voltage derivatives from the bidomain model.
        dstate[:, 0] = dVm_dt

        # Return all state derivatives as one flattened vector
        return dstate.ravel()

    # ---------------------------------------------------------
    # Output sampling times
    # ---------------------------------------------------------

    n_output = int(np.ceil(total_time / dt))

    times = np.linspace(
        0.0,
        total_time,
        n_output + 1
    )

    voltage_history = np.zeros(
        (len(times), cell_number)
    )

    voltage_history[0] = states[:, 0]

    # ---------------------------------------------------------
    # Split integration at stimulus switching times
    # ---------------------------------------------------------

    # This ensures that the adaptive solver does not step across
    # the abrupt start or end of a stimulus pulse.

    boundaries = [0.0, total_time]

    if duration > 0:

        n_periods = int(np.ceil(total_time / 1000.0))

        for k in range(n_periods + 1):

            pulse_start = k * 1000.0
            pulse_end = pulse_start + duration

            if 0 < pulse_start < total_time:
                boundaries.append(pulse_start)

            if 0 < pulse_end < total_time:
                boundaries.append(pulse_end)

    boundaries = sorted(set(boundaries))

    # Initial condition for the complete tissue
    y0 = states.ravel()

    # ---------------------------------------------------------
    # Integrate each interval
    # ---------------------------------------------------------

    for a, b in zip(boundaries[:-1], boundaries[1:]):

        # Determine whether the stimulus is on during this interval
        midpoint = 0.5 * (a + b)

        stimulus_on = (
            duration > 0 and
            (midpoint % 1000.0) < duration
        )

        solution = solve_ivp(
            tissue_rhs,
            (a, b),
            y0,
            args=(stimulus_on,),
            method="BDF",
            dense_output=True,
            max_step=max_step,
            rtol=rtol,
            atol=atol
        )

        if not solution.success:
            raise RuntimeError(
                f"solve_ivp failed between {a} and {b} ms: "
                f"{solution.message}"
            )

        # Select output times belonging to this interval
        tolerance = 1e-12

        indices = np.where(
            (times >= a - tolerance) &
            (times <= b + tolerance)
        )[0]

        # Extract only the membrane voltages.
        # We do not need to store all 93 states at every output time.
        for j in indices:

            y_at_time = solution.sol(times[j])

            state_at_time = y_at_time.reshape(
                cell_number,
                n_states
            )

            voltage_history[j] = state_at_time[:, 0]

        # Pass the final state of this interval into the next one
        y0 = solution.y[:, -1]

    return times, voltage_history
