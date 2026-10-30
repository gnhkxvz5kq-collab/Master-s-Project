import Trovato_Purkinje_Model_corrected as tpf
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
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

'''Run the Trovato Purkinje functions to check for any errors and see the outputs'''
time = 0 #ms
# Collect all currents at every saved time point

# Obtain the 46 initial conditions
y0 = tpf.initial_state()

# Simulate for 1000 ms
solution = solve_ivp(
    fun=lambda t, y: tpf.Trovato_Purkinje_Model(t, y),
    t_span=(0.0, 1000.0),
    y0=y0,
    method="BDF",
    max_step=0.1,
    rtol=1e-6,
    atol=1e-9,
)

if not solution.success:
    raise RuntimeError(solution.message)

# Membrane potential over time
time = solution.t
voltage = solution.y[tpf.STATE_NAMES.index("v"), :]
# Example: maximum membrane voltage
print("Maximum voltage:", np.max(voltage), "mV")
plt.plot(time, voltage)
plt.xlabel("time / ms")
plt.ylabel("Voltage / mV")
plt.plot()
plt.savefig('HPC_V_vs_T_test_fig.png')
plt.close()

current_data = {}

for i, t in enumerate(solution.t):
    state = solution.y[:, i]

    _, diagnostics = tpf.Trovato_Purkinje_Model(
        t, state, return_diagnostics=True
    )

    for name, value in diagnostics["currents"].items():
        current_data.setdefault(name, []).append(float(value))

# Store results in a DataFrame
currents_df = pd.DataFrame(current_data, index=solution.t)
currents_df.index.name = "time_ms"

# Save for later analysis
currents_df.to_csv("purkinje_currents.csv")

print(currents_df.head())
print("Currents collected:", list(currents_df.columns))




plot_dir = Path("current_plots")
plot_dir.mkdir(exist_ok=True)

for current_name in currents_df.columns:
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.plot(
        currents_df.index,
        currents_df[current_name],
        linewidth=1.2
    )

    ax.set_xlabel("Time /ms")
    ax.set_ylabel(f"{current_name} / mA ")
    ax.set_title(f"{current_name} vs time")

    fig.tight_layout()
    fig.savefig(plot_dir / f"HPC_{current_name}_vs_time.png", dpi=200)
    plt.close(fig)
