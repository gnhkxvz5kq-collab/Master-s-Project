#testing the T World translation from Cellml to Python to validate that it works as expected 
"""Minimal example for running and plotting the generated T-World model."""
import matplotlib.pyplot as plt
import TWorld_model_functions as tw

# CellML default pacing: -53 uA/uF stimulus, 1 ms duration, 1000 ms period.
sol, traces = tw.simulate(
    t_span=(0.0, 1000.0),
    dt=0.1,
    method="BDF",
    rtol=1e-6,
    atol=1e-8,
    max_step=0.5,
)

tw.plot_standard_traces(traces)
plt.savefig('TWorld_Trial_1.png')

# All model outputs are available by name, for example:
# t = traces["time"]
# voltage_mV = traces["v"]
# cytosolic_calcium_mM = traces["cai"]
# active_tension_model_output = traces["Ta"]