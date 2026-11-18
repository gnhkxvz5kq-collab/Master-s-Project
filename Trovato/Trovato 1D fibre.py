import numpy as np
import matplotlib.pyplot as plt

from Purkinje_functions_tissues import (
    initial_state,
    STATE_NAMES,
    purkinje_bidomain
)


def run_and_plot_purkinje(
    cell_number,
    total_time,
    dt,
    sigma_i,
    sigma_e,
    beta,
    Cm,
    amp=-40.0,
    duration=1.0,
    L=0.0164,
    radius=0.00175,
    cell_length=None,
    threshold=-30.0
):

    # --------------------------------------------------
    # 1. Initialise all Purkinje cells
    # --------------------------------------------------

    single_cell_state = initial_state()

    states = np.tile(
        single_cell_state,
        (cell_number, 1)
    )

    # --------------------------------------------------
    # 2. Run bidomain simulation
    # --------------------------------------------------

    times, voltage_history, I_ion, I_Na, I_Ca= (
        purkinje_bidomain(
            total_time,
            cell_number,
            states,
            dt=dt,
            sigma_i=sigma_i,
            sigma_e=sigma_e,
            beta=beta,
            Cm=Cm,
            cell_length=cell_length,
            amp=amp,
            duration=duration,
            L=L,
            radius=radius,
            #return_diagnostics=True
        )
    )

    # Spatial positions
    if cell_length is None:
        cell_length = L * 1e-2

    positions = (
        np.arange(cell_number) + 0.5
    ) * cell_length

    positions_mm = positions * 1e3

    # Select cells distributed along the fibre
    selected_cells = np.arange(0, cell_number, 5)

    # --------------------------------------------------
    # 3. Plot voltage traces
    # --------------------------------------------------

    plt.figure(figsize=(10, 6))

    for i in selected_cells:
        plt.plot(
            times,
            voltage_history[:, i],
            label=f"Cell {i}"
        )

    plt.xlabel("Time (ms)")
    plt.ylabel("Membrane potential (mV)")
    plt.title("Purkinje fibre: membrane potential at selected cells")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # --------------------------------------------------
    # 4. Plot space-time voltage propagation
    # --------------------------------------------------

    plt.figure(figsize=(10, 6))

    plt.pcolormesh(
        positions_mm,
        times,
        voltage_history,
        shading="auto"
    )

    plt.colorbar(label="Membrane potential (mV)")
    plt.xlabel("Position along fibre (mm)")
    plt.ylabel("Time (ms)")
    plt.title("Purkinje fibre: action potential propagation")
    plt.tight_layout()
    plt.show()

    # --------------------------------------------------
    # 6. Plot ionic currents
    # --------------------------------------------------

    # I_ion = diagnostics["I_ion"]
    # I_Na = diagnostics["I_Na"]
    # I_Ca = diagnostics["I_Ca"]
    plt.figure(figsize=(10, 5))

    for i in selected_cells:
        plt.plot(
            times[:-1],
            I_ion[:, i],
            label=f"Cell {i}"
        )

    plt.xlabel("Time (ms)")
    plt.ylabel("I_ion (micro A)")
    plt.title("Trovato I_ion at selected cells")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(10, 5))

    for i in selected_cells:
        plt.plot(
            times[:-1],
            I_Na[:, i],
            label=f"Cell {i}"
        )

    plt.xlabel("Time (ms)")
    plt.ylabel("I_Na (microA)")
    plt.title("Trovato I_Na at selected cells")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 5))

    for i in selected_cells:
        plt.plot(
            times[:-1],
            I_Ca[:, i],
            label=f"Cell {i}"
        )

    plt.xlabel("Time (ms)")
    plt.ylabel("I_Ca (micro A)")
    plt.title("Trovato I_Ca at selected cells")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()
    

    # --------------------------------------------------
    # 7. Calculate activation times
    # --------------------------------------------------

    activation_times = np.full(cell_number, np.nan)

    for i in range(cell_number):

        V = voltage_history[:, i]

        crossings = np.where(
            (V[:-1] < threshold) &
            (V[1:] >= threshold)
        )[0]

        if len(crossings) > 0:

            j = crossings[0]

            # Linear interpolation between time points
            t1 = times[j]
            t2 = times[j + 1]

            V1 = V[j]
            V2 = V[j + 1]

            activation_times[i] = (
                t1
                + (threshold - V1)
                * (t2 - t1)
                / (V2 - V1)
            )

    # --------------------------------------------------
    # 8. Plot activation time against position
    # --------------------------------------------------

    activated = np.isfinite(activation_times)

    plt.figure(figsize=(9, 5))

    plt.scatter(
        positions_mm[activated],
        activation_times[activated]
    )

    plt.xlabel("Position along fibre (mm)")
    plt.ylabel("Activation time (ms)")
    plt.title(f"Activation times (threshold = {threshold} mV)")
    plt.grid()
    plt.tight_layout()
    plt.show()

    # --------------------------------------------------
    # 9. Estimate conduction velocity
    # --------------------------------------------------

    # Exclude the stimulated cell from the fit
    fit_mask = activated.copy()
    fit_mask[0] = False

    if np.sum(fit_mask) >= 2:

        slope, intercept = np.polyfit(
            positions_mm[fit_mask],
            activation_times[fit_mask],
            1
        )

        if slope > 0:

            # 1 mm/ms = 1 m/s
            conduction_velocity = 1.0 / slope

            print(
                f"Estimated conduction velocity: "
                f"{conduction_velocity:.3f} m/s"
            )

            plt.figure(figsize=(9, 5))

            plt.scatter(
                positions_mm[fit_mask],
                activation_times[fit_mask],
                label="Measured activation times"
            )

            plt.plot(
                positions_mm[fit_mask],
                slope * positions_mm[fit_mask] + intercept,
                label="Linear fit"
            )

            plt.xlabel("Position along fibre (mm)")
            plt.ylabel("Activation time (ms)")
            plt.title("Conduction velocity estimation")
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.show()

        else:
            print("Activation-time fit has a non-positive slope.")

    else:
        print(
            "Insufficient activated cells to estimate "
            "conduction velocity."
        )

    return times, voltage_history, activation_times

results = run_and_plot_purkinje(
    cell_number= 200,
    total_time= 1000,
    dt=0.005,
    sigma_i=2 * 0.6,
    sigma_e=6 * 0.6,
    beta=115500,
    Cm=0.01,
    amp=-40.0,
    duration=2.2
)
### these values actually give graphs that work now and give good conduction values 
