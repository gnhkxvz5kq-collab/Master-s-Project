import numpy as np 
import matplotlib.pyplot as plt

from Ventricle_func_tissues import (initial_state, STATE_NAMES, ventricle_bidomain_1D)

def run_and_plot_ventricle(cell_number, 
                           total_time, 
                           dt, 
                           sigma_i, 
                           sigma_e, 
                           beta, 
                           Cm, 
                           amp=-40.0, 
                           duration=1.0,
                           L=500, 
                           cell_length=None, 
                            threshold=-30.0
):
    # --------------------------------------------------
    # 1. Initialise all ventricular cells
    # --------------------------------------------------

    single_cell_state = initial_state()

    states = np.tile(
        single_cell_state,
        (cell_number, 1)
    )

    # --------------------------------------------------
    # 2. Run bidomain simulation
    # --------------------------------------------------

    times, voltage_history = (
        ventricle_bidomain_1D(
            total_time,
            cell_number,
            states,
            dt=dt,
            sigma_i=sigma_i,
            sigma_e=sigma_e,
            beta=beta,
            Cm=Cm,
            dx=None,
            amp=amp,
            duration=duration,
            L=L,
        )
    )
        # Spatial positions
    if cell_length is None:
        cell_length = L * 1e-6

    positions = (
        np.arange(cell_number) + 0.5
    ) * cell_length

    positions_mm = positions * 1e3

    # Select cells distributed along the fibre
    selected_cells = np.arange(0, cell_number, 1)

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
    plt.title("Ventricle fibre: membrane potential at selected cells")
    plt.legend(loc = 'upper right', fontsize = 'small')
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
    plt.title("Ventricle fibre: action potential propagation")
    plt.tight_layout()
    plt.show()


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
    return times, voltage_history

results = run_and_plot_ventricle(
    cell_number=10,
    total_time=400,
    dt=0.005,
    sigma_i=0.2,
    sigma_e=0.6,
    beta=215122,
    Cm=0.02,
    amp=-53.5,
    duration=1.0
)
