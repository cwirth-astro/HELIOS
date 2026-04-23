import pandas as pd
import numpy as np

# temperatures and compositions
temps = np.arange(200, 1800, 200)
comps = ["Bath(H4)", "EH3", "EL3", "Krymka(LL3)", "R3", "asplund_m3"]

#dayside results folder
dir_dayside = "output/grid_day"
#heating profile results folder
dir_heating = "input/heating/grid_night"

#gaussian cutoff parameters (in bar)
g_mean = 4e-1
g_sigma = 1.6e-1

#heat redistribution parameter (TODO: iterate this)
f = 0.35

for T in temps:
    for comp in comps:
        dir_output = f'{dir_dayside}/{T}K_{comp}_day'
        # load in T-P profile
        tp = pd.read_csv(f'{dir_output}/{T}K_{comp}_day_tp.dat', sep='\\s+', skiprows=1)
        #load in integrated flux per layer
        int_flux = pd.read_csv(f'{dir_output}/{T}K_{comp}_day_integrated_flux.dat', sep='\\s+', skiprows=2)
        #load in column mass
        colmass = pd.read_csv(f'{dir_output}/{T}K_{comp}_day_colmass_mu_cp_kappa_entropy.dat', sep='\\s+', skiprows=1)

        #calculate the total flux to the night side
        flux_to_night = int_flux['F_down'][71] * (0.5-f)/f

        #pressure column
        press_col_bar = (tp['press.[10^-6bar]'][1:]/10**6).reset_index(drop=True)  # in bar
        #gaussian factor for wind efficiency
        prop_to_gauss = np.exp(-0.5 * ((press_col_bar - g_mean) / g_sigma)**2) / np.sum(np.exp(-0.5 * ((press_col_bar - g_mean) / g_sigma)**2))

        #mass per layer column
        prop_to_mass = colmass['delta_col.mass[g'] / np.sum(colmass['delta_col.mass[g'])

        #multiplication of the two
        prop_to_mass_gauss = prop_to_mass * prop_to_gauss / np.sum(prop_to_mass * prop_to_gauss)

        #turn this into the heating rate per layer
        heating_rate = flux_to_night * prop_to_mass_gauss


        heating_file_df = pd.DataFrame({
            'pressure': press_col_bar,
            'heating_rate': heating_rate
        })
        #print heating file df to file
        heating_file_df.to_csv(f'{dir_heating}/{T}K_{comp}_night_heating_rate.dat', sep='\t', index=False)