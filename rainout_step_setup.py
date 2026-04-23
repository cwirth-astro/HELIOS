import numpy as np
import pandas as pd
import os
from source import tools as tls

def print_config_file(path, profile, calculation, chem_output, cond_output, monitor, abundance):
    with open(path, 'w') as f:
        f.write(f'#Atmospheric profile input file\n{profile}\n\n')
        f.write(f"#Chemistry calculation type (gas phase only = g, equilibrium condensation = ce, rainout condensation = cr)\n{calculation}\n\n")
        f.write(f"#Chemistry output file\n{chem_output} {cond_output}\n\n")
        f.write(f"#Monitor output file\n{monitor}\n\n")
        f.write(f"FastChem console verbose level\n1\n\n")
        f.write(f"#Output mixing ratios (MR) or particle number densities (ND, default)\nMR\n\n")
        f.write(f"#Element abundance file\n{abundance}\n\n")
        f.write(f"#Species data files\n/home/cwirth/FastChem/input/logK/logK.dat /home/cwirth/FastChem/input/logK/logK_condensates.dat\n\n")
        f.write(f"#Accuracy of chemistry iteration\n1.0e-4\n\n")
        f.write(f"#Accuracy of element conservation\n1.0e-3\n\n")
        f.write(f"#Max number of chemistry iterations\n80000 \n\n")
        f.write(f"#Max number internal solver iterations\n20000\n\n")

# set list of abundances to use
abundance_list = ['Bath(H4)', 'EH3', 'EL3', 'Krymka(LL3)', 'R3']
escaped_abundance_list = ['Bath\(H4\)', 'EH3', 'EL3', 'Krymka\(LL3\)', 'R3']
# set list of equilibrium temperatures to use
t_eq_list = np.arange(200, 1800, 200)

print_tp = False
print_config= True
run_fastchem = True

for abundance, escaped_abundance in zip(abundance_list, escaped_abundance_list):
    for t_eq in t_eq_list:
        # set up paths
        helios_tp_path = f"output/grid_day/{t_eq}K_{abundance}_day/{t_eq}K_{abundance}_day_tp.dat"
        config_path = f"output/grid_day/{t_eq}K_{abundance}_day/chemistry/config.dat"
        chem_output_path = f"output/grid_day/{t_eq}K_{abundance}_day/chemistry/chem.dat"
        tp_output_path = f"output/grid_day/{t_eq}K_{abundance}_day/chemistry/tp.dat"
        cond_output_path = f"output/grid_day/{t_eq}K_{abundance}_day/chemistry/condensates.dat"
        monitor_path = f"output/grid_day/{t_eq}K_{abundance}_day/chemistry/monitor.dat"
        if abundance == 'asplund_m3':
            abundance_path = "/home/cwirth/FastChem/input/element_abundances/asplund_2020_times1000.dat"
        else:
            abundance_path = f"/home/cwirth/FastChem/input/element_abundances/bulk_comp/{abundance}.dat"

    
        # create chemistry directory if it doesn't exist
        if not os.path.exists(f"output/grid_day/{t_eq}K_{abundance}_day/chemistry"):
            os.makedirs(f"output/grid_day/{t_eq}K_{abundance}_day/chemistry")
        
        # read helios tp
        press, temp, cpress0, ctemp0, cpress1, ctemp1, cpress2, ctemp2, cpress3, ctemp3 = tls.read_helios_tp(helios_tp_path)

        # print tp profile to chemistry dir
        if print_tp:
            arr = np.column_stack((np.array(press), np.array(temp)))
            np.savetxt(tp_output_path, arr, header='pressure(bar) temperature(K)', fmt='%.6e %.6f', comments='')
            print(f"Wrote {arr.shape[0]} rows to {tp_output_path}")

        # print config file
        if print_config:
            print_config_file(config_path, tp_output_path, 'cr', chem_output_path, cond_output_path, monitor_path, abundance_path)
            print(f"Wrote config file to {config_path}")

        if run_fastchem:
            # run fastchem
            os.system(f"/home/cwirth/FastChem/fastchem {config_path}")