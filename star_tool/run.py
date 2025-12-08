"""
Documentation for this file is at https://heliosexo.readthedocs.io/en/latest/sections/tutorial.html#creating-a-stellar-spectrum-file

But in brief, to download & interpolate PHOENIX spectra for ~any star:
1) Create an entry for your desired star, following the example of gj1214 below.
2) change the first input argument of fc.main_loop() to be the new star dict you created.
3) Set the output file, as desired.
4) At the terminal, run:  `python3 run.py`

"""
import sys
sys.path.append("..")
import functions as fc

# Units are [temp]=K, [log_g]=log(cm s^-2), m = log[n(Fe)/n(H)]_star - log[n(Fe)/n(H)]_sun, with n(A) being number density of A

gj1214 = {
    "data_format": "phoenix",
    "name": "gj1214",
    "temp": 3026,
    "log_g": 4.944,
    "m": 0.39
} # reference: Harpsoe et al. (2013)

sun = {
    "data_format": "ascii",
    "source_file": "./input/ascii/sun_gueymard_2003.txt",
    "name": "sun",
    "w_conversion_factor": 1e-7,
    "flux_conversion_factor": 1e10,
    "temp": 5772
}

star_of_interest = {
    "data_format": "muscles",
    "name": "star_of_interest",
    "source_file" : "./input/downloaded_muscles_spectrum.fits",
    "w_conversion_factor": 1e-8,
    "flux_conversion_factor": 1e8,
    "distance_from_Earth": 4.67517, # in pc
    "R_star": 0.366999654557, # in R_sun
    "temp": 3293.7
}

WASP15 = {
    "data_format": "ascii",
    "source_file": "../input/spectra/WASP15_spectrum.txt",
    "name": "WASP15",
    "w_conversion_factor": 1e-7,
    "flux_conversion_factor": 1e7,
    "temp": 6372
}

KELT7 = {
    "data_format": "ascii",
    "source_file": "../input/spectra/KELT7_spectrum.txt",
    "name": "KELT7",
    "w_conversion_factor": 1e-7,
    "flux_conversion_factor": 1e7,
    "temp": 6768
}

HATP30 = {
    "data_format": "ascii",
    "source_file": "../input/spectra/HATP30_spectrum.txt",
    "name": "HATP30",
    "w_conversion_factor": 1e-7,
    "flux_conversion_factor": 1e7,
    "temp": 6304
}

NGTS2 = {
    "data_format": "ascii",
    "source_file": "../input/spectra/NGTS2_spectrum.txt",
    "name": "NGTS2",
    "w_conversion_factor": 1e-7,
    "flux_conversion_factor": 1e7,
    "temp": 6478
}

TrES4 = {
    "data_format": "ascii",
    "source_file": "../input/spectra/TrES4_spectrum.txt",
    "name": "TrES4",
    "w_conversion_factor": 1e-7,
    "flux_conversion_factor": 1e7,
    "temp": 6295
}

WASP94A = {
    "data_format": "ascii",
    "source_file": "../input/spectra/WASP94A_spectrum.txt",
    "name": "WASP94A",
    "w_conversion_factor": 1e-7,
    "flux_conversion_factor": 1e7,
    "temp": 6170
}

WASP17 = {
    "data_format": "ascii",
    "source_file": "../input/spectra/WASP17_spectrum.txt",
    "name": "WASP17",
    "w_conversion_factor": 1e-7,
    "flux_conversion_factor": 1e7,
    "temp": 6650
}

HD149026 = {
    "data_format": "ascii",
    "source_file": "../input/spectra/HD149026_spectrum.txt",
    "name": "HD149026",
    "w_conversion_factor": 1e-7,
    "flux_conversion_factor": 1e7,
    "temp": 6147
}

# run the thing

for star in [WASP15, KELT7, HATP30, NGTS2, TrES4, WASP94A, WASP17, HD149026]:
    fc.main_loop(star,
                convert_to='r50_kdistr',
                opac_file_for_lambdagrid="../ktable/output/r50_kdistr/mixed/mixed_opac_kdistr.h5",
                output_file=f"{star['name']}.h5",
                plot_and_tweak='automatic',
                save_ascii='no',
                save_in_hdf5='yes')
