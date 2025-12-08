import os

BASE_DIR = "/project/diana8/cwirth/HELIOS/opacity"
OUTPUT_FILE = "input/mixed_species_table.dat"

# species that get "yes" for scattering
SCATTERING_YES = {"H2", "He", "H", "H2O", "CO", "CO2", "O2", "N2"}

species_dirs = []

# Collect directory names
for name in sorted(os.listdir(BASE_DIR)):
    full_path = os.path.join(BASE_DIR, name)
    if os.path.isdir(full_path):
        species_dirs.append(name)

with open(OUTPUT_FILE, "w") as f:
    f.write("species     absorbing  scattering  mixing_ratio\n\n")

    for name in species_dirs:
        absorbing = "yes"
        scattering = "yes" if name in SCATTERING_YES else "no"
        mixing_ratio = "FastChem"

        # Adjust spacing for neat columns
        f.write(f"{name:<12} {absorbing:<10} {scattering:<10} {mixing_ratio}\n")

print(f"Wrote {len(species_dirs)} species to {OUTPUT_FILE}")
