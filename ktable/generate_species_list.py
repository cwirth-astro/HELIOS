import os

BASE_DIR = "/project/diana8/cwirth/HELIOS/opacity"   # change if needed
OUTPUT_FILE = "input/individual_species_1.dat"

species_dirs = []

# Loop through directories in BASE_DIR
for name in sorted(os.listdir(BASE_DIR)):
    full_path = os.path.join(BASE_DIR, name)
    if os.path.isdir(full_path):
        species_dirs.append((name, full_path))

# Write output
with open(OUTPUT_FILE, "w") as f:
    f.write("name_of_species_to_include     path_to_files\n")
    for name, path in species_dirs:
        f.write(f"{name:<20} {path}\n")

print(f"Wrote {len(species_dirs)} species to {OUTPUT_FILE}")