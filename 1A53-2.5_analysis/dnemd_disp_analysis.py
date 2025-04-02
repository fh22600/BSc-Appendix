import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np

repeats = range(1,11) ##range(1,11)  
start_times = [50, 60, 70, 75, 80, 90, 100, 110, 120, 125, 130, 140, 150, 160, 170, 175, 180, 190, 200,210, 220, 225, 230, 240, 250, 260, 270, 275, 280, 290, 300,310, 320, 325, 330, 340, 350, 360, 370, 375, 380, 390, 400, 410, 420, 425, 430, 440, 450, 460, 470, 475, 480, 490,500 ]

# Load the base reference structure
base_ref = mda.Universe('ref.pdb')
base_ca = base_ref.select_atoms('name CA')

# Create a new Universe object with only the C-alpha atoms
ca_universe = mda.Merge(base_ca)

# Save the positions of base C-alpha atoms to a file
np.savetxt("base_ca.txt", base_ca.positions)

time_points = [170] ## change to whatever timepoints add some to 70 80 90 

for tp in time_points:
    vectors = []

    for repeat in repeats:
        for time in start_times: ## check compared to dump
           
            eq_coords = mda.Universe(f'Repeat_{repeat}/ns_{time}/equil/md{time}_{tp}.pdb')
            noneq_coords = mda.Universe(f'Repeat_{repeat}/ns_{time}/md{time}_{tp}.pdb')
            align.alignto(noneq_coords, eq_coords, select='name CA')
            ca_eq = eq_coords.select_atoms('name CA').positions
            ca_noneq = noneq_coords.select_atoms('name CA').positions
            diff = ca_noneq - ca_eq
            vectors.append(diff)

    avg_vectors = np.mean(np.array(vectors), axis=0)

    # Calculate stats
    avg_disp = np.linalg.norm(avg_vectors, axis=1)
    sd = np.std(vectors, axis=0)
    se = sd / np.sqrt(len(vectors))
    avg_x = avg_vectors[:, 0]
    avg_y = avg_vectors[:, 1]
    avg_z = avg_vectors[:, 2]

    # SE of magnitude
    se_mag = (1 / avg_disp) * np.sqrt((avg_x**2)*(se[:, 0]**2) + (avg_y**2)*(se[:, 1]**2) + (avg_z**2)*(se[:, 2]**2))
    # relaxed to 1 x std, then relaxed to 0
    adjusted_avg_disp = np.where(np.abs(avg_disp) - 1*se_mag > 0, avg_disp, 0)
    # adjusted_avg_disp = np.where(np.abs(avg_disp) - 2*se_mag > 0, avg_disp, 0)
    adjusted_avg_vectors = np.where(np.abs(avg_vectors) - 1*se > 0, avg_vectors, 0)

    # Map the norms of the C-alpha vectors onto the C-alpha atoms in the new Universe object
    for atom, bfactor in zip(ca_universe.atoms, adjusted_avg_disp):
        atom.tempfactor = bfactor

    ## Write the C-alpha-only structure with mapped B-factors
    pdb_filename = f'pdbs/550_1SD_1A5{tp}.pdb'
    ca_universe.atoms.write(pdb_filename)
