import MDAnalysis as mda
from MDAnalysis.analysis import align

import numpy as np

runs = [1,2,3,4,5,6,7,8,9,10] ## 1 to 10 
times = [50, 60, 70, 75, 80, 90, 100, 110, 120, 125, 130, 140, 150, 160,170, 175, 180, 190, 200,210, 220, 225, 230, 240, 250, 260, 270, 275, 280, 290, 300,310, 320, 325,330, 340, 350, 360, 360,370, 375, 380, 390, 400, 410, 420, 425, 430, 440, 450, 460, 470, 475, 480, 490,500 ]


# Load the base reference structure
base_ref = mda.Universe('ref.pdb')
base_ca = base_ref.select_atoms('name CA')

# Create a new Universe object with only the C-alpha atoms
ca_universe = mda.Merge(base_ca)

# Save the positions of base C-alpha atoms to a file
np.savetxt("base_ca.txt", base_ca.positions)

time_points = [0] ## as a test, ask which are best


for tp in time_points:
    vectors = []
    
    for run in runs:
        for time in times:
            try:
                ref = mda.Universe(f'Repeat_{run}/ns_{time}/equil/md{time}_{tp}.pdb') ## generate eq dumps for new 
                ## ^^ equilibrium dumps
                other = mda.Universe(f'Repeat_{run}/ns_{time}/null/md{time}_{tp}.pdb')
                ## ^^ non equilibrium dumps
                align.alignto(other, ref, select='name CA')
                ca_ref = ref.select_atoms('name CA').positions
                ca_other = other.select_atoms('name CA').positions
                diff = ca_other - ca_ref
                vectors.append(diff)
            except Exception as e:
                print(f"Error processing {run}/{time}: {e}")
                continue

    
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

    adjusted_avg_disp = np.where(np.abs(avg_disp) - 1*se_mag > 0, avg_disp, 0)

    adjusted_avg_vectors = np.where(np.abs(avg_vectors) - 1*se > 0, avg_vectors, 0)

    # Write summary stats
    with open(f"vec_norm_stats/550_1SD_1A5_298_null_{tp}.txt", 'w') as f:
       # Write header
       f.write("#" * 128 + "\n")
       f.write("#                 Average Displacement                Sample Size                                   SE (1SD)                SE (1SD) of              SE (2SD)                SE (2SD) of       #\n")
       f.write("#              --------------------------       ---------------------      Average        --------------------------          Average        --------------------------          Average       #\n")
       f.write("#       CA     x-axis    y-axis    z-axis       x         y         z    Displacement     x-axis    y-axis    z-axis    Displacement vector  x-axis    y-axis    z-axis    Displacement vector #\n")
       f.write("#" * 128 + "\n")

       # Write data
       for i in range(len(avg_vectors)):
           f.write(f"{i+1:>9} {avg_vectors[i][0]:>8.3f} {avg_vectors[i][1]:>8.3f} {avg_vectors[i][2]:>8.3f} ") # removed adjusted vectors
           f.write(f"{len(vectors):>8} {len(vectors):>8} {len(vectors):>8} {avg_disp[i]:>12.3f} ") #removed adjusted displacement
           f.write(f"{se[i][0]:>10.3f} {se[i][1]:>10.3f} {se[i][2]:>10.3f} {se_mag[i]:>12.3f} ")
           f.write(f"{2*se[i][0]:>10.3f} {2*se[i][1]:>10.3f} {2*se[i][2]:>10.3f} {2*se_mag[i]:>12.3f}\n")
           
    # Map the norms of the C-alpha vectors onto the C-alpha atoms in the new Universe object
    for atom, bfactor in zip(ca_universe.atoms, adjusted_avg_disp):
        atom.tempfactor = bfactor

    ## Write the C-alpha-only structure with mapped B-factors
    pdb_filename = f'pdbs/550_1A5_null_vec_norm_{tp}.pdb'
    ca_universe.atoms.write(pdb_filename)
