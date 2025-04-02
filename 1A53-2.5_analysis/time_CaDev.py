import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis import distances



repeats = [1,2,3,4,5,6,7,8,9,10] #range(1,11)
start_times = [50,60, 75, 80, 100, 120, 125, 140, 150, 160, 175, 180, 200, 220, 225, 250, 260, 275, 280, 300, 320, 325, 340, 350, 360, 375, 380, 400, 420, 425, 440, 450, 460, 475, 480, 500]

time_points = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]


for tp in time_points:

    devs = []

    for rep in repeats:
        for time in start_times:

            ref = mda.Universe(f'Repeat_{rep}/ns_{time}/equil/md{time}_{tp}.pdb') ## equilibrium
            other = mda.Universe(f'Repeat_{rep}/ns_{time}/md{time}_{tp}.pdb') ##non equilibrium
            
            align.alignto(other, ref, select='name CA')
         
            # Calculate deviation for C-alpha atoms
            ref_ca = ref.select_atoms('name CA')
            other_ca = other.select_atoms('name CA')     
           
            dev = distances.distance_array(ref_ca.positions, other_ca.positions)    
            diag_dev = np.diagonal(dev) * 0.1                    

            devs.append(diag_dev)

            avg_devs = np.mean(devs, axis=0)
             
            std_devs = np.std(devs, axis=0)
            std_errs = std_devs / np.sqrt(len(devs))
    column_name = 'rmsd'       
    with open(f'time_CaDeviations/rmsd_{tp}.txt', 'w') as f:
        f.write(f'{column_name}\n')
        for val in avg_devs:
            
            f.write(f' {val}\n')
     
    with open(f'time_CaDeviations/dev_stats_{tp}.txt', 'w') as f:
      for i in range(len(avg_devs)):
        mean = avg_devs[i]
        std = std_devs[i]
        SE = std_errs[i]
        f.write(f'{mean} {std} {SE}\n')
