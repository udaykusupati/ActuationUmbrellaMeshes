import os
from images import generate_stresses_1D, generate_heightsEnergies_1D
from tools import img_to_gif_stress_1D,\
                  img_to_gif_heights_1D,\
                  img_to_gif_energies,\
                  gif_to_img_duration

def compare(paths, name, nb_steps, stresses=['VonMises'], gif_duration=5, loop=0, deployments=None):
    if len(paths)<2 : return
    head, tail = os.path.split(os.path.normpath(paths[0]))
    save_dir = f'{head}/compare_{name}'
    
    if deployments==None: deployments = ['linear']*len(paths)
    
    duration = gif_to_img_duration(gif_duration, nb_steps)
    
    # heights and energies
    generate_heightsEnergies_1D(paths, deployments, save_dir)
    img_to_gif_heights_1D(save_dir, duration=duration, loop=loop)
    img_to_gif_energies(save_dir, duration=duration, loop=loop)
    
    # stresses
    for stress in stresses:
        generate_stresses_1D(paths, deployments, save_dir, stress_type=stress)
        img_to_gif_stress_1D(save_dir, stress, duration=duration, loop=loop)
        
        