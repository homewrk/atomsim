from neutral_atom_imaging_simulation.Camera import EMCCDCamera
from neutral_atom_imaging_simulation.Experiment import TweezerArray
from neutral_atom_imaging_simulation.ImageGenerator import ImageGenerator
import matplotlib.pyplot as plt
import math
import scipy.constants as constants
import numpy as np

def get_ixon_ultra_888_camera() -> EMCCDCamera:
    cam = EMCCDCamera(resolution=(1024, 1024),
                    dark_current_rate=0.00011, 
                    cic_chance=0.00114, #value of inducing a current during analog to digital conversion
                    quantum_efficiency=0.86,  # at 780 nm, probability of a photon being turned into an electron
                    numerical_aperture=0.75,
                    physical_pixel_size=13.,
                    magnification=62.5, #two lenses on the camera, is a ratio of the two focal lengths
                    bias_clamp=100, #basically biases it so the output is only positive numbers
                    preampgain=5.29,
                    scic_chance=0,
                    number_gain_reg=10,
                    p0=1,
                    exposure_time=1E-3,
                      binning=1,
                      )

    return cam

def get_tweezer_array():
    twz = TweezerArray(stray_light_rate=10,
                 imaging_wavelength=0.780,
                 scattering_rate=6.E6,
                 survival_probability=0.99,
                 fill_rate=1)
    return twz

def check_image_generation():
    cam = get_ixon_ultra_888_camera()
    twz = get_tweezer_array()
    twz.set_atom_sites_camera_space([(0.5, 0.5), (0.75, 0.75), (0.1, 0.1)])
    imgen = ImageGenerator()
    imgen.set_camera(cam)
    imgen.set_experiment(twz)
    im = imgen.create_image()
    fig, ax = plt.subplots(1, 1)
    cax = ax.imshow(im[0], vmin=0)
    fig.colorbar(cax)
    plt.show()

    
def intensity_of_atom(atom_location, xx, yy):
    x = atom_location[0]
    y = atom_location[1]
    z = np.linspace(0, atom_location[2], 250)
    # Numerical aperature (NA) is 0.75
    # wavelength is 780nm
    # 6 mil photons / second coming out of the atom and probably like 15% of that gets to the objective
    wavelength = 780 * (10 ** -9)
    numerical_aperture = 0.75
    refractive_index = 1
    radial_distance = np.sqrt(np.power(xx - x, 2) + np.power(yy - y, 2))
    # derivied by relating the numerical aperature --> raleigh range and gaussian beam waist --> raleigh range, setting the equations equal to each other
    # and solving for the waist radius given numerical aperature

    waist_radius = wavelength / (np.pi * numerical_aperture)
    rayleigh_range = (np.pi * (np.power(waist_radius, 2)) * refractive_index) / wavelength

    radius_at_z = waist_radius * (np.sqrt(1 + np.power((z / rayleigh_range), 2)))

    frequency = constants.speed_of_light / wavelength
    #realistically only get 15%
    power_from_atom = (6000000 * constants.Planck * frequency) * 0.15
    initial_intensity = (2 * power_from_atom) / (math.pi * (waist_radius ** 2)) 
    zz = initial_intensity * np.power((np.power((waist_radius / radius_at_z), 2)), 2) * np.power(np.e, (-2 * (np.power(radial_distance, 2)) / (np.power(radius_at_z, 2))))
    return zz
    

if __name__ == "__main__":
    atoms = [
        (4, 1 , 1),
        (5, 5, 1),
        (6, 4, 1)
    ]
    x_axis = np.linspace(0, 10, 250)
    y_axis = np.linspace(0, 10, 250)

    xx, yy = np.meshgrid(x_axis, y_axis)
    totalintensity = 0
    for atom in atoms:
        zz = intensity_of_atom(atom, xx, yy)
        totalintensity += zz
    plt.imshow(totalintensity)
    plt.show()
