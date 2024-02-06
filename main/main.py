from neutral_atom_imaging_simulation.Camera import EMCCDCamera
from neutral_atom_imaging_simulation.Experiment import TweezerArray
from neutral_atom_imaging_simulation.ImageGenerator import ImageGenerator
import matplotlib.pyplot as plt

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


if __name__ == "__main__":
    check_image_generation()
