from libraries import *
from variables import *

# Defining a class for the simulation parameters
# Class for holding info about the simulation params
class SIM_config:
    def __init__(self,unit,N,dt,wavelength0):
        self.number_of_points=N
        speed_of_lightdictionary = {
            'SI': speed_of_light,                                                                                   # Speed of light value in SI units to calculate with                                                                                              # Speed of light value in atomic units to calculate with
            'OTHER' : 300,                                                                                          # Speed of light value in nm / fs units to calculate with
            'DIMENSIONLESS': 1
        }
        self.speed_of_light = speed_of_lightdictionary[unit]
        wavelength0dictionary = {
            'SI': wavelength0,                                                                                      # Pulse central wavelengt [m]                                                                       # Central wavelength value in atomic units to calculate with
            'OTHER' : wavelength0 * 1e9,
            'DIMENSIONLESS': 1                                                                                                                                                          # Central wavelength value in nm units to calculate with
        }
        self.wavelength0 = wavelength0dictionary[unit]
        timestepdictionary = {
            'SI': dt,                                                                                               # Time resolution [s]                                                                                  # Time resolution value in atomic units to calculate with
            'OTHER': dt * 1e15,
            'DIMENSIONLESS': dt / T                                                                                 # Dimensionless time step                                                                                                                                                                           # Time resolution value in FS units to calculate with
        }
        self.time_step = timestepdictionary[unit]
        t = np.linspace(0,N * dt,N)                                                                                 # Time step array
        t = t - np.mean(t)
        tdictionary = {
            'SI': t,
            'OTHER': t * 1e15,
            'DIMENSIONLESS': t / T                                                                                  # Dimensionless time array
        }
        self.t = tdictionary[unit]
        frequency0 = speed_of_light / wavelength0
        frequency0dictionary = {
            'SI': frequency0,                                                                                       # Pulse central frequency [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm                                                                     # Central frequency value in atomic units to calculate with
            'OTHER': frequency0 / 1e15,                                                                             # Central frequency value in pHz units to calculate with
            'DIMENSIONLESS': 1
        }
        self.frequency0 = frequency0dictionary[unit]
        f = fftshift(fftfreq(N,d=dt))                                                                               # Frequency step array
        fdictionary = {
            'SI': f,
            'OTHER': f / 1e15,
            'DIMENSIONLESS': f / frequency0
        }
        self.f = fdictionary[unit]
        f_rel = f + frequency0                                                                                      # Relative frequency step array 
        f_reldictionary = {
            'SI': f_rel,                                                                                                                                                                           
            'OTHER': f_rel / 1e15,
            'DIMENSIONLESS': f_rel / frequency0
        }
        self.f_rel = f_reldictionary[unit]
        wavelength = speed_of_light / f                                                                             # Wavelength step array 
        wavelengthdictionary = {
            'SI': wavelength,
            'OTHER': wavelength * 1e9,
            'DIMENSIONLESS': wavelength / wavelength0
        }
        self.wavelength = wavelengthdictionary[unit]
        wavelength_rel = wavelength + wavelength0                                                                   # Relative wavelength step array 
        wavelength_reldictionary = {
            'SI': wavelength_rel,
            'OTHER': wavelength_rel * 1e9,
            'DIMENSIONLESS': wavelength_rel / wavelength0
        }
        self.wavelength_rel = wavelength_reldictionary[unit]
        durationdictionary = {
            'SI': duration,                                                                                          # Pulse duration (FWHM) [s]                                                                           
            'OTHER': duration * 1e15,                                                                                # Pulse duration (FWHM) in fs to calculate with
            'DIMENSIONLESS': duration / T
        }
        self.duration = durationdictionary[unit]                                                                     # Repetition frequency value in atomic units to calculate with
        pulse_energy = average_power / repetition_frequency                                                          # Pulse energy [J]
        peak_power = pulse_energy / duration                                                                         # Pulse peak power [W]                                                              
        amplitudedictionary = {
            'SI': 1,                                                                                                 # Electrical field strength amplitude in units of sqrt(W)                                                                          
            'OTHER': 1,
            'DIMENSIONLESS': np.sqrt(peak_power)                                                                    
        }
        self.amplitude = amplitudedictionary[unit]

# Class for holding info about the plot params
class PLOT_config:
    def __init__(self,unit,N,dt,wavelength0,duration):
        self.number_of_points=N
        speed_of_lightdictionary = {
            'SI': speed_of_light,                                                                                   # Speed of light value in SI units to calculate with                                                                                             # Speed of light value in atomic units to calculate with
            'OTHER' : 300,                                                                                          # Speed of light value in nm / fs units to calculate with
            'DIMENSIONLESS': 1
        }
        self.speed_of_light = speed_of_lightdictionary[unit]
        wavelength0dictionary = {
            'SI': wavelength0,                                                                                      # Pulse central wavelengt [m]                                                                       # Central wavelength value in atomic units to calculate with
            'OTHER' : wavelength0 * 1e9,                                                                            # Central wavelength value in nm units to calculate with
            'DIMENSIONLESS': 1
        }
        self.wavelength0 = wavelength0dictionary[unit] 
        t = np.linspace(0,N * dt,N)                                                                                 # Time step array
        t = t - np.mean(t)
        tdictionary = {
            'SI': t,
            'OTHER': t * 1e15,
            'DIMENSIONLESS': t / T                                                                                  # Dimensionless time array
        }
        frequency0 = speed_of_light / wavelength0
        frequency0dictionary = {
            'SI': frequency0,                                                                                       # Pulse central frequency [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm                                                                     # Central frequency value in atomic units to calculate with
            'OTHER': frequency0 / 1e15,                                                                             # Central frequency value in pHz units to calculate with
            'DIMENSIONLESS': 1
        }
        self.frequency0 = frequency0dictionary[unit]
        self.t = tdictionary[unit]
        f = fftshift(fftfreq(N,d=dt))                                                                               # Frequency step array
        fdictionary = {
            'SI': f,
            'OTHER': f / 1e15,
            'DIMENSIONLESS': f / frequency0
        }
        self.f = fdictionary[unit]
        f_rel = f + frequency0                                                                                      # Relative frequency step array 
        f_reldictionary = {
            'SI': f_rel,                                                                                                                                                                           
            'OTHER': f_rel / 1e15,
            'DIMENSIONLESS': f_rel / frequency0
        }
        self.f_rel = f_reldictionary[unit]
        wavelength = speed_of_light / f                                                                             # Wavelength step array 
        wavelengthdictionary = {
            'SI': wavelength,
            'OTHER': wavelength * 1e9,
            'DIMENSIONLESS': wavelength / wavelength0
        }
        self.wavelength = wavelengthdictionary[unit]
        wavelength_rel = wavelength + wavelength0                                                                   # Relative wavelength step array 
        wavelength_reldictionary = {
            'SI': wavelength_rel,
            'OTHER': wavelength_rel * 1e9,
            'DIMENSIONLESS': wavelength_rel / wavelength0
        }
        self.wavelength_rel = wavelength_reldictionary[unit]
        durationdictionary = {
            'SI': duration,                                                                                         # Pulse duration (FWHM) [s]                                                                            # Pulse duration (FWHM) in atomic units to calculate with  
            'OTHER': duration * 1e15,
            'DIMENSIONLESS': duration / T                                                                           # Dimensionless time duration                                                                                                                                                                # Pulse duration (FWHM) in fs to calculate with
        }
        self.duration = durationdictionary[unit]

# Class for holding info about the fiber
class Fiber_config:
    def __init__(self,unit,nsteps,L,gamma,beta2,alpha_dB_per_m):
        self.nsteps=nsteps
        self.ntraces=nsteps+1                                                                                       # NOTE: If we want to do 100 steps, we will get 101 calculated pulses
        #dz = L / nsteps                                                                                            # Spatial resolution [m]
        dz = 1*1e-6
        spacestepdictionary = {
            'SI': dz,                                                                                                                                   
            'DIMENSIONLESS' : dz / Z
        }
        self.dz = spacestepdictionary[unit]                                                                         # Spatial step array
        #zlocs = np.linspace(0,L,nsteps)                                                                            # Locations of each calculated pulse                                                                                         
        zlocs = np.linspace(0,dz*nsteps,nsteps)
        zlocsdictionary = {
            'SI': zlocs,                                                    
            'DIMENSIONLESS': zlocs / Z                                      
        } 
        self.zlocs_array = zlocsdictionary[unit]
        beta2dictionary = {
            'SI': beta2,
            'DIMENSIONLESS': 1
        }              
        self.beta2 = beta2dictionary[unit]                                                                          # GVD [s^2/m] of fused silica @ 800nm
        gammadictionary = {
            'SI': gamma,
            'DIMENSIONLESS': 1
        } 
        self.gamma = gammadictionary[unit]                                                                          # Nonlinear parameter [1/(W*m)]
        alpha_dB_per_mdictionary = {
            'SI': alpha_dB_per_m,                                                                                   # Power attenuation coeff [dB/m]. Usual value @ 1550 nm is 0.2 dB/km                                                                                                                     
            'DIMENSIONLESS': alpha_dB_per_m * Z                                      
        } 
        self.alpha_dB_per_m = alpha_dB_per_mdictionary[unit]                   
        alpha_Np_per_m = alpha_dB_per_m*np.log(10)/10.0                                                             # Power attenuation coeff [Nepers/km]
        alpha_Np_per_mdictionary = {
            'SI': alpha_Np_per_m,
            'DIMENSIONLESS': alpha_Np_per_m * Z
        }
        self.alpha_Np_per_m = alpha_Np_per_mdictionary[unit]                                                        # NOTE: The loss coefficient is usually specified in dB/km, but Nepers/km is more useful for calculations