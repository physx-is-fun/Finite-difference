from libraries import *

# constants
pi = np.pi
speed_of_light=3*1e8 # Speed of light [m/s]

# Units
simunit = 'DIMENSIONLESS' # SI or OTHER or DIMENSIONLESS
fiberunit = 'DIMENSIONLESS' # SI or DIMENSIONLESS
plotunit = 'DIMENSIONLESS' # SI or OTHER or DIMENSIONLESS

# Labels
timedictionary = {
    'SI': '[s]',
    'OTHER': '[fs]',
    'DIMENSIONLESS': '[arbitrary unit]'
}
frequencydictionary = {
    'SI': '[Hz]',
    'OTHER': '[PHz]',
    'DIMENSIONLESS': '[arbitrary unit]'
}
wavelengthdictionary = {
    'SI': '[m]',
    'OTHER': '[nm]',
    'DIMENSIONLESS': '[arbitrary unit]'
}

distancedictionary = {
    'SI': '[m]',
    'OTHER': '[m]',
    'DIMENSIONLESS': '[arbitrary unit]'
}

# Defining parameters for the simulation
# Initialize Gaussian pulse parameters (OCTAVIUS-85M-HP from THORLABS) https://www.thorlabs.com/thorproduct.cfm?partnumber=OCTAVIUS-85M-HP
wavelength0=800*1e-9                                        # Pulse central wavelengt [m]
frequency0=speed_of_light/wavelength0                       # Pulse central frequency [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm
duration=8*1e-15                                            # Pulse duration in FWHM [s]
repetition_frequency=85*1e6                                 # Pulse repetition frequency [Hz]
average_power=600*1e-3                                      # Pulse average power [W]
pulse_energy=average_power/repetition_frequency             # Pulse energy [J]
peak_power=pulse_energy/duration                            # Pulse peak power [W]
amplitude = np.sqrt(peak_power)                             # Electrical field strength amplitude in units of sqrt(W)
N= 2**15 #2**15                                             # Number of points                                                    
dt = 0.1*1e-15                                              # Time resolution [s]
alpha = 0.5                                                 # fractional order

# Defining the parameters of the fiber
Length=1*1e-3 #5*1e-2                                                              # Fiber length [m]
nsteps=2**10                                                                        # Number of steps we divide the fiber into
effective_mode_diameter=5e-6                                                       # Effective mode diameter [um] from https://www.thorlabs.com/thorproduct.cfm?partnumber=780HP
effective_mode_area=(pi/4)*effective_mode_diameter**2                              # Effective mode area [m^2]
nonlinear_refractive_index=2.7*1e-20                                               # Nonlinear refractive index [m^2/W] of fused silica @ 800 nm from https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-26-37940&id=424534
gammaconstant=(2*pi*nonlinear_refractive_index)/(wavelength0*effective_mode_area)  # Nonlinear parameter [1/(W*m)]
beta2=36.16                                                                        # GVD in fs^2/mm (units typically used when referring to beta2) of fused silica @ 800nm from https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
beta2*=(1e-27)                                                                     # Convert GVD to s^2/m so everything is in SI units of fused silica @ 800nm
alpha_dB_per_m=0.2*1e-3                                                            # Power attenuation coeff in decibel per m. Usual value @ 1550 nm is 0.2 dB/km
# NOTE: beta2>0 is normal dispersion with red light pulling ahead, causing a negative leading chirp
# NOTE: beta2<0 is anomalous dispersion with blue light pulling ahead, causing a positive leading chirp

# Some useful parameters
nonlinear_length=1/(gammaconstant*peak_power)
dispersion_length=(duration**2)/(np.abs(beta2))

# Characteristic time and distance scales
T = np.sqrt(beta2)                                                                 # Characteristic time scale
Z = gammaconstant * peak_power                                                     # Characteristic propagation distance