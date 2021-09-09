import numpy as np

rocket_thrust_force = 6*1.4667467764157737e6 # N
fuel_consumption = 6*1.5382671116999938e-1 # kg/s
initial_rocket_mass = 1100.0 # kg

# time = 314 # s
# dt = 0.0001 # s
# N = int(np.ceil(time/dt))

def fuel_burn(rocket_thrust_force, fuel_consumption, initial_rocket_mass, speed_boost):
        dv = speed_boost
        
