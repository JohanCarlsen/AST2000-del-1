import numpy as np

rocket_thrust_force = 6*1.4667467764157737e6 # N
fuel_consumption = 6*1.5382671116999938e-1 # kg/s
initial_rocket_mass = 1100.0 # kg

def fuel_burn(rocket_thrust_force, fuel_consumption, initial_rocket_mass, speed_boost):
        fuel_consumed = 0.
        F = rocket_thrust_force
        for i in range(1000):
            m = initial_rocket_mass-fuel_consumption    # mass whith mass loss
            a = F/m # acceleration
