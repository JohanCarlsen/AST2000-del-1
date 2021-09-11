import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
import numpy as np
import random
import scipy.constants as sc
import matplotlib.pyplot as plt
import time
from numba import njit
from numba_test_B import *

seed = utils.get_seed('antonabr')
system = SolarSystem(seed)
mission = SpaceMission(seed)

planet_radius = system.radii[0]*1000
planet_mass = system.masses[0]*1.989e30
spacecraft_mass = mission.spacecraft_mass
initial_fuel = 50000
number_of_boxes = 1e13
thrust_force =




def rocket_launch(time, N):
    N = N
    dt = time / N
    T = 0
    AU = 149597871
    theta = 0
    e_r = np.array([np.cos(theta), np.sin(theta)])
    v = np.zeros((time, 2))
    r = np.zeros((time, 2))
    v[0] = [0,0]
    r[0] = [planet_radius, 0]
    planet_mass = planet_mass
    system_M = initial_fuel + spacecraft_mass
    F = 2.0612e7
    a = (F - const.G*planet_mass*r[0]/np.linalg.norm(r[0])**3)/system_M
    for i in range(time):
        T += dt
        v[i+1] = v[i] + a*dt
        r[i+1] = r[i] + v[i+1]*dt
        system_M -= engine.mass_loss*1e13*dt
        a = (F - const.G*planet_mass*r[i+1]/np.linalg.norm(r[i+1])**3)/system_M
        if np.linalg.norm(v[i+1]) >= v_escape:
            print(f'Launch success at time {T}')
            return v, r
    print('No success')
    return v, r

v, r = rocket_launch(600, 10000)
print(v)
