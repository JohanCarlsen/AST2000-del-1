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

"""
Egen kode

Task B, C, D
"""

seed = utils.get_seed('antonabr')
system = SolarSystem(seed)
"""
print('My system has a {:g} solar mass star with a radius of {:g} kilometers.'
      .format(system.star_mass, system.star_radius))

for planet_idx in range(system.number_of_planets):
    print('Planet {:d} is a {} planet with a semi-major axis of {:g} AU.'
          .format(planet_idx, system.types[planet_idx], system.semi_major_axes[planet_idx]))

# print(help(system))
print(system.has_moons)     # True
system.print_info()
"""
mission = SpaceMission(seed)
"""
home_planet_idx = 0 # The home planet always has index 0
print('My mission starts on planet {:d}, which has a radius of {:g} kilometers.'
      .format(home_planet_idx, mission.system.radii[home_planet_idx]))

print('My spacecraft has a mass of {:g} kg and a cross-sectional area of {:g} m^2.'
      .format(mission.spacecraft_mass, mission.spacecraft_area))
"""

# Dekorator for å fange tid kode tar å kjøre (reell tid)
def timer(f):
    def wrapper(*args, **kwargs):
        t1 = time.time()
        result = f(*args, **kwargs)
        t2 = time.time() - t1
        print(f'{f.__name__} ran in {t2} seconds')
        return result

    return wrapper

class Particle:
    def __init__(self, m):
        self.mass = m

H2 = Particle(2*1.6735575e-27)  # H2 particle

class Box:
    def __init__(self, len_box, T, num_particles, nozzle_side_len):
        self.Temp = T; self.N = num_particles
        self.nozzle_side = nozzle_side_len
        self.L = len_box

    @timer
    def simulate(self):
        N = self.N; L = self.L; nozzle_side = self.nozzle_side
        H2_mass = H2.mass; Temp = self.Temp
        dt = 1e-12
        k = sc.k
        r_init = np.zeros((N, 1, 3), dtype=np.float64)
        v_init = np.zeros((N, 1, 3), dtype=np.float64)
        random.seed(2001)
        for i in range(N):
            for j in range(3):
                r_init[i, 0, j] = random.uniform(-L/2, L/2)
                v_init[i, 0, j] = random.gauss(0, np.sqrt(k*Temp / H2_mass))

        @njit
        def run():
            T = 0
            particles_out = 0
            rocket_p = 0
            rocket_F = 0
            rocket_P = 0
            r = r_init.copy()
            v = v_init.copy()
            for i in range(1000):
                T += dt
                # ax.scatter(r[0,0,0], r[0,0,1], r[0,0,2], s=0.2, color='#1f77b4')
                for j in range(N):
                    r[j,0,:] = r[j,0,:] + v[j,0,:]*dt
                    if abs(r[j, 0, 0]) >= L/2:     # Tester x-koordinat for vegg-kollisjon
                        v[j, 0, 0] = -v[j, 0, 0]
                        rocket_P += 2*H2_mass*(-v[j, 0 ,0]) / (L**2)
                    if abs(r[j, 0, 1]) >= L/2:     # Tester y-koordinat for vegg-kollisjon
                        v[j, 0, 1] = -v[j, 0, 1]
                        rocket_P += 2*H2_mass*(-v[j, 0 ,1]) / (L**2)
                    if r[j, 0, 2] >= L/2:     # Tester z-koordinat for tak-kollisjon
                        v[j, 0, 2] = -v[j, 0, 2]
                        rocket_P += 2*H2_mass*(-v[j, 0 ,2]) / (L**2)
                    # Tester z-koordinat for kollisjon med gulvet, evt. ut gjennom hullet
                    if r[j, 0, 2] <= -L/2:
                        if abs(r[j, 0, 0]) < nozzle_side/2 and abs(r[j, 0, 1]) < nozzle_side/2:
                            particles_out += 1
                            rocket_p += H2_mass*(-v[j, 0, 2])
                            rocket_F += 2*H2_mass*(-v[j, 0 ,2]) / dt
                            # v[j, 0, 2] = -v[j, 0 , 2]
                            r[j, 0, :] = [0, 0, L/2]
                        else:
                            v[j, 0, 2] = -v[j, 0 , 2]
                            rocket_P += 2*H2_mass*(-v[j, 0 ,2]) / (L**2 - nozzle_side**2)

            return v, r, rocket_p, rocket_F, rocket_P, particles_out, T

        self.v, self.r, self.rocket_p, self.rocket_F, self.rocket_P, self.particles_out, self.T = run()
        analytic_P = self.N*sc.k*self.Temp / 1e-18
        self.mass_loss = self.particles_out * H2.mass / self.T
        print(analytic_P)
        print(self.rocket_P)
        print(f'Particles out : {self.particles_out}')
        print(f'Total time : {self.T}')
        print(f'rocket_p : {self.rocket_p}')
        print(f'rocket_F : {self.rocket_F}')
        print(f'mass loss : {self.mass_loss}')

if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    engine = Box(len_box=1e-6, T=3e3, num_particles=100000, nozzle_side_len=0.5e-6)  # Må kjøre på 100 000 !!!!!
    engine.simulate()

    print(mission.spacecraft_mass)

# system.print_info()

G = const.G
v_escape = np.sqrt(2*G*system.masses[0]*1.989e30 / (system.radii[0]*1000))
print(mission.spacecraft_area)
print(f'v_escape = {v_escape}m/s')

# Testkjøring resultater:
"""
len_box=1e-6, T=3e3, num_particles=100000, nozzle_side_len=0.5e-6

Particles out : 46314
Total time : 1.000000000000004e-09
rocket_p : 7.487250746872403e-19
rocket_F : 1.4974501493744745e-09
mass loss : 1.5501828410999938e-13
simulate ran in 17.815258979797363 seconds
1100.0
v_escape = 12249.367264147773m/s
F = 22457.17331760425N
"""
planet_radius = system.radii[0]*1000
planet_mass = system.masses[0]*1.989e30
spacecraft_mass = mission.spacecraft_mass
initial_fuel = 50000
number_of_boxes = 1.6e13
thrust_force = engine.rocket_F * number_of_boxes
mass_loss_rocket = engine.mass_loss * number_of_boxes
G = const.G
AU = 149597871*1000
speed_factor = AU / (365*24*60*60)
v_initial = system.initial_velocities[:,0]*speed_factor
r_initial = system.initial_positions[:,0]*AU
print(v_initial)

@njit
def rocket_boost(dv):
    total_mass = initial_fuel + spacecraft_mass
    a = thrust_force / total_mass
    v = 0
    dt = 1e-3
    while v <= dv:
        v = v + a*dt
        total_mass -= mass_loss_rocket*dt
        a = thrust_force / total_mass
    consumption = initial_fuel + spacecraft_mass - total_mass
    return consumption

mass_left = rocket_boost(1000)
print(mass_left)

"""
1100.0      # spacecraft_mass
16.0        # spacecraft_area
v_escape = 12249.367264147773m/s
Launch success at time 424.02000000001453
51034.56860757253
[[0.         0.        ]
 [1.73239477 0.        ]
 [3.46478985 0.        ]
 ...
 [0.         0.        ]
 [0.         0.        ]
 [0.         0.        ]]
"""
# print(system.initial_positions[:,0])
#
# pos = (system.initial_positions[0,0] + system.radii[0]*6.68458712e-9, system.initial_positions[1,0])
#
# a = const.G*system.masses[0]*1.989e30/(system.radii[0]*1000)**2
# print(a)
# mission.set_launch_parameters(54*engine.rocket_F*1e12, 54*engine.mass_loss*1e12, 5000, 2*600, pos, 0)
# mission.launch_rocket(0.01)
# mission.verify_launch_result()
