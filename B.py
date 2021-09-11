import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission
import numpy as np
import random
import scipy.constants as sc
import matplotlib.pyplot as plt
import time

"""
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
random.seed(1)
#is just a prank bro
class Box:
    def __init__(self, len_box, T, num_particles, nozzle_side_len):
        self.Temp = T; self.N = num_particles
        self.nozzle_side = nozzle_side_len
        self.p_pos = np.zeros((num_particles, 1, 3))
        self.p_vel = np.zeros((num_particles, 1, 3))
        self.L = len_box

    @timer
    def set_initials(self):
        for i in range(self.N):
            for j in range(3):
                self.p_pos[i, 0, j] = random.uniform(-self.L/2, self.L/2)
                self.p_vel[i, 0, j] = random.gauss(0, np.sqrt(sc.k*self.Temp / H2.mass))
    @timer
    def simulate(self):
        v = self.p_vel; r = self.p_pos
        dt = 1e-12
        self.T = 0
        self.particles_out = 0
        self.rocket_p = 0
        self.rocket_F = 0
        self.mass_loss = 0
        for i in range(1000):
            r[:,0,:] = r[:,0,:] + v[:,0,:]*dt
            self.T += dt
            # ax.scatter(r[0,0,0], r[0,0,1], r[0,0,2], s=0.2, color='#1f77b4')
            for j in range(self.N):
                    if abs(r[j, 0, 0]) >= self.L/2:     # Tester x-koordinat for vegg-kollisjon
                        v[j, 0, 0] = -v[j, 0, 0]
                    if abs(r[j, 0, 1]) >= self.L/2:     # Tester y-koordinat for vegg-kollisjon
                        v[j, 0, 1] = -v[j, 0, 1]
                    if r[j, 0, 2] >= self.L/2:     # Tester z-koordinat for tak-kollisjon
                        v[j, 0, 2] = -v[j, 0, 2]
                    # Tester z-koordinat for kollisjon med gulvet, evt. ut gjennom hullet
                    if r[j, 0, 2] <= -self.L/2:
                        if abs(r[j, 0, 0]) < self.nozzle_side/2 and abs(r[j, 0, 1]) < self.nozzle_side/2:
                            self.particles_out += 1
                            self.rocket_p += H2.mass*(-v[j, 0, 2])
                            self.rocket_F += 2*H2.mass*(-v[j, 0, 2]) / dt
                            r[j, 0, :] = [0, 0 , self.L/2]
                        else:
                            v[j, 0, 2] = -v[j, 0 , 2]
        self.mass_loss = self.particles_out * H2.mass / self.T
        print(f'Particles out : {self.particles_out}')
        print(f'Total time : {self.T}')
        print(f'rocket_p : {self.rocket_p}')
        print(f'rocket_F : {self.rocket_F}')
        print(f'mass loss : {self.mass_loss}')

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
engine = Box(len_box=1e-6, T=3e3, num_particles=100000, nozzle_side_len=0.5e-6)  # Må kjøre på 100 000 !!!!!
# engine.set_initials()
# engine.simulate()
print(mission.spacecraft_mass)

# system.print_info()

G = const.G
v_escape = np.sqrt(2*G*system.masses[0]*1.989e30 / (system.radii[0]*1000))
F = mission.spacecraft_mass*v_escape / 600
print(f'v_escape = {v_escape}m/s')
print(f'F = {F}N')
# plt.show()

# Testkjøring resultater:
"""
len_box=1e-6, T=3e3, num_particles=100000, nozzle_side_len=0.5e-6

set_initials ran in 2.1236491203308105 seconds
Particles out : 45958
Total time : 1.000000000000004e-09
rocket_p : 7.333733882078984e-19
rocket_F : 1.4667467764157737e-06
mass loss : 1.5382671116999938e-13
simulate ran in 303.56196689605713 seconds
1100.0
v_escape = 12245.420725509939m/s
F = 22449.93799676822N
"""
print(system.initial_positions[:,0])

pos = (system.initial_positions[0,0] + system.radii[0]*6.68458712e-9, system.initial_positions[1,0])

a = const.G*system.masses[0]*1.989e30/(system.radii[0]*1000)**2
print(a)
mission.set_launch_parameters(6*1.4667467764157737e6, 6*1.5382671116999938e-1, 200000, 2*600, pos, 0)
print(314*1.5382671116999938e-1)
mission.launch_rocket(0.01)
# mission.verify_launch_result()
