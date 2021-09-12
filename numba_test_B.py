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
Egen kode !!!

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
        random.seed(1)
        # Setter initialbetingelser for hastighet og posisjon fordelt gaussisk og uniformt
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
            # Bruker Eulers metode for integrasjon
            for i in range(1000):
                T += dt
                for j in range(N):
                    r[j,0,:] = r[j,0,:] + v[j,0,:]*dt
                    if abs(r[j, 0, 0]) >= L/2:     # Tester x-koordinat for vegg-kollisjon
                        v[j, 0, 0] = -v[j, 0, 0]
                        rocket_P += 2*H2_mass*(abs(v[j, 0 ,0])) / (L**2)
                    if abs(r[j, 0, 1]) >= L/2:     # Tester y-koordinat for vegg-kollisjon
                        v[j, 0, 1] = -v[j, 0, 1]
                        rocket_P += 2*H2_mass*(abs(v[j, 0 ,1])) / (L**2)
                    if r[j, 0, 2] >= L/2:     # Tester z-koordinat for tak-kollisjon
                        v[j, 0, 2] = -v[j, 0, 2]
                        rocket_P += 2*H2_mass*(abs(v[j, 0 ,2])) / (L**2)
                    # Tester z-koordinat for kollisjon med gulvet, evt. ut gjennom hullet
                    if r[j, 0, 2] <= -L/2:
                        if abs(r[j, 0, 0]) < nozzle_side/2 and abs(r[j, 0, 1]) < nozzle_side/2:
                            particles_out += 1
                            rocket_p += H2_mass*(-v[j, 0, 2])
                            # v[j, 0, 2] = -v[j, 0 , 2]
                            r[j, 0, :] = [0, 0, L/2 - 1e-12]    # Setter 1e-12 under så ikke koden skal registrere ny posisjon som kollisjon og snu fortegnet
                        else:
                            v[j, 0, 2] = -v[j, 0 , 2]
                            rocket_P += 2*H2_mass*(abs(v[j, 0 ,2])) / (L**2 - nozzle_side**2)

            return v, r, rocket_p/T, rocket_P/(6*T), particles_out/T, T

        self.v, self.r, self.rocket_p, self.rocket_P, self.particles_out, self.T = run()
        analytic_P = self.N*sc.k*self.Temp / self.L**3
        self.mass_loss = self.particles_out * H2.mass
        self.rocket_F = 2*self.rocket_p
        print('-------------------------------------------')
        print(f'| Units per. second:')
        print(f'| Analytic pressure : {analytic_P}Pa')
        print(f'| Numerical pressure : {self.rocket_P}Pa')
        print(f'| Particles out : {self.particles_out}')
        print(f'| Total time : {self.T}s')
        print(f'| rocket_p : {self.rocket_p}kgm/s')
        print(f'| rocket_F : {self.rocket_F}N')
        print(f'| mass loss : {self.mass_loss}kg/s')
        print('-------------------------------------------\n')

if __name__ == '__main__':
    engine = Box(len_box=1e-6, T=3.5e3, num_particles=100000, nozzle_side_len=0.6e-6)  # Må kjøre på 100 000 !!!!!
    engine.simulate()
    print(f'Spacecraft mass : {mission.spacecraft_mass}kg')
    print(f'Spacecraft area : {mission.spacecraft_area}m^2\n')


# Testkjøring resultater:
"""
len_box=1e-6, T=3.5e3, num_particles=100000, nozzle_side_len=0.6e-6

Analytic pressure : 4832.2715Pa
Numerical pressure : 3622.973344556907Pa
Particles out : 78772
Total time : 1.000000000000004e-09s
rocket_p : 1.3734987349945878e-18kgm/s
rocket_F : 2.7469974699891646e-09N
mass loss : 2.636589427799989e-13kg/s
simulate ran in 18.163471937179565 seconds
Spacecraft mass : 1100.0kg
Spacecraft area : 16.0m^2
"""

G = const.G
v_escape = np.sqrt(2*G*system.masses[0]*1.989e30 / (system.radii[0]*1000))
print(f'v_escape = {v_escape}m/s\n')
planet_radius = system.radii[0]*1000
planet_mass = system.masses[0]*1.989e30
spacecraft_mass = mission.spacecraft_mass
initial_fuel = 4000
number_of_boxes = 1.6e13
thrust_force = engine.rocket_F * number_of_boxes
mass_loss_rocket = engine.mass_loss * number_of_boxes
AU = 149597871*1000
speed_factor = AU / (365*24*60*60)
v_initial = system.initial_velocities[:,0]*speed_factor
r_initial = system.initial_positions[:,0]*AU

# @njit(cache=True)
def rocket_boost(dv, fuel, gravity=True):
    if gravity == True:
        gamma = G
        print('gravity is on')
    if gravity == False:
        gamma = 0
        print('gravity is off')
    total_mass = fuel + spacecraft_mass
    print(f'gamma={gamma}')
    fuel_consumption = 0
    time = 9000
    dt = 1e-3
    N = int(np.ceil(time/dt))
    v = np.zeros((N,2))
    r = np.zeros((N,2))
    a = np.zeros((N,2))
    v_esc = np.zeros((N,2))
    r[0,:] = r_initial
    r_vec0 = r[0,:]
    r_norm0 = np.linalg.norm(r_vec0)
    r_unit0 = r_vec0/r_norm0
    v[0,:] = v_initial
    a[0,:] = ((thrust_force - (gamma*planet_mass*total_mass / r_norm0**2)) / total_mass)*r_unit0
    print(f'acceleration vector at t=0 s = {a[0,:]}')
    v_rad = np.dot(v[0,:], r_unit0)
    T = np.zeros(N)
    i = 0
    T_success = 0
    while v_rad <= dv:
        r_vec = r[i,:]
        r_norm = np.linalg.norm(r_vec)
        r_unit = r_vec/r_norm
        v[i+1,:] = v[i,:] + a[i,:]*dt
        r[i+1,:] = r[i,:] + v[i+1]*dt
        T[i+1] = T[i] + dt
        v_esc[i+1,:] = dv*r_unit
        v_rad = np.dot(v[i+1,:], r_unit)
        total_mass -= mass_loss_rocket*dt
        fuel_consumption += mass_loss_rocket*dt
        a[i+1,:] = ((thrust_force - (gamma*planet_mass*total_mass / r_norm**2))) / total_mass * r_unit
        i += 1
        T_success = T[i]
        if total_mass <= spacecraft_mass:
            print('not enough fuel')
            print('r=', r[i,:])
            print('v=', v[i,:])
            return fuel_consumption, total_mass - spacecraft_mass, T[:i-1], r[:i-1,:], v[:i-1,:]
    print('Succesful boost in', T_success, 'seconds')
    return fuel_consumption, total_mass - spacecraft_mass, T[:i-1], r[:i-1,:], v[:i-1,:]

mass_consumed_boost, fuel_left, T, r, v = rocket_boost(v_escape, initial_fuel, gravity=True)
print(f'consumption : {mass_consumed_boost}\nfuel left : {fuel_left}')
print(f'Position:{r[-1,:]}')
print(f'Velocity:{v[-1,:]}')

r_planet = np.array([1.8589133169 - (planet_radius*6.68458712e-9), 0])   # vektoren peker i sentrum av planeten
R = r_planet+r[-1,:]    # vektoren peker fra stjernen på raketten
print(R)

plt.plot(T, r[:,0], label='Gravity on')

mass_consumed_boost, fuel_left, T, r, v = rocket_boost(v_escape, initial_fuel, gravity=False)

print(f'consumption : {mass_consumed_boost}\nfuel left : {fuel_left}')
print(f'Position:{r[-1,:]}')
print(f'Velocity:{v[-1,:]}')

plt.plot(T, r[:,0], label='Gravity off')
plt.legend()
plt.ticklabel_format(style='plain', useOffset=None)
plt.title('Distance from planet in x-direction')
plt.xlabel('Time')
plt.ylabel('Km')
# plt.show()
"""
v_escape = 12249.367264147773m/s

gravity is on
gamma=6.6743e-11
acceleration vector at t=0 s = [8.61803127 0.        ]
Succesful boost in 835.7659999871487 seconds
consumption : 3525.7148795100384
fuel left : 474.28512044418994
Position:[2.78093612e+11 2.54765157e+07]
Velocity:[12246.54613139 30483.39596771]
[2.78093612e+11 2.54765157e+07]
gravity is off
gamma=0
acceleration vector at t=0 s = [8.61803128 0.        ]
Succesful boost in 835.7659999871487 seconds
consumption : 3525.7148795100384
fuel left : 474.28512044418994
Position:[2.78093612e+11 2.54765157e+07]
Velocity:[12246.54613756 30483.39596771]
"""

launch_position = np.array([1.8589133169+(planet_radius*6.68458712e-12), 0])
mission.set_launch_parameters(thrust_force, mass_loss_rocket, initial_fuel+3500, 1800, launch_position, 0)
mission.launch_rocket()
mission.verify_launch_result(R)
