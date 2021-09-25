import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
from ast2000tools.space_mission import SpaceMission

seed = utils.get_seed('antonabr')
system = SolarSystem(seed)
mission = SpaceMission(seed)

# print(mission.spacecraft_mass)
# system.print_informations()
# print(system.star_temperature)
# print(system.aphelion_angles)
# print(system.eccentricities)
# print(system.semi_major_axes)

'''
Planets in order of closest to star:
0, 1, 6, 7, 2, 5, 3, 4
'''


'''
Star temperature: 7985.019468972192 K
Sun temperature: 5778 K
'''

'''
Aphelion angels: [0, 4.15430044, 2.36900623, 0.84972349, 0.72126164, 4.03998827, 4.72664666, 1.41925568] in radians

Eccentricities: [0.00562084, 0.01669338, 0.04303139, 0.08387873, 0.02863436, 0.03421009, 0.05231877, 0.02263108]

Semi-major axes: [ 1.84852307,  2.79833111,  7.75328387, 13.5615841,  21.12461887, 10.00416126, 3.63010736,  4.77438331] in AU
'''

'''
Spacecraft mass: 1100.0 kg

Information about the solar system with seed 12497:

Number of planets: 8
Star surface temperature: 7985.02 K
Star radius: 1.05932e+06 km
Star mass: 1.95532 solar masses

Individual planet information. Masses in units of m_sun, radii in km,
atmospheric densities in kg/m^3, rotational periods in Earth days.
Planet |      Mass      |     Radius     |   Atm. dens.   |  Rot. period   |
     0 | 4.29832228e-06 |   7605.7520967 |    1.078127856 |     1.14007756 |
     1 | 4.78574031e-06 |   7802.4617203 |    0.931610754 |     1.23165684 |
     2 | 8.63339133e-08 |   1956.2179196 |    1.035659804 |     8.48022285 |
     3 | 6.12889201e-05 |  31720.5340165 |   22.890517718 |     0.23785162 |
     4 | 5.90408327e-09 |    924.4812718 |    1.372329079 |    18.40430715 |
     5 | 9.71237425e-05 |  35956.4210889 |   23.885096642 |     0.85762836 |
     6 | 8.29175119e-08 |   1861.7435943 |    1.048075284 |    11.14961304 |
     7 |  0.00474102794 | 121492.0095111 |   19.407881947 |     0.82165937 |

Individual planet initial positions (in AU) and velocities (in AU/yr).
Planet |       x        |       y        |       vx       |       vy       |
     0 |   1.8589133169 |   0.0000000000 |   0.0000000000 |   6.4259156226 |
     1 |   0.4490322100 |  -2.7973842411 |   5.1121279534 |   0.8789685580 |
     2 |   2.1818413255 |  -7.1311951343 |   3.1149297504 |   1.0213346778 |
     3 |  -9.2883580520 |   9.9006298663 |  -1.5952694828 |  -1.7707095761 |
     4 |  21.1961564258 |   4.3905636813 |  -0.3517341360 |   1.8314963881 |
     5 | -10.0336061581 |   1.6165049612 |  -0.5164741841 |  -2.6848009086 |
     6 |   2.3772472238 |   2.5529518086 |  -3.6209823467 |   3.1433858518 |
     7 |  -4.4606581566 |  -1.5492522085 |   1.4095527892 |  -3.8131016460 |
'''
