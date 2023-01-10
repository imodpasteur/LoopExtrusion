import numpy as np
L=600
# Select random initial loading of extrusion
start = np.random.randint(277,325)
anchor_left=276
anchor_right=325
print("variable j equal", start-1)
print("variable i equal", start+1)
print("variable L equal", L)
print("variable Left_anch equal", anchor_left)
print("variable Right_anch equal", anchor_right)

import random
seed_velocity=random.randint(1000,10000000)
seed_langevin=random.randint(1000,10000000)
print("variable v equal", seed_velocity)
print( "variable l equal", seed_langevin)

L = []
for i in range(6):
    L.append(np.random.randint(200,100000))
    print("variable rand_%s equal"%(i), L[i])

# Draws the length of the closed state from an exponential distribution
mean_exp = 1400  # Mean number of frames in closed state if distribution is not cut
total_nbr_closed_max = 2500   # max number of frames in the closed state
for i in range(10000):
    h = int(np.random.exponential(mean_exp))
    if h < total_nbr_closed_max:
        death = h
        break
frames_closed = death*1000
print("variable frc equal", frames_closed)

size_nuc = 18
print("variable Radius equal", size_nuc)
