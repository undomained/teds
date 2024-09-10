# Orbit sim notes
import numpy as np


# Calculate gm parameters for new scene
#   input: new lat/lat of scene
#   output: new ltan, start_time, end_time

# note: manual tweaking needed to cover the powerplant plume
# note2: Matimba a bit off, something with southern hemisphere

######
# input
#######
new_lon = 19.828912600536132
new_lat = 51.44156094162394
#####

###
#  Janschwalde input, move relative to this
lon = 14.29
lat = 51.79
ltan = 13.823
time_start = 13.73
time_end = 13.80


delta_lon = new_lon - lon
delta_lat = new_lat - lat



# constants
h = 500.0 #km
Re = 6378.137 #km
mu_E = 3.986005E14 #m^3 s^-2
sidereal_day = 86164.0905 # s
solar_day = 86400.0 # s

# radius orbit
R = (h+Re)*1e3 # m

# orbital period
T = 2*np.pi*np.sqrt(R**3 / mu_E) # s


# sun-synchronous inclination

i  = np.arccos(-1 * (R / (12352.0*1e3))**(7./2.))  # radians


# rotational motion

v_sat = 360/T # [deg/s]

v_earth = 360/sidereal_day # [deg/s]


# rate of sub-satellite point lat/lon
# only on day-side, for i > 90

# lon_rate = - (v_sat*np.sin(i-np.pi/2) + v_earth)  # deg / s
lon_rate = - (v_sat*np.sin(i-np.pi/2) + v_earth)  # deg / s

lat_rate = v_sat*np.cos(i-np.pi/2)  # deg / s


# ltan: local time ascending note. ltan changes with lon
lt_rate = solar_day/360 # s / deg

delta_time = delta_lat/lat_rate

# cannot go back in time in gm
if delta_time < 0 :
    delta_time = T + delta_time
    delta_ltan = lt_rate*delta_lon + delta_time

else:

    delta_ltan = lt_rate*(delta_lon-delta_time*lon_rate)




ltan_new = ltan + delta_ltan/3600 
time_start_new = time_start + delta_time/60
time_end_new = time_end + delta_time/60


print(f'new ltan (h)   = {ltan_new:.3f}')
print(f'new time start (min) = {time_start_new:.2f}')
print(f'end time start (min) = {time_end_new:.2f}')

breakpoint()