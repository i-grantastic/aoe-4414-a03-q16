# lh_to_rz.py
#
# Usage: python3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km
#   converts SEZ coordinates to ECEF coordinates

# Parameters:
#   o_lat_deg: origin latitude in degrees
#   o_lon_deg: origin longitude in degrees
#   o_hae_km: origin height above the ellipsoid in km
#   s_km: s-component
#   e_km: e-component
#   z_km: z-component

# Output:
#   ECEF coordinates
#
# Written by Grant Chapman
# Other contributors: None

# import Python modules
import math # math module
import sys # argv
import numpy as np

# constants
R_E_KM = 6378.137
E_E    = 0.081819221456

## calculate demoninator
# (eccentricity, latitude in radians)
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-ecc**2.0 * math.sin(lat_rad)**2.0)

# initialize script arguments
o_lat_deg = float('nan')
o_lon_deg = float('nan')
o_hae_km  = 0.0
s_km = float('nan')
e_km = float('nan')
z_km = float('nan')

# parse script arguments
if len(sys.argv) == 7:
  o_lat_deg = float(sys.argv[1])
  o_lon_deg = float(sys.argv[2])
  o_hae_km  = float(sys.argv[3])
  s_km = float(sys.argv[4])
  e_km = float(sys.argv[5])
  z_km = float(sys.argv[6])
else:
  print(\
    'Usage: '\
    'python3 sez_to_ecef.py o_lat_deg o_lon_deg o_hae_km s_km e_km z_km'\
  )
  exit()

### script below this line ###

# latitude and longitude to radians
o_lat_rad = o_lat_deg*math.pi/180.0
o_lon_rad = o_lon_deg*math.pi/180.0

# rotation matrices
r_y = [[math.sin(o_lat_rad), 0, math.cos(o_lat_rad)],
       [0, 1, 0],
       [-math.cos(o_lat_rad), 0, math.sin(o_lat_rad)]]
r_z = [[math.cos(o_lon_rad), -math.sin(o_lon_rad), 0],
       [math.sin(o_lon_rad), math.cos(o_lon_rad), 0],
       [0, 0, 1]]

# SEZ vector
sez_vector = [s_km, e_km, z_km]

# first rotation
first_rotation = np.matmul(r_y, sez_vector)

# second rotation
second_rotation = np.matmul(r_z, first_rotation)

# finding c_E and s_E
c_E = R_E_KM / calc_denom(E_E, o_lat_rad)
s_E = R_E_KM*(1-E_E**2) / calc_denom(E_E, o_lat_rad)

# finding ECEF components
r_x_km = (c_E+o_hae_km)*math.cos(o_lat_rad)*math.cos(o_lon_rad)
r_y_km = (c_E+o_hae_km)*math.cos(o_lat_rad)*math.sin(o_lon_rad)
r_z_km = (s_E+o_hae_km)*math.sin(o_lat_rad)

# finding ECEF vector of SEZ target
ecef_x_km = second_rotation[0] + r_x_km
ecef_y_km = second_rotation[1] + r_y_km
ecef_z_km = second_rotation[2] + r_z_km

# print
print(ecef_x_km)
print(ecef_y_km)
print(ecef_z_km)