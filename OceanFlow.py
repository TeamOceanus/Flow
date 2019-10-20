# -*- coding: utf-8 -*-
from netCDF4 import Dataset
dataset = Dataset('world_oscar_vel_5d2019.nc.gz.nc4')

"""Read the `uf` data to numpy array."""

uf = dataset.variables['uf'][0][0].data
vf = dataset.variables['vf'][0][0].data

lon = dataset.variables['longitude']
lat = dataset.variables['latitude']

import matplotlib.pyplot as plt
import numpy as np
import random as rd

x_axis = np.arange(0, 360, 1)
y_axis = np.arange(139, -1, 1)

factor = 1 # 상수

def findNext(latitude, longitude): # 위도, 경도
    global factor, uf, vf, lon, lat
    lat2 = 139 - (grid(latitude, lat[0]) % 140)
    lon2 = grid(longitude, lon[0]) % 360
    u2 = uf[lat2][lon2]
    v2 = vf[lat2][lon2]
    if np.isnan(u2) or np.isnan(v2):
        return False
    finPhi = latitude + factor * v2
    finTheta = longitude + factor * u2
    if finTheta < 20.5:
        amount = 20.5 - finTheta
        finTheta = 379.5 - amount
    if finTheta > 379.5:
        amount = finTheta - 379.5
        finTheta = 20.5 + amount
    return (finPhi, finTheta) # 위도, 경도

def grid(A, org):
    return int(round(A - 0.5) + 0.5 - org) 
    
def findStart():
    global uf, vf
    longitude = rd.uniform(20.5, 379.5)
    latitude = rd.uniform(-69.5, 69.5)
    Theta = 139 - (grid(latitude, lat[0]) % 140)
    Phi = grid(longitude, lon[0]) % 360
    while np.isnan(uf[Theta][Phi]):
        longitude = rd.uniform(20.5, 379.5)
        latitude = rd.uniform(-69.5, 69.5)
        Theta = 139 - (grid(latitude, lat[0]) % 140)
        Phi = grid(longitude, lon[0]) % 360
    return latitude, longitude # 위도, 경도

def euclidD(A, B):
    return pow((A[0] - B[0])**2 + (A[1] - B[1])**2, 0.5)

plt.figure(figsize=(12.8, 9.6))
plt.quiver(lon[::2], lat[::2], uf[::2, ::2], vf[::2, ::2])
for k in range(1000):
    prevposition = findStart() # lat & lon
    routeX = [prevposition[1]] # lon
    routeY = [prevposition[0]] # lat
    print(f"Start Point : {prevposition}")
    for i in range(1000) :
        position = findNext(prevposition[0], prevposition[1])
        if not position or euclidD(position, prevposition) <= 1e-5:
            break
        routeY.append(position[0]) # lat
        routeX.append(position[1]) # 
        prevposition = position
    plt.scatter(routeX, routeY, 0.1)
plt.show()
print("E")

