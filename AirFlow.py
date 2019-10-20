# -*- coding: utf-8 -*-
from netCDF4 import Dataset
import numpy as np
dataset = Dataset('pentad_20111227_v11l35flk.nc.gz.nc4')

"""Read the `uf` data to numpy array."""

uf = dataset.variables['uwnd'][0].data * 0.00152597204
vf = dataset.variables['vwnd'][0].data * 0.00152597204

uf[uf < -50] = np.nan
uf[uf > 50] = np.nan

lon = dataset.variables['lon']
lat = dataset.variables['lat']

vf[vf < -50] = np.nan
vf[vf > 50] = np.nan

import matplotlib.pyplot as plt
import random as rd

factor = 0.1 # 상수

def findNext(latitude, longitude): # 위도, 경도
    global factor, uf, vf, lon, lat
    lat2 = grid(latitude, lat[0])  % len(lat)
    lon2 = grid(longitude, lon[0])  % len(lon)
    u2 = uf[lat2][lon2]
    v2 = vf[lat2][lon2]
    if np.isnan(u2) or np.isnan(v2):
        return False
    finPhi = latitude + factor * v2
    finTheta = longitude + factor * u2
    if finTheta < lon[0]:
        amount = lon[0] - finTheta
        finTheta = lon[-1] - amount
    if finTheta > lon[-1]:
        amount = finTheta - lon[-1]
        finTheta = lon[0] + amount
    return (finPhi, finTheta) # 위도, 경도

def grid(A, org):
    return 4 * int(round(4*(A % 1)) / 4 + np.floor(A) - org)
   
def findStart():
    global uf, vf
    longitude = rd.uniform(lon[0], lon[-1])
    latitude = rd.uniform(lat[0], lat[-1])
    Theta = grid(latitude, lat[0]) % len(lat)
    Phi = grid(longitude, lon[0]) % len(lon)
    while np.isnan(uf[Theta][Phi]):
        longitude = rd.uniform(lon[0], lon[-1])
        latitude = rd.uniform(lat[0], lat[-1])
        Theta = grid(latitude, lat[0]) % len(lat)
        Phi = grid(longitude, lon[0]) % len(lon)
    return latitude, longitude # 위도, 경도

def euclidD(A, B):
    return pow((A[0] - B[0])**2 + (A[1] - B[1])**2, 0.5)

plt.figure(figsize=(12.8, 9.6))
plt.quiver(lon[::8], lat[::8], uf[::8, ::8], vf[::8, ::8], scale = 3, scale_units = 'x')

for k in range(2):
    prevposition = findStart() # lat & lon
    routeX = [prevposition[1]] # lon
    routeY = [prevposition[0]] # lat
    print(f"Start Point : {prevposition}")
    for i in range(1000) :
        position = findNext(prevposition[0], prevposition[1])
        if not position or euclidD(position, prevposition) <= 1e-5:
            break
        routeY.append(position[0]) # lat
        routeX.append(position[1]) # lon
        prevposition = position
    plt.scatter(routeX, routeY, 0.1)

plt.show()
print("E")

