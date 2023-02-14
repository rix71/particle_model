import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------
PARTICLE_DENSITY = 920.0
PARTICLE_RADIUS = 0.0005

# ------------------------------------
BOX_LIMIT = 30
# ------------------------------------
topo = xr.open_dataset("../topo/topo.nc")
topo.bathymetry.values[topo.bathymetry.values < -5.] = np.nan
seamask = np.isnan(
    topo.bathymetry.values[BOX_LIMIT:-BOX_LIMIT, BOX_LIMIT:-BOX_LIMIT])
seamask_f = seamask.flatten()
# ------------------------------------
x, y = np.meshgrid(
    topo.lon.values[BOX_LIMIT:-BOX_LIMIT], topo.lat.values[BOX_LIMIT:-BOX_LIMIT])
x = x.flatten()
y = y.flatten()

# ------------------------------------
lon_ini = x[~seamask_f]
lat_ini = y[~seamask_f]

npart = len(lon_ini)

# with open("init_particles_full.dat", "w") as f:
#     f.write(f"{npart}\n")
#     for i in range(npart):
#         f.write(
#             f"{lon_ini[i]} {lat_ini[i]} 0.0 1.0 259200.0 {PARTICLE_DENSITY} {PARTICLE_RADIUS}\n")


# plt.pcolormesh(topo.lon.values, topo.lat.values, topo.bathymetry,
#                alpha=0.5, cmap='gray', ec='k', lw=0.01)
# plt.scatter(lon_ini, lat_ini, s=1.5, c='r')
# plt.show()
# ------------------------------------
# Reduce number of particles
lon_ini = lon_ini[::50]
lat_ini = lat_ini[::50]

npart = len(lon_ini)

with open("init_particles_50.dat", "w") as f:
    f.write(f"{npart}\n")
    for i in range(npart):
        f.write(
            f"{lon_ini[i]} {lat_ini[i]} 0.0 1.0 259200.0 {PARTICLE_DENSITY} {PARTICLE_RADIUS}\n")


plt.pcolormesh(topo.lon.values, topo.lat.values, topo.bathymetry,
               alpha=0.5, cmap='gray', ec='k', lw=0.01)
# plt.scatter(lon_ini, lat_ini, s=1.5, c='r')
plt.plot(-13.555633544921875, 59.86676788330078, 'ro')
plt.show()
