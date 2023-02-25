import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

print("Plotting trajectory...")

df = xr.open_dataset("./out/eddy.out.nc")
currents = xr.open_dataset("./data/data/eddy.nc")
topo = xr.open_dataset("./data/eddy_topo.nc")

fig = plt.figure()
ax = fig.add_subplot(111)
ax.quiver(topo.lon[::10], topo.lat[::10],
          currents.u[0, 0, ::10, ::10], currents.v[0, 0, ::10, ::10])
# Plot every 5th particle
ax.plot(df.isel(particle=slice(0, None, 5)).lon, df.isel(
    particle=slice(0, None, 5)).lat, 'go-', linewidth=0.5, markersize=0.5)
ax.plot(df.isel(particle=slice(0, None, 5), time=0).lon,
        df.isel(particle=slice(0, None, 5), time=0).lat, 'r.')

ax.set_aspect('equal')
fig.savefig("eddy.png", dpi=300, bbox_inches='tight')

print("Saved to eddy.png")
