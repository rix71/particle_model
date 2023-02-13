import numpy as np
import xarray as xr


print("Creating initial particle positions...")

topo = xr.open_dataset('eddy_topo.nc')


nx = topo.lon.size
ini_x = topo.lon.values[nx//2-150:nx//2+150:40]
ny = topo.lat.size
ini_y = topo.lat.values[ny//2-150:ny//2+150:40]

ini_x, ini_y = np.meshgrid(ini_x, ini_y)

ini_x = ini_x.flatten()
ini_y = ini_y.flatten()

nparticles = ini_x.size

print(f"Number of particles: {nparticles}")

# Save to file
with open("ini_particles.dat", "w") as f:
    f.write(f"{nparticles}\n")
    for i in range(ini_x.size):
        f.write(f"{ini_x[i]} {ini_y[i]} 0.0 1.0 9999999.0 1000.0 0.0005\n")

print("Saved to ini_particles.dat")
