import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

# 2D simulation doesn't actually need bathymetry, but it does need the grid and land mask
# So we'll just use the mask from the CMEMS dataset and make a fake bathymetry (land = -10, ocean = 10)

df = xr.open_dataset("./../hydro/cmems_hydro.nc")


mask = np.isnan(df.uo.values[0, 0, :, :])


topo = np.ones_like(mask, dtype=np.float64)
topo[mask] *= -10.
topo[~mask] *= 10.


# Save to netCDF
ds = xr.Dataset(
    data_vars={
        'bathymetry': (('lat', 'lon'), topo),
    },
    coords={
        'lat': df.latitude.values,
        'lon': df.longitude.values,
    },
)
ds.to_netcdf('./topo.nc')

# plt.pcolormesh(topo)
# plt.colorbar()
# plt.show()
