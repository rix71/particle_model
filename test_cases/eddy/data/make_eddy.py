import xarray as xr
import numpy as np

print("Generating eddy data...")

# Make a fake dataset
x = np.linspace(-10, 10, 600)
y = np.linspace(-10, 10, 500)

xx, yy = np.meshgrid(x, y)

u = np.ones(shape=(100, 1, 500, 600))*-np.cos(np.pi/180.)*yy*10.
v = np.ones(shape=(100, 1, 500, 600))*np.cos(np.pi/180.)*xx*10.

# Save to netCDF
ds = xr.Dataset(
    data_vars={
        'u': (('time', 'level', 'lat', 'lon'), u),
        'v': (('time', 'level', 'lat', 'lon'), v),
    },
    coords={
        'time': np.arange(100)*600,
        'level': np.arange(1),
        'lat': y,
        'lon': x,
    },
)
ds.time.attrs['units'] = 'seconds since 2018-01-01 00:00:00'
encoding = {
    'time': {'dtype': 'float64'},
    'lon': {'dtype': 'float64'},
    'lat': {'dtype': 'float64'},
    'level': {'dtype': 'float64'},
}

ds.to_netcdf('./data/eddy.nc', encoding=encoding)
print("Saved to ./data/eddy.nc")


# Save topo
print("Generating eddy topo...")

bathy = np.ones_like(xx)*10.
ds = xr.Dataset(
    data_vars={
        'bathymetry': (('lat', 'lon'), bathy),
    },
    coords={
        'lat': y,
        'lon': x,
    },
)
ds.to_netcdf('./eddy_topo.nc')
print("Saved to ./eddy_topo.nc")
