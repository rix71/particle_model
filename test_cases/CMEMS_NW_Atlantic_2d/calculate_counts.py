#!/usr/bin/env python3
# ---------------------
from glob import glob
import ast
import zipfile
import shutil
import tempfile
import warnings
from netCDF4 import Dataset, num2date, date2num
import numpy as np
import datetime
import os
from multiprocessing import Pool
import matplotlib as mpl
warnings.filterwarnings("ignore")


class NcReader:
    SIZE_LIMIT = 30.0  # GB

    def __init__(self, fname, zipname=None):
        self.fname = fname
        self.zipname = zipname
        self.tmpdir = None
        if (self.zipname is not None):
            self.zf = zipfile.ZipFile(self.zipname)
            if (self.get_file_size() > self.SIZE_LIMIT):
                print(
                    f"File {self.fname} is too large. Extracting to temporary directory")
                self.tmpdir = os.getcwd()+"/"+tempfile.mkdtemp()
                self.fname = self.zf.extract(self.fname, path=self.tmpdir)
                self.f = Dataset(self.fname)
            else:
                self.f = Dataset(
                    "dummy", memory=self.zf.open(self.fname).read())
        else:
            self.f = Dataset(self.fname)

    def get_file_size(self):
        # Size of the zipped file in GB
        return self.zf.getinfo(self.fname).file_size/(1024*1024*1024)

    # Destructor
    def __del__(self):
        print(f"Closing {self.fname}")
        self.f.close()
        if (self.zipname is not None):
            print(f"Closing {self.zipname}")
            self.zf.close()
            if (self.tmpdir is not None):
                print(f"Removing {self.tmpdir}")
                shutil.rmtree(self.tmpdir)

    def var_exists(self, vname):
        return vname in self.f.variables

    def get_var(self, vname, indices=None):
        if (indices is not None):
            return np.array(self.f.variables[vname][indices])
        return np.array(self.f.variables[vname][:])

    def get_time(self, vname='time'):
        timein = self.f.variables[vname][:]
        timeunits = self.f.variables[vname].units
        return np.array([datetime.datetime(val.year, val.month, val.day, val.hour, val.minute, val.second) for val in num2date(timein, timeunits)])


def get_counts_by_state(x, y, state, binx, biny):
    x = x.flatten()
    y = y.flatten()
    state = state.flatten()

    u_states = np.unique(state[state > 0]).astype(int)

    nx = len(binx)
    ny = len(biny)
    counts = np.zeros((N_STATES, ny-1, nx-1), dtype=int)
    positions = np.concatenate((x[:, np.newaxis], y[:, np.newaxis]), axis=1)
    nanmask = np.isnan(positions)
    nanmask = nanmask[:, 0] + nanmask[:, 1]
    positions = positions[~nanmask]
    state = state[~nanmask]
    bins = np.array([binx, biny])

    for istate in state_vals.keys():
        if (istate not in u_states):
            continue
        kkstate = np.squeeze(np.where(state == istate))
        if(kkstate.size <= 1):
            kkstate = np.array([kkstate])
        c, _ = np.histogramdd(positions[kkstate], bins)
        counts[istate-1, :, :] = c.T
    return counts


def get_scaled_concentration(x, y, state, ids, binx, biny):
    x = x.flatten()
    y = y.flatten()
    state = state.flatten()
    ids = ids.flatten()

    u_states = np.unique(state[state > 0]).astype(int)
    u_ids = np.unique(ids).astype(int)

    nx = len(binx)
    ny = len(biny)
    counts = np.zeros((N_STATES, ny-1, nx-1), dtype=float)
    positions = np.concatenate((x[:, np.newaxis], y[:, np.newaxis]), axis=1)
    nanmask = np.isnan(positions)
    nanmask = nanmask[:, 0] + nanmask[:, 1]
    positions = positions[~nanmask]
    state = state[~nanmask]
    ids = ids[~nanmask]

    bins = np.array([binx, biny])

    # For each id, get the number of particles in each state
    for istate in state_vals.keys():
        if (istate not in u_states):
            continue
        mstate = state == istate
        for iid in scale_dict.keys():
            if (iid not in u_ids):
                continue
            mid = ids == iid
            c, _ = np.histogramdd(positions[mstate*mid], bins)
            counts[istate-1, :, :] += c.T*scale_dict[iid]
    return counts/cell_area


def save_counts(fname, lon, lat, time_dict, loc_dict=None, **kwargs):
    ncout = Dataset(fname, 'w', "NETCDF4")
    ncout.createDimension('time', None)
    ncout.createDimension('lon', len(lon))
    ncout.createDimension('lat', len(lat))
    timeout = ncout.createVariable('time', "float64", ("time"))
    timeout.units = time_dict["time_units"]
    timeout[:] = date2num(time_dict["time"], time_dict["time_units"])
    lonout = ncout.createVariable('lon', "float64", ("lon"))
    lonout[:] = lon
    latout = ncout.createVariable('lat', "float64", ("lat"))
    latout[:] = lat

    if (loc_dict is not None):
        attr_name = loc_dict.pop("name")
        ncout.setncattr(attr_name, str(loc_dict))

    for arg in kwargs:
        varout = ncout.createVariable(arg, "float64", ("time", "lat", "lon"))
        varout[0, :, :] = kwargs[arg]

    ncout.close()


def merge_files(fnames, outfname, overwrite=False):
    import xarray as xr
    print(f"{len(fnames)} files to merge")
    if ((os.path.exists(outfname)) and (not overwrite)):
        print("Not overwriting. Use -O or --overwrite to overwrite")
        return
    df = xr.open_mfdataset(fnames, combine='by_coords')
    df.to_netcdf(outfname)
    print(f"merged files to {outfname}")


def file_loop(itime):
    fnameout = f"{tmpdir}/counts_{itime:05d}.nc"
    if (resolution is not None):
        fnameout = fnameout.replace(".nc", f".res{int(resolution)}.nc")
    if (os.path.exists(fnameout)):
        return
    xi = x[itime, :]
    yi = y[itime, :]
    statei = state[itime, :]
    timei = time[itime]

    # ! For information only !
    # u_states = sorted(list(set(statei[statei>0])))
    # nstates = len(u_states)
    # print(f"{itime}: {u_states} ({nstates} states)")

    # Get histogram
    counts = get_counts_by_state(
        xi, yi, statei, bin_lon, bin_lat).astype(float)

    if (scale_dict is not None):
        concentration = get_scaled_concentration(
            xi, yi, statei, ids, bin_lon, bin_lat).astype(float)
    else:
        concentration = counts/cell_area

    if(resolution is None):
        counts[:, np.isnan(topo[:-1, :-1])] = np.nan
        concentration[:, np.isnan(topo[:-1, :-1])] = np.nan

    # Write to netcdf
    out_dict = {f"counts_{state_vals[i]}": counts[i-1, :, :]
                for i in state_vals.keys()}
    out_dict["counts_sum"] = np.sum(counts, axis=0)
    out_dict = dict(
        out_dict, **{f"concentration_{state_vals[i]}": concentration[i-1, :, :] for i in state_vals.keys()})

    time_dict = {"time": timei,
                 "time_units": "seconds since 1970-01-01 00:00:00"}

    loc_dict = {f"counts{i}": s for i, s in state_vals.items()}
    loc_dict["name"] = "state_vals"
    save_counts(fnameout, bin_lon[:-1],
                bin_lat[:-1], time_dict, loc_dict=loc_dict, **out_dict)
    print(f"Saved {fnameout}")


state_vals = {
    1: "BEACHED",
    2: "ON_BOUNDARY",
    3: "ACTIVE",
    4: "BOTTOM"
}
N_STATES = 4


def main():
    global tmpdir, resolution, topo, cell_area, bin_lon, bin_lat, x, y, state, ids, time, scale_dict
    import argparse

    TOPOFILE = f"{os.getcwd()}/topo.nc"
    COUNTS_OUT = f"{os.getcwd()}/counts_alltime.nc"
    SCALE_FILE = f"{os.getcwd()}/scale.dat"

    parser = argparse.ArgumentParser(description="Plot particle distribution")
    parser.add_argument(
        "-O", "--overwrite", help=f"Overwrite files", action='store_true')
    parser.add_argument(
        "-s", "--source", help="Path to data file. If <source> is contained in .zip folder, the path is inside to .zip", type=str, required=True)
    parser.add_argument(
        "--zip-file", help="Path to zipped data file if <source> is contained in .zip folder", type=str, required=False)
    parser.add_argument(
        "-o", "--output", help="Output file", type=str, required=False, default=COUNTS_OUT)
    parser.add_argument(
        "-dx", "--resolution", help="Resolution of grid (meters)", type=float)
    parser.add_argument(
        "--ini-file", help="Initial particle position file", type=str, required=False)
    parser.add_argument(
        "--topo-file", help="Topography file", type=str, required=False, default=TOPOFILE)
    parser.add_argument(
        "-np", "--nproc", help="Number of processors", type=int, default=1)
    parser.add_argument("--id-scale", help="Dictionary of scale factors for each particle id",
                        type=str, required=False)
    
    

    args = parser.parse_args()
    overwrite = args.overwrite
    filename = args.source
    resolution = args.resolution
    ini_file = args.ini_file
    final_out = args.output
    num_proc = args.nproc
    topo_file = args.topo_file
    zip_file = args.zip_file
    scale_file = args.id_scale

    print(f"Overwrite: {overwrite}")
    print(f"File name: {filename}")
    print(f"Output file: {final_out}")

    # First check if the file exists
    if (overwrite is False and os.path.exists(final_out)):
        print(f"File {final_out} exists. Use -O to overwrite")
        return

    print(f"Topo file: {topo_file}")
    print(f"Scales file: {scale_file}")
    print(f"Number of processors: {num_proc}")

    # Read topo file
    topo_reader = NcReader(topo_file)

    # Read scales file
    if (scale_file is not None):
        with open(scale_file, 'r') as f:
            scale_dict = ast.literal_eval(f.read())
    else:
        scale_dict = None

    if (zip_file is not None):
        print(f"Zip file: {zip_file}")
        data_reader = NcReader(filename, zip_file)
    else:
        data_reader = NcReader(filename)
    if (resolution is not None):
        print(f"Resolution: {resolution} m")
    if (ini_file is not None):
        print(f"ini file: {ini_file}")

    # Check if the 'state' variable exists before proceeding
    if (not data_reader.var_exists('state')):
        print(f"State variable not found in {filename}")
        return

    # read topo
    lon = topo_reader.get_var('lon')
    lat = topo_reader.get_var('lat')
    topo = topo_reader.get_var('bathymetry')
    topo[topo < -5.] = np.nan

    m2deg = 1./(1852.*60.)
    deg2m = 1852.*60.

    if (resolution is not None):
        dlon_bin = resolution*m2deg/np.cos(np.nanmean(lat) * np.pi/180.)
        bin_lon = np.arange(
            np.nanmin(lon), np.nanmax(lon) + dlon_bin, dlon_bin)
        dlat_bin = resolution*m2deg
        bin_lat = np.arange(
            np.nanmin(lat), np.nanmax(lat) + dlat_bin, dlat_bin)
    else:
        bin_lon = lon
        bin_lat = lat

    dlon = np.diff(bin_lon)
    dlat = np.diff(bin_lat)
    dLon, dLat = np.meshgrid(dlon, dlat)
    dx = dLon * deg2m * np.cos(np.nanmean(lat) * np.pi/180.)
    dy = dLat * deg2m
    # Lucky for us, it has to be one smaller... so no stacking :)
    cell_area = dx * dy

    # Read particle data
    if (scale_file is not None):
        # Don't need id if we have no scale file
        print("Reading 'id' variable...")
        ids = data_reader.get_var("id")
        ids[ids > 9.e20] = np.nan

    print("Reading 'state' variable...")
    state = data_reader.get_var("state")
    state[state < 0] = -1

    # ! For info only !
    u_states = np.unique(state[state > 0]).astype(int)
    nstates = len(u_states)
    print(f"States ({nstates}): {u_states}")

    print("Reading 'lon' and 'lat' variables...")
    x = data_reader.get_var('lon')
    x[x > 999.] = np.nan
    y = data_reader.get_var('lat')
    y[y > 999.] = np.nan
    print("Reading 'time' variable...")
    time = data_reader.get_time()

    ntimes, nparticles = x.shape

    # Make temporary directory
    # tmpdir = f"{os.getcwd()}/tmp"
    tmpdir = os.getcwd() + tempfile.mkdtemp()
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    # Loop over time
    # par_args = [(i, x[i], y[i], state[i], time[i]) for i in range(ntimes)]
    if (num_proc > 1):
        with Pool(num_proc) as pool:
            # res = pool.starmap(file_loop, par_args)
            res = pool.map(file_loop, range(ntimes))
    else:
        for i in range(ntimes):
            file_loop(i)

    # Merge files
    print("Merging files...")
    files_to_merge = glob(
        f"{tmpdir}/*.nc") if resolution is None else glob(f"{tmpdir}/*.res{int(resolution)}.nc")
    merge_files(files_to_merge,
                outfname=final_out, overwrite=overwrite)

    # Delete temporary directory
    shutil.rmtree(tmpdir)


if __name__ == "__main__":
    from time import time as tictoc
    print("="*80)
    print("Calculating particle counts")
    print("="*80)
    ts = tictoc()
    main()
    te = tictoc()
    print("="*80)
    print(f"Done [{te - ts} seconds]")
    print("="*80)
