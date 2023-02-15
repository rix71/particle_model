import xarray as xr
import matplotlib.pyplot as plt
import argparse
import warnings
warnings.filterwarnings("ignore")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--source", help="Source files", required=True)
    parser.add_argument("-o", "--output-path",
                        help="Output path", required=False, default=".")

    args = parser.parse_args()

    print("="*80)
    print("Plotting counts")
    print("="*80)

    df = xr.open_mfdataset(args.source)

    df.counts_BEACHED.attrs['units'] = '# of particles in cell'
    df.counts_ON_BOUNDARY.attrs['units'] = '# of particles in cell'
    df.counts_ACTIVE.attrs['units'] = '# of particles in cell'
    df.counts_BOTTOM.attrs['units'] = '# of particles in cell'
    df.counts_sum.attrs['units'] = '# of particles in cell'

    df.concentration_BEACHED.attrs['units'] = 'particles/m2'
    df.concentration_ON_BOUNDARY.attrs['units'] = 'particles/m2'
    df.concentration_ACTIVE.attrs['units'] = 'particles/m2'
    df.concentration_BOTTOM.attrs['units'] = 'particles/m2'

    for varname in df.data_vars:
        if (varname.startswith("concentration_")):
            continue
        print(varname)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        vmax = df[varname].max()
        vmax = 10
        df[varname].mean(dim="time").plot(ax=ax, cmap=plt.cm.rainbow, vmin=0., vmax=vmax, cbar_kwargs={
            "pad": 0.01, "label": f"{varname} [{df[varname].attrs['units']}]"})
        fig.savefig(f"{args.output_path}/fig_{varname}.png")

    print("="*80)
    print("Done")
    print("="*80)


if __name__ == "__main__":
    main()
