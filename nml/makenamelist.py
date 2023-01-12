#!/usr/bin/env python3
# ---------------------
import os
import sys
import ast
import f90nml
import argparse
# ---------------------


def getVarType(varin):
    try:
        convertedVar = ast.literal_eval(varin)
        return convertedVar
    except:
        # ast doesn't detect strings?!
        return varin


def main():
    parser = argparse.ArgumentParser(
        description='Generate namelist from environment variables')
    parser.add_argument(
        "-i", "--input", help="Namelist template file", default="./nmltemplate")
    parser.add_argument(
        "-o", "--output", help="Namelist output file", default="./input.inp")
    args = parser.parse_args()

    nmltempfile = args.input
    outfile = args.output
    
    nmltemplate = f90nml.read(nmltempfile)

    for namelist in nmltemplate.keys():
        for varname in nmltemplate[namelist].keys():
            varout = os.getenv(varname)
            if varout is None:
                continue  # Default from template
            else:
                nmltemplate[namelist][varname] = getVarType(varout)

    with open(outfile, 'w') as nmlout:
        f90nml.patch(nmltempfile, nmltemplate, nmlout)

    return 0


# ---------------------
if __name__ == "__main__":
    sys.exit(main())
