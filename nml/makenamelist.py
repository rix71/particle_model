#!/usr/bin/env python3
#---------------------
import os, sys, ast
import f90nml
#---------------------
def getVarType(varin):
  try:
    convertedVar = ast.literal_eval(varin)
    return convertedVar
  except:
    # ast doesn't detect strings?!
    return varin
  
def main(args):
  templatelocation = f"{os.getcwd()}/nmltemplate"
  if len(args) > 1:
      templatelocation = args[1]
  nmltempfile = templatelocation
  nmltemplate = f90nml.read(nmltempfile)

  for namelist in nmltemplate.keys():
    for varname in nmltemplate[namelist].keys():
      varout = os.getenv(varname)
      if varout is None: continue # Default from template
      else:
        nmltemplate[namelist][varname] = getVarType(varout)

  with open('input.inp', 'w') as nmlout:
    f90nml.patch(nmltempfile, nmltemplate, nmlout)
  
  return 0
#---------------------
if __name__ == "__main__":
  sys.exit(main(sys.argv))
