import os 
import json
import h5py
import numpy as np
import sympy as syms
from scipy.sparse import diags
from sympy import lambdify
from numpy import linalg as LA
import matplotlib.pyplot as plt
path = os.path.dirname(os.path.realpath(__file__))

# Recursive function that traverses the json file and 
# converts it to hdf5
def json_to_hdf5(key,val,group):
  if isinstance(val,(dict)):
    grp = group.create_group(str(key))
    data = val
    for key,val in data.items():
      json_to_hdf5(key,val,grp)
  elif isinstance(val,str):
    dset = group.create_dataset(str(key),(1,),dtype = "S"+str(len(val)))
    dset[0] =  np.string_(val) 
  elif isinstance(val,int):
    dset = group.create_dataset(str(key),(1,),dtype = "i4")
    dset[0] =  val 
  elif isinstance(val,float):
    dset = group.create_dataset(str(key),(1,),dtype = "f8")
    dset[0] =  val
  elif isinstance(val,list):
    if all(isinstance(x,int) for x in val):
      dset = group.create_dataset(str(key),(len(val),),dtype = "i4")
      dset[:] =  val[:] 
    elif all(isinstance(x,float) for x in val):
      dset = group.create_dataset(str(key),(len(val),),dtype = "f8")
      dset[:] =  val[:]      
    elif isinstance(val[0],dict):
      for i in range(len(val)):
        grp = group.create_group(str(key)+str(i))
        for k,v in val[i].items():
          json_to_hdf5(k,v,grp)
    else:
      print("Innaproprate datatype present in dict")
  else:
    print("Innaproprate datatype present in input.json")

# Gets the current working directory
cwd = os.getcwd()
# If an old version of the parameters.h5 file exits remove it
# since we are are reading in a new input.json
os.system('rm parameters.h5')
# Create a new parameters.h5 file to write to
params = h5py.File('parameters.h5','w')

# Read the json file and write it to hdf5.
with open(cwd+"/input.json", 'r') as f:
  data = json.load(f)
  for key,val in data.items():
    json_to_hdf5(key,val,params)




dset = params.create_dataset("install_directory",(1,),dtype = "S"+str(len(path)))
dset[0] = np.string_(path)
params.close()


# Compute basis
if data["EPS"]["compute"] == 1 :
  if data["mpi"]["assign_host"] == 1:
    print("basisf90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " --host " + str(data["mpi"]["host"]) + " "+path+"/basisf90	-eps_view_values :eigenvalues.m:ascii_matlab -eps_view_vectors :myeigenvectors.m:ascii_matlab -eps_monitor_conv")
  else:
    print("basisf90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " "+path+"/basisf90 -eps_view_values :myeigenvalues.m:ascii_matlab	-eps_view_vectors :eigenvectors.m:ascii_matlab -eps_monitor_conv")

