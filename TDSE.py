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
light = 137.035999084
def potential(c0,zc,a,b,c,r,cap_present,rot_present,gobbler,cap_eta=None,cap_re_alpha1=None,cap_im_alpha1=None,cap_re_alpha2=None,cap_im_alpha2=None,cap_beta=None,rot_theta=None):
  if rot_present == 1:
    r = np.exp(1j*rot_theta)*r
  V =  np.zeros((2,len(r)))
  V[0,:] += np.real(-c0/r)
  V[0,:] += np.real(-zc*np.exp(-c*r)/r)
  V[1,:] += np.imag(-c0/r)
  V[1,:] += np.imag(-zc*np.exp(-c*r)/r)
  for i in range(len(a)):
    V[0,:] += np.real(-a[i]*np.exp(-b[i]*r))
    V[1,:] += np.imag(-a[i]*np.exp(-b[i]*r))
  if cap_present == 1:
    iabs = max(int(round(gobbler*Rmax/dr)) - 1,0)
    if cap_eta != None:
      V[1,iabs:] += -cap_eta*np.sin(np.pi*(r[iabs:]/Rmax-gobbler)/(2.0*(1.0-gobbler)))**2.0      
    elif cap_re_alpha1 != None and cap_im_alpha1 != None and cap_re_alpha1 != None and cap_im_alpha2 != None and cap_beta != None:
      V[0,:] += -0.5*((cap_re_alpha1/cap_beta**2)/(np.abs(np.exp(-(r[:]-r[-1])/(2*cap_beta))))-(cap_im_alpha2/cap_beta**2)/(np.abs(np.exp(-(r[:]-r[-1])/(2*cap_beta))))**2)
      V[1,:] += -0.5*((cap_im_alpha1/cap_beta**2)/(np.exp(-(r[:]-r[-1])/(2*cap_beta)))+(cap_re_alpha2/cap_beta**2)/(np.exp(-(r[:]-r[-1])/(2*cap_beta)))**2)
  return V

def electric_field(lam,I,shape,frequency_shift,t,num_cycles,t_start,custom_envalope_phase,cep,ellipticity,normal_vec,perp_vec):
  normal_vec = normal_vec/LA.norm(normal_vec,2)
  E =  np.zeros((len(t),3))
  conEV = 27.2110
  w = 1239.8006/(lam*conEV)
  conWpCM2 = 3.5101e16
  E0 = np.sqrt(I/conWpCM2)
  if shape == 'sin2': 
    if frequency_shift == 1:
      mu = 4.0*np.arcsin(np.exp(-0.25))**2.0
    else:
      mu = 0.0
    wa = 2.0*w/(1.0 + np.sqrt(1.0 + mu/num_cycles**2.0))
    T0 = num_cycles*2.0*np.pi/wa
    t_cep = T0/2.0 + t_start
    if custom_envalope_phase == 0:
      t_cep = T0/2.0 + t_start
      cep  = 0.0
    t_end = t_start + T0

    if ellipticity == 0.0:
      A, tt = syms.symbols('A tt')
      A = light*E0*syms.sin(np.pi*(tt-t_start)/T0)**2*syms.sin(wa*(tt-t_cep)+cep)/wa
      F = lambdify(tt,-syms.diff(A,tt)/light, 'numpy')
      for i in range(len(t)):
        if t_start <= t[i] < t_end:
          AMP = F(t[i])
          E[i,0] = normal_vec[0]*AMP
          E[i,1] = normal_vec[1]*AMP
          E[i,2] = normal_vec[2]*AMP
    else:
      A0, A1, tt = syms.symbols('A0 A1 tt')
      perp_vec_1 = np.cross(normal_vec,perp_vec)
      A0 = (light*E0/wa)*syms.sin(np.pi*(tt-t_start)/T0)**2*syms.sin(wa*(tt-t_cep)+cep)/np.sqrt(1+ellipticity**2)
      A1 = (light*E0/wa)*ellipticity*syms.sin(np.pi*(tt-t_start)/T0)**2*syms.cos(wa*(tt-t_cep)+cep)/np.sqrt(1+ellipticity**2)
      F0 = lambdify(tt,-syms.diff(A0,tt)/light, 'numpy')
      F1 = lambdify(tt,-syms.diff(A1,tt)/light, 'numpy')
      for i in range(len(t)):
        E[i,0] = perp_vec[0]*F0(t[i])+perp_vec_1[0]*F1(t[i])
        E[i,1] = perp_vec[1]*F0(t[i])+perp_vec_1[1]*F1(t[i])
        E[i,2] = perp_vec[2]*F0(t[i])+perp_vec_1[2]*F1(t[i])

  elif shape == 'gaussian': 
    if frequency_shift == 1:
      mu = 8.0*np.log(2.0)/np.pi**2.0
    else:
      mu = 0.0
    wa = 2.0*w/(1.0 + np.sqrt(1.0 + mu/num_cycles**2.0))
    T0 = num_cycles*2.0*np.pi/wa
    t_cep = T0/2.0 + t_start
    if custom_envalope_phase == 0:
      t_cep = 5.0*T0 + t_start
      cep  = 0.0
    t_end = t_start + 10.0*T0
          
    if ellipticity == 0.0:
      A, tt = syms.symbols('A tt')
      A = light*E0*2.0**(-(4.0*(-tt + t_cep)**2.0)/T0**2.0)*syms.sin(wa*(tt-t_cep)+cep)/wa
      F = lambdify(tt,-syms.diff(A,tt)/light, 'numpy')
      for i in range(len(t)):
        if t_start <= t[i] < t_end:
          AMP = F(t[i])
          E[i,0] = normal_vec[0]*AMP
          E[i,1] = normal_vec[1]*AMP
          E[i,2] = normal_vec[2]*AMP
    else:
      A0, A1, tt = syms.symbols('A0 A1 tt')
      perp_vec_1 = np.cross(normal_vec,perp_vec)
      A0 = (light*E0/wa)*2.0**(-(4.0*(-tt + t_cep)**2.0)/T0**2.0)*syms.sin(wa*(tt-t_cep)+cep)/np.sqrt(1+ellipticity**2)
      A1 = (light*E0/wa)*ellipticity*2.0**(-(4.0*(-tt + t_cep)**2.0)/T0**2.0)*syms.cos(wa*(tt-t_cep)+cep)/np.sqrt(1+ellipticity**2)
      F0 = lambdify(tt,-syms.diff(A0,tt)/light, 'numpy')
      F1 = lambdify(tt,-syms.diff(A1,tt)/light, 'numpy')
      for i in range(len(t)):
        E[i,0] = perp_vec[0]*F0(t[i])+perp_vec_1[0]*F1(t[i])
        E[i,1] = perp_vec[1]*F0(t[i])+perp_vec_1[1]*F1(t[i])
        E[i,2] = perp_vec[2]*F0(t[i])+perp_vec_1[2]*F1(t[i])
  return E

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
Rmax = data["EPS"]["R_max"]
dr = data["EPS"]["delta_x"]
num_points = int(round(Rmax/dr))
r = np.zeros((num_points,))
for i in range(num_points):
  r[i] = dr*(i+1)
a = data["EPS"]["nuclei"]["a"]
b = data["EPS"]["nuclei"]["b"]
c = data["EPS"]["nuclei"]["c"]
c0 = data["EPS"]["nuclei"]["c0"]
zc = data["EPS"]["nuclei"]["Zc"]
cap_present = data["EPS"]["cap_present"]
rot_present = data["EPS"]["rot_present"]
rot_theta = data["EPS"]["rot_theta"]
eta = data["EPS"]["cap_eta"]
beta = data["EPS"]["cap_beta"]
re_a1 = data["EPS"]["cap_re_alpha1"]
im_a1 = data["EPS"]["cap_im_alpha1"]
re_a2 = data["EPS"]["cap_re_alpha2"]
im_a2 = data["EPS"]["cap_im_alpha2"]
gobbler = data["EPS"]["gobbler"]

if rot_present == 1:
  V = potential(c0,zc,a,b,c,r,cap_present,rot_present,gobbler,rot_theta=rot_theta)
else:
  if eta > 1e-10: 
    V = potential(c0,zc,a,b,c,r,cap_present,rot_present,gobbler,cap_eta = eta)
  elif re_a1 > 1e-10 or im_a1 > 1e-10 or re_a2 > 1e-10 or im_a2 > 1e-10:
    V = potential(c0,zc,a,b,c,r,cap_present,rot_present,gobbler,cap_re_alpha1 = re_a1,cap_im_alpha1 = im_a1,cap_re_alpha2 = re_a2,cap_im_alpha2 = im_a2,cap_beta=beta)
  else:
    V = potential(c0,zc,a,b,c,r,cap_present,rot_present,gobbler)

plt.figure()
plt.plot(r,V[0,:])
plt.plot(r,V[1,:])
plt.savefig('potential.png')

dV =((diags(np.ones(num_points-1),1)-diags(np.ones(num_points-1),-1))/(2*r[1])).dot(np.transpose(V))

dV = np.transpose(dV)


dset = params["EPS"].create_dataset("V",(2,num_points),dtype = "f8")
dset[:,:] = V[:,:]
dset = params["EPS"].create_dataset("r",(num_points,),dtype = "f8")
dset[:] = r[:]
dset = params["EPS"].create_dataset("dV",(2,num_points),dtype = "f8")
dset[:,:] = dV[:,:]
conEV = 27.2110
dt = np.zeros((len(data["laser"]["pulse"]),))
L = np.zeros((len(data["laser"]["pulse"]),))
t_start = np.zeros((len(data["laser"]["pulse"]),))
t_end = np.zeros((len(data["laser"]["pulse"]),))
for i in range(len(data["laser"]["pulse"])):
  if data["laser"]["pulse"][i]["pulse_shape"] == 'sin2': 
    if data["laser"]["frequency_shift"] == 1:
      mu = 4.0*np.arcsin(np.exp(-0.25))**2.0
    else:
      mu = 0.0
  elif data["laser"]["pulse"][i]["pulse_shape"] == 'gaussian': 
    if data["laser"]["frequency_shift"] == 1:
      mu = 8.0*np.log(2.0)/np.pi**2.0
    else:
      mu = 0.0

  w = 1239.8006/(data["laser"]["pulse"][i]["wavelength"]*conEV)
  wa = 2.0*w/(1.0 + np.sqrt(1.0 + mu/data["laser"]["pulse"][i]["cycles"]**2.0))
  T0 = data["laser"]["pulse"][i]["cycles"]*2.0*np.pi/wa
  dt[i] = data["laser"]["pulse"][i]["relative_delay"]
  if data["laser"]["pulse"][i]["pulse_shape"] == 'sin2':
    L[i] = T0
  elif data["laser"]["pulse"][i]["pulse_shape"] == 'gaussian':
    L[i] = 2.0*data["laser"]["pulse"][i]["gaussian_length"]*T0
  if i == 0:
    Lfirst = L[0]
    Llast  = L[0]
    dtfirst = dt[0]
    dtlast = dt[0]
  for j in range(i+1):
    if L[i]/2.0-dt[i] > L[j]/2.0-dt[j]:
      Lfirst = L[i]
      dtfirst = dt[i]
    if L[i]/2.0+dt[i] > L[j]/2.0+dt[j]:
      Llast = L[i]
      dtlast = dt[i]

max_time = 0.5*(Lfirst+Llast) + dtlast - dtfirst 
for i in range(len(data["laser"]["pulse"])):
  t_start[i] = (Lfirst/2.0 - dtfirst) - (L[i]/2.0-dt[i])
  t_end[i] = t_start[i] + L[i]
num_steps = int(np.ceil(max_time/data["TDSE"]["delta_t"]))+1
t = np.zeros((num_steps,))
dset = params["TDSE"].create_dataset("t",(num_steps,),dtype="f8")
for i in range(num_steps):
  t[i] = i*data["TDSE"]["delta_t"]
dset[:] = t[:]
E = np.zeros((t.shape[0],3))
for i in range(len(data["laser"]["pulse"])):
  lam = data["laser"]["pulse"][i]["wavelength"]
  I = data["laser"]["pulse"][i]["intensity"]
  shape = data["laser"]["pulse"][i]["pulse_shape"]
  frequency_shift = data["laser"]["frequency_shift"]
  num_cycles = data["laser"]["pulse"][i]["cycles"]
  custom_envalope_phase = data["laser"]["pulse"][i]["custom_envalope_phase"]
  cep = data["laser"]["pulse"][i]["cep"]
  ellipticity = data["laser"]["pulse"][i]["ellipticity"]
  normal_vec = data["laser"]["pulse"][i]["normal_vec"]
  perp_vec = data["laser"]["pulse"][i]["perp_vec"]
  Ei = electric_field(lam,I,shape,frequency_shift,t,num_cycles,t_start[i],custom_envalope_phase,cep,ellipticity,normal_vec,perp_vec)
  dset = params["laser"]["pulse"+str(i)].create_dataset("E"+str(i),(num_steps,3),dtype = "f8")
  dset[:,:] = Ei[:,:]
  E[:,:] += Ei[:,:]
dset = params["laser"].create_dataset("E",(num_steps,3),dtype = "f8")
dset[:,:] = E[:,:]
dset = params["laser"].create_dataset("E_length",(1,),dtype = "i4")
dset[0] = num_steps

dset = params.create_dataset("install_directory",(1,),dtype = "S"+str(len(path)))
dset[0] = np.string_(path)
params.close()


# Compute basis
if data["EPS"]["compute"] == 1 :
  if data["mpi"]["assign_host"] == 1:
    print("basisf90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " --host " + str(data["mpi"]["host"]) + " "+path+"/basisf90")
  else:
    print("basisf90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " "+path+"/basisf90")

params = h5py.File('parameters.h5','r+')
if data["EPS"]["local"] == 0:
  atom = h5py.File(str(data["EPS"]["location"]) +"/" + str(data["EPS"]["label"])+".h5","r")
else:
  atom = h5py.File(str(data["EPS"]["label"])+".h5","r")
block_n = np.array([])
block_l = np.array([])
if data["block_state"]["sfa"] == 1:
  for l in range(data["TDSE"]["l_max"]+1):
    for n in range(l+1,data["TDSE"]["n_max"]+1):
      if n == 1 and l == 0:
        continue
      if atom["Energy_l"+str(l)][0][n-l-1] < 0.0 and \
        np.abs(atom["Energy_l"+str(l)][1][n-l-1]) < np.log(1.0/data["block_state"]["sfa_tol"])/max_time:
        block_n = np.append(block_n,n)
        block_l = np.append(block_l,l)
  del params["block_state"]["num_block"] 
  dset = params["block_state"].create_dataset("num_block",(1,),dtype = "i4")
  dset[0] = len(block_n)
  del params["block_state"]["n_index"]
  dset = params["block_state"].create_dataset("n_index",(len(block_n),),dtype = "i4")
  dset[:] = block_n[:]
  del params["block_state"]["l_index"]
  dset = params["block_state"].create_dataset("l_index",(len(block_n),),dtype = "i4")
  dset[:] = block_l[:]

params.close()


# Compute operators  
if data["operators"]["compute"] == 1:
  # If manually assigning host (done when running on multiple nodes)
  if data["mpi"]["assign_host"] == 1:
    print("generateH0f90")
    os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
      " --host " + str(data["mpi"]["host"]) + \
        " "+path+"/generateH0f90")
    if data["TDSE"]["l_max"] > 0:
      print("generateH1f90")
      os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"]))+ \
        " --host " + str(data["mpi"]["host"]) + \
          " "+path+"/generateH1f90")
      print("generateDipoleAccelerationf90")
      os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
        " --host " + str(data["mpi"]["host"]) + \
          " "+path+"/generateDipoleAccelerationf90")
  else:
    print("generateH0f90")
    os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
      " "+path+"/generateH0f90")
    if data["TDSE"]["l_max"] > 0:
      print("generateH1f90")
      os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
        " "+path+"/generateH1f90")
      print("generateDipoleAccelerationf90")
      os.system("mpirun -np " + str(min(data["mpi"]["np"],data["TDSE"]["l_max"])) + \
        " "+path+"/generateDipoleAccelerationf90")

# Propagate 
if data["TDSE"]["propagate"] == 1:
  # If running on multiple nodes 
  if data["mpi"]["assign_host"] == 1:
    print("propagatef90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " --host " + str(data["mpi"]["host"]) + " "+path+"/propagatef90")
  else:
    print("propagatef90")
    os.system("mpirun -np " + str(data["mpi"]["np"]) + \
      " "+path+"/propagatef90")
    #end for 
