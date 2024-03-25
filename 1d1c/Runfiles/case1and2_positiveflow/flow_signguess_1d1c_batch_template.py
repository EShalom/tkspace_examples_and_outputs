#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 13:07:27 2021

@author: pyess
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from time import gmtime, strftime

from tkspace.utilities import mkdir_p as mkdir
from tkspace.forward_models.onedim.sampling import downsampler
from tkspace.forward_models.onedim.sampling import gaussian_noise as noise
from tkspace.forward_models.onedim.systeminfo import geometry
from tkspace.forward_models.onedim.systeminfo import multires_1d as multires
from tkspace.forward_models.onedim.loadcase import load_P

from tkspace.forward_models.onedim.flow import onedim_onecomp_flow_diff as fmod
from tkspace.forward_models.onedim.fm_utilities import fmodsetup_1d1c as fmod_setup
from tkspace.forward_models.onedim.fm_utilities import interp_linear as interp
from tkspace.forward_models.onedim.cases import case1_1D1C as gen_case
from tkspace.inversion_methods.guess import guess_1d1c as gen_guess
from tkspace.forward_models.onedim.loadcase import load_flow_1d1c as load_kinetics
from tkspace.inversion_methods.gradient_descent import gradient_descent_1d1c_flow_setup as inverse_method
from tkspace.forward_models.onedim.flow import onedim_onecomp_flow_diff as fitmod
from tkspace.inversion_methods.gradient_descent import simple_gradient_descent as method
from numpy import heaviside as H

from tkspace.inversion_methods.cost_functions import rms_psp_cost as cost

dims        = 1 # what dimension to be used
pic         = 1 # 1 = normal local conc, 2= tissue conc
comp        = 1 # true compartment number
fitcomp     = 1   # number of compartments considered in the fit

guess_no    = PH1 # guess type (1 - all optimsable params is 1 stack, 2 - multiple Pfix columns for multi stage fits)
case        = PH2 # case number for system type (detailed in testsystems.py)
SNR         = PH3 # SNR for run
Dt          = PH4  # sampling time (s)
maxiter     = PH5
updatemethod = 'adam'
step        = PH6
fdir = 'PH7'#'server_runs_17102022/'+'case{}_guess{}_{}snr_{}step_{}Dt_{}maxiter'.format(case,guess_no,SNR,step,Dt,maxiter)
mkdir(fdir)

int_res, meas_res, sysdim, voxels, meas_voxels, max_vals = geometry()

Lx, Ly, Lz, T = sysdim
Nx, Ny, Nz = meas_voxels
dx, dy, dz = meas_res
dxyz = dx * dy * dz 

f, D, v, Jp, Jn, dt, Nt, t = load_kinetics(dims,comp,case)


C_tiss = fmod(Lx, dx, t, dt, f, D,v, Jp,Jn)
C_ds, tsamp = downsampler(C_tiss, T, dt, Dt, 0, 0)

Nt_samp = len(tsamp)
print("Parameters allocated \n")

# Signal to noise ratio ( 0 = noiseless) calculated as: noise stdev = mean(signal)/SNR

C_meas = noise(C_ds, SNR) # add guassian noise at scale of CNR

print("Downsample and noise addition complete \n")

# downsample known parameters to directly compare with inversion outputs
# allocate multi res options if using

Nx, dx, xc_fine, x_fine = multires(dx, Lx)

print("Forwards model complete \n")

print("Starting Optimisation.... \n")

# create output folder with open files for cost function results and Param lists
Pnorm, Pmax, dt, Nt, t = load_P(dims, comp, case)
 
args = [C_meas, fitmod, fmod_setup, meas_res, sysdim, tsamp, comp, pic, step, dx]
#%%
vguess = C_meas[:,-1]/max(C_meas[:,-1])

f1_guess= [11,7,4,9,5,2]#,10,8,6,4,-11,-7,-4]

f1_init = f1_guess[PH1-1]

D_init = np.zeros(Nx+1)
v_init = vguess

skip = int(Dt/dt)
start_time = 12
Jp_init = C_meas[0,:start_time-1]*H(f1_init,0)#np.zeros(start_time-1)
Jn_init = C_meas[-1,:start_time-1]*H(-f1_init,0)#np.zeros(start_time)#C_meas[-1,:start_time]*abs(ua_init[-1])*dxyz*1e-3#

init_og = [f1_init, D_init,v_init, Jp_init, Jn_init]
init = [f1_init, D_init,v_init, Jp_init, Jn_init]
fix_vals = [1,0,1,1] #(f,D,v, J)
lb = [-1,0,0,0] #(f,D,v, J)
ub = [1,1,1,1] #(f,D,v, J)
exit_change = [1e-4,1e-9,1e-4,0] # still working on how to enforce this
max_vals = [14,4e-6,1,15]#[3.5,4e-6,0.8 * (Lx*Ly*Lz) * 1e-3] # umax 3.5(cm/s), Kva_max (1/s),Dmax (cm^2/s), Jmax (mmol/s)

np.savez_compressed('{}/args_in_opt'.format(fdir), init = init, sysdim=sysdim,x_fine=x_fine,xc_fine=xc_fine,C_meas=C_meas, tsamp=tsamp, C_tiss=C_tiss,t=t, fitmod=fitmod, Pnorm=Pnorm,Pmax=Pmax,int_res=int_res,meas_res=meas_res, dims=str(dims),comp=comp,case=case)

for i in range(start_time, Nt_samp+1):
    args = [C_meas[:,0:i], fitmod, fmod_setup, meas_res, sysdim, tsamp[0:i], comp, pic, step, dx]

    value_tuple, Precon, fun, Pmax = inverse_method(method,updatemethod,init, fix_vals, lb, ub, exit_change, cost, maxiter, max_vals, fdir, args)
    f1_inv, D_inv,v_inv,Jp_inv, Jn_inv = value_tuple
    
    if i != Nt_samp:
        Jp_init = np.append(Jp_inv,Jp_inv[-1])#Jpa_inv[-1])#np.append(Jpa_inv,Jpa[i*skip])#
        Jn_init =np.append(Jn_inv,Jn_inv[-1])#Jna_inv[-1])
    init = [f1_inv,  D_inv,v_inv, Jp_init, Jn_init]
    print("Finshed fitting to {}s".format(tsamp[i-1]))
    print (strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    
    np.savez_compressed('{}/FinalOptimisationValues_{}s'.format(fdir,tsamp[i-1]), fun=fun, Precon=Precon,Pmax=Pmax, value_tuple=value_tuple, i=i)
f1_inv,  D_inv, v_inv, Jp_inv, Jn_inv = value_tuple
print("finished! \n")
