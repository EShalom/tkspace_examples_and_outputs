#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 10:09:18 2023

@author: pyess
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from TKfunctions.forward_models.onedim.flow import onedim_onecomp_flow_diff as fmod_f
from TKfunctions.forward_models.onedim.loadcase import load_velocity_1d1c as load_kinetics_v
from TKfunctions.utilities import mkdir_p as mkdir
import seaborn as sns

def convertTtoL_oneComp(value_tuple):

    from numpy import zeros_like

    f1, D, v, Bl, Br = value_tuple

    u = zeros_like(D)

    for i in range(0,len(D)):
        if f1 < 0:
            if i == len(D)-1:
                u[i] = f1/v[i-1]
            else:
                u[i] = f1/v[i]
        elif f1 > 0:
            if i == 0:
                u[i] = f1/v[i]
            else:
                u[i] = f1/v[i-1]
        else:
            u[i]=0

    output = [u, D, Bl, Br]

    return output

timepoints = np.arange(80,82,2)
ouput_time = '80'

cases = [1,2,3]
guess = [1,2,3,4,5,6]
colors = sns.color_palette("tab10", n_colors=12).as_hex()

colors = colors[4:10]
        
#%% Different guess case plot
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable 
from scipy import stats

timepoints = np.arange(80,82,2)
ouput_time = '80'

cases = [1,2,3]
guess = [1,2,3,4,5,6]
colors = sns.color_palette("tab10", n_colors=12).as_hex()

colors = colors[4:10]
k=0

m=0

cmap = mpl.cm.viridis
plt.rcParams.update({'font.size':30})

cm = 1/2.54  # centimeters in inches
two_column = 17.56 #cm

fig_height = (10*3)*cm*2
fig_width = (10*3.5)*cm*2
row=3
col=3
fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height),sharey=True,sharex=True)


for c in cases:  
    k = 0   
    for g in guess:
        opt_vals = np.load("1d1c_case{}_guess{}_0snr_0.0001step_2Dt_10000maxiter/FinalOptimisationValues_{}.0s.npz".format(c,g,ouput_time),allow_pickle=True)
        args = np.load("1d1c_case{}_guess{}_0snr_0.0001step_2Dt_10000maxiter/args_in_opt.npz".format(c,g),allow_pickle=True)
        
        value_tuple = opt_vals['value_tuple']
        x_fine = args['x_fine']
        xc_fine = args['xc_fine']
        tsamp=args['tsamp']
        i = opt_vals['i']
        t=args['t']
        C_tiss = args['C_tiss']
        C_meas = args['C_meas']
        meas_res = args['meas_res']
        int_res = args['int_res']
        T = tsamp[-1]
        Lx = x_fine[-1]
        dx,dy,dz=meas_res
        dxyz=dx*dy*dz
        
        x_res = np.linspace(0,Lx,int(Lx/dx)+1)
        
        f_inv, D_inv, v_inv, Cp_inv, Cn_inv = value_tuple
        tissue_params = convertTtoL_oneComp(value_tuple)
        u_inv, D_inv, Cp_inv, Cn_inv = tissue_params
        
        u, D, Jp, Jn, dt, Nt, t = load_kinetics_v(1,1,c)
        
        Jp_inv = Cp_inv * (f_inv/(dx))
        Jn_inv = Cn_inv * (abs(f_inv)/(dx))

        C_inv = fmod_f(Lx, dx, tsamp[:i], dt, f_inv, D_inv, v_inv,  np.append(Jp_inv,Jp_inv[-1]), np.append(Jn_inv,Jn_inv[-1]))
        Dt = tsamp[1]-tsamp[0]
        skip = int(Dt/0.2)
        C_inv= C_inv[:,0::skip]
        
        diff = 100*((C_inv-C_meas)/np.max(C_meas))
        print('1d1c: Case {}, Guess {}, max {}, min {}, mean {}, stdev {}'.format(c,g,np.max(abs(diff)), np.min(abs(diff)), np.mean(abs(diff)),np.std(abs(diff))))
        
        if g == 1:

            ax[m,0].set(xlabel='t (s)', ylabel='x (cm)')
            ax[m,1].set(xlabel='t (s)')
            ax[m,2].set(xlabel='t (s)')
            
            im1=ax[m,0].pcolor(tsamp,x_fine/10,C_inv,vmin=0,vmax=np.max((C_meas,C_inv)),cmap=cmap)
            divider = make_axes_locatable(ax[m,0])
            cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
            fig.add_axes(cax)
            cb=fig.colorbar(im1, cax = cax, orientation = 'vertical')
            cb.ax.tick_params(labelsize=28)
            cb.set_label('$C~ (mM)$', rotation=90)
            
            im2=ax[m,1].pcolor(tsamp,x_fine/10,C_meas,vmin=0,vmax=np.max((C_meas,C_inv)),cmap=cmap)
            divider = make_axes_locatable(ax[m,1])
            cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
            fig.add_axes(cax)
            cb=fig.colorbar(im2, cax = cax, orientation = 'vertical')
            cb.ax.tick_params(labelsize=28)
            cb.set_label('$C ~(mM)$', rotation=90)
            
            im3=ax[m,2].pcolormesh(tsamp,x_fine/10,100*((C_inv-C_meas)/np.max(C_meas)),vmin=-2.5,vmax=2.5,cmap=mpl.cm.bwr.copy())
             
            cmapE = im3.get_cmap()
            cmapE.set_under(color='#030057')
            cmapE.set_over(color='#9E0000')
            im3.set_cmap(cmapE)
        
            divider = make_axes_locatable(ax[m,2])
            cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
            fig.add_axes(cax)
            cb = fig.colorbar(im3, cax = cax, orientation = 'vertical',extend='both')
            cb.ax.tick_params(labelsize=28)
            cb.set_label('$\Delta C ~(\%)$', rotation=90)
            
            if m==0:
                ax[m,0].title.set_text('$C_{inv}$ \n')
                ax[m,1].title.set_text('$C_{meas}$ \n')
                ax[m,2].title.set_text('$\Delta C/C_{max}$ \n')
            
            m=m+1
plt.figtext(0.95,0.65, '(a)',fontsize=36)
plt.figtext(0.95,0.35, '(b)',fontsize=36)
plt.figtext(0.95,0.05, '(c)',fontsize=36)        
fig.tight_layout()
fig.savefig('1d1c_guess1_conc_allcase.png',dpi=300)
