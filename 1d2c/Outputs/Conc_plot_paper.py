#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 10:09:18 2023

@author: pyess
"""

from mpl_toolkits.axes_grid1 import make_axes_locatable

#%% Different guess case plot
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from TKfunctions.forward_models.onedim.flow import onedim_twocomp_flow_diff as fmod_f
from TKfunctions.forward_models.onedim.loadcase import load_velocity_1d2c as load_kinetics_v
from TKfunctions.forward_models.onedim.loadcase import load_flow_1d2c as load_kinetics_f
from TKfunctions.utilities import mkdir_p as mkdir
from TKfunctions.forward_models.onedim.fm_utilities import convertTtoL
import seaborn as sns

from matplotlib.colors import TwoSlopeNorm

cases = [6,11,3]
guess = [1,2,3,4,5]
m=0

cmap = mpl.cm.viridis


for g in [4]:  
    m = 0 
    plt.rcParams.update({'font.size':30})

    cm = 1/2.54  # centimeters in inches
    two_column = 17.56 #cm
    
    fig_height = (10*3)*cm*2
    fig_width = (10*3.5)*cm*2
    row=3
    col=3
    fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height),sharey=True,sharex=True)
    for c in cases:
        opt_vals = np.load("case{}_guess{}_0snr_0.00005step_2Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,g),allow_pickle=True)
        args = np.load("case{}_guess{}_0snr_0.00005step_2Dt_10000maxiter/args_in_opt.npz".format(c,g),allow_pickle=True)
        
        value_tuple = opt_vals['value_tuple']
        x_fine = args['x_fine']
        xc_fine = args['xc_fine']
        tsamp=args['tsamp']
        i = opt_vals['i']
        t=args['t']
        C_tiss = args['C_tiss']
        C_meas = args['C_meas'][:,:-1]
        Ca = args['Ca_ds']
        Cv=args['Cv_ds']
        meas_res = args['meas_res']
        int_res = args['int_res']
        T = tsamp[-1]
        Lx = x_fine[-1]
        dx,dy,dz=meas_res
        dxyz=dx*dy*dz
        
        x_res = np.linspace(0,Lx,int(Lx/dx)+1)
        
        fa1_inv, fv1_inv,fa_inv,fv_inv, F_inv, Da_inv, Dv_inv,va_inv,vaf_inv, vv_inv, v_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = value_tuple
        tissue_params = convertTtoL(value_tuple)
        ua_inv, uv_inv, Kva_inv, Da_inv, Dv_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = tissue_params
   
        fa, fv, Da, Dv, F, va, va_frac, vv, v, Cpa, Cpv, Cna, Cnv, dt, Nt, t = load_kinetics_f(1,2,c)
        ua, uv, Da, Dv, Kva, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_v(1,2,c)

        Jpa_inv = Cpa_inv * (fa_inv[0]/(dx))
        Jna_inv = Cna_inv * (abs(fa_inv[-1])/(dx))
        Jpv_inv = Cpv_inv * (abs(fv_inv[0])/(dx))
        Jnv_inv = Cnv_inv * (fv_inv[-1]/(dx))

        
        C_inv, Ca_inv, Cv_inv = fmod_f(Lx, dx, tsamp[:i], dt, fa_inv, fv_inv, F_inv, Da_inv, Dv_inv, va_inv, vv_inv,  np.append(Jpa_inv,Jpa_inv[-1]), np.append(Jpv_inv,Jpv_inv[-1]), np.append(Jna_inv,Jna_inv[-1]), np.append(Jnv_inv,Jnv_inv[-1]))
        
        Dt = tsamp[1]-tsamp[0]
        if c==11:
            skip = int(Dt/0.02)
        else:
            skip = int(Dt/0.2)
        C_inv= C_inv[:,0::skip]
        
        diff = 100*((C_inv-C_meas)/np.max(C_meas))
        if c==6 and g==1:
            print('Case, Guess, max, min, mean, stdev')
        print('{}, {}, {}, {}, {}, {}'.format(c,g,np.max((diff)), np.min((diff)), np.mean(abs(diff)),np.std(abs(diff))))

        ax[m,0].set(xlabel='t (s)', ylabel='x (cm)')
        ax[m,1].set(xlabel='t (s)')
        ax[m,2].set(xlabel='t (s)')
        
        if c ==11:
            im1=ax[m,0].pcolormesh(tsamp,x_fine,C_inv,vmin=0,vmax=np.max((C_meas,C_inv)),cmap=cmap)
            divider = make_axes_locatable(ax[m,0])
            cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
            fig.add_axes(cax)
            cb=fig.colorbar(im1, cax = cax, orientation = 'vertical')
            cb.ax.tick_params(labelsize=28)
            cb.set_label('$C~ (mM)$', rotation=90)
            
            im2=ax[m,1].pcolormesh(tsamp,x_fine,C_meas,vmin=0,vmax=np.max((C_meas,C_inv)),cmap=cmap)
            divider = make_axes_locatable(ax[m,1])
            cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
            fig.add_axes(cax)
            cb=fig.colorbar(im2, cax = cax, orientation = 'vertical')
            cb.ax.tick_params(labelsize=28)
            cb.set_label('$C ~(mM)$', rotation=90)

            im3=ax[m,2].pcolormesh(tsamp,x_fine,100*((C_inv-C_meas)/np.max(C_meas)),vmin=-5,vmax=5,cmap=mpl.cm.bwr.copy())
        else:
            im1=ax[m,0].pcolormesh(tsamp,x_fine/10,C_inv,vmin=0,vmax=np.max((C_meas,C_inv)),cmap=cmap)
            divider = make_axes_locatable(ax[m,0])
            cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
            fig.add_axes(cax)
            cb=fig.colorbar(im1, cax = cax, orientation = 'vertical')
            cb.ax.tick_params(labelsize=28)
            cb.set_label('$C~ (mM)$', rotation=90)
            
            im2=ax[m,1].pcolormesh(tsamp,x_fine/10,C_meas,vmin=0,vmax=np.max((C_meas,C_inv)),cmap=cmap)
            divider = make_axes_locatable(ax[m,1])
            cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
            fig.add_axes(cax)
            cb=fig.colorbar(im2, cax = cax, orientation = 'vertical')
            cb.ax.tick_params(labelsize=28)
            cb.set_label('$C ~(mM)$', rotation=90)

            im3=ax[m,2].pcolormesh(tsamp,x_fine/10,100*((C_inv-C_meas)/np.max(C_meas)),vmin=-5,vmax=5,cmap=mpl.cm.bwr.copy())
        
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
    fig.savefig('Figures/1d2c_guess{}_conc_allcase_newcolors.pdf'.format(g),dpi=300)
    

#%% conc plots Ca
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from TKfunctions.forward_models.onedim.flow import onedim_twocomp_flow_diff as fmod_f
from TKfunctions.forward_models.onedim.loadcase import load_velocity_1d2c as load_kinetics_v
from TKfunctions.forward_models.onedim.loadcase import load_flow_1d2c as load_kinetics_f
from TKfunctions.utilities import mkdir_p as mkdir
from TKfunctions.forward_models.onedim.fm_utilities import convertTtoL
import seaborn as sns
cases = [6,11,3]
guess = [1,2,3,4,5]
m=0

cmap = mpl.cm.viridis


for g in guess:  
    m = 0 
    plt.rcParams.update({'font.size':30})

    cm = 1/2.54  # centimeters in inches
    two_column = 17.56 #cm
    
    fig_height = (10*3)*cm*2
    fig_width = (10*3.5)*cm*2
    row=3
    col=3
    fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height),sharey=True,sharex=True)
    for c in cases:
        opt_vals = np.load("case{}_guess{}_0snr_0.00005step_2Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,g),allow_pickle=True)
        args = np.load("case{}_guess{}_0snr_0.00005step_2Dt_10000maxiter/args_in_opt.npz".format(c,g),allow_pickle=True)
        
        value_tuple = opt_vals['value_tuple']
        x_fine = args['x_fine']
        xc_fine = args['xc_fine']
        tsamp=args['tsamp']
        i = opt_vals['i']
        t=args['t']
        C_tiss = args['C_tiss']
        C_meas = args['C_meas'][:,:-1]
        Ca = args['Ca_ds'][:,:-1]
        Cv=args['Cv_ds'][:,:-1]
        meas_res = args['meas_res']
        int_res = args['int_res']
        T = tsamp[-1]
        Lx = x_fine[-1]
        dx,dy,dz=meas_res
        dxyz=dx*dy*dz
        
        x_res = np.linspace(0,Lx,int(Lx/dx)+1)
        
        fa1_inv, fv1_inv,fa_inv,fv_inv, F_inv, Da_inv, Dv_inv,va_inv,vaf_inv, vv_inv, v_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = value_tuple
        tissue_params = convertTtoL(value_tuple)
        ua_inv, uv_inv, Kva_inv, Da_inv, Dv_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = tissue_params
        
        fa, fv, Da, Dv, F, va, va_frac, vv, v, Cpa, Cpv, Cna, Cnv, dt, Nt, t = load_kinetics_f(1,2,c)
        ua, uv, Da, Dv, Kva, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_v(1,2,c)
        
        Jpa_inv = Cpa_inv * (fa_inv[0]/(dx))
        Jna_inv = Cna_inv * (abs(fa_inv[-1])/(dx))
        Jpv_inv = Cpv_inv * (abs(fv_inv[0])/(dx))
        Jnv_inv = Cnv_inv * (fv_inv[-1]/(dx))
        
        C_inv, Ca_inv, Cv_inv = fmod_f(Lx, dx, tsamp[:i], dt, fa_inv, fv_inv, F_inv, Da_inv, Dv_inv, va_inv, vv_inv,  np.append(Jpa_inv,Jpa_inv[-1]), np.append(Jpv_inv,Jpv_inv[-1]), np.append(Jna_inv,Jna_inv[-1]), np.append(Jnv_inv,Jnv_inv[-1]))
        Dt = tsamp[1]-tsamp[0]
        if c==11:
            skip = int(Dt/0.02)
        else:
            skip = int(Dt/0.2)
        Ca_inv= Ca_inv[:,0::skip]

        ax[m,0].set(xlabel='t (s)', ylabel='x (cm)')
        ax[m,1].set(xlabel='t (s)')
        ax[m,2].set(xlabel='t (s)')
        
        maxCa = np.max((np.max(Ca),np.max(Ca_inv)))
        if c==11:
            x_plot = x_fine*1
        else:
            x_plot = x_fine/10
        
        im1=ax[m,0].pcolor(tsamp,x_plot,Ca_inv,vmin=0,vmax=maxCa,cmap=cmap)
        divider = make_axes_locatable(ax[m,0])
        cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
        fig.add_axes(cax)
        cb=fig.colorbar(im1, cax = cax, orientation = 'vertical')
        cb.ax.tick_params(labelsize=28)
        cb.set_label('$C^{a}~ (mM)$', rotation=90)
        
        im2=ax[m,1].pcolor(tsamp,x_plot,Ca,vmin=0,vmax=maxCa,cmap=cmap)
        divider = make_axes_locatable(ax[m,1])
        cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
        fig.add_axes(cax)
        cb=fig.colorbar(im2, cax = cax, orientation = 'vertical')
        cb.ax.tick_params(labelsize=28)
        cb.set_label('$C^{a} ~(mM)$', rotation=90)
        
        im3=ax[m,2].pcolor(tsamp,x_plot,100*((Ca_inv-Ca)/np.max(Ca)),vmax=np.max(100*(Ca_inv-Ca)/np.max(Ca)),cmap=cmap)
        divider = make_axes_locatable(ax[m,2])
        cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
        fig.add_axes(cax)
        cb = fig.colorbar(im3, cax = cax, orientation = 'vertical')
        cb.ax.tick_params(labelsize=28)
        cb.set_label('$\Delta C^{a} ~(\%)$', rotation=90)
        
        if m==0:
            ax[m,0].title.set_text('$C^{a}_{inv}$ \n')
            ax[m,1].title.set_text('$C^{a}_{meas}$ \n')
            ax[m,2].title.set_text('$\Delta C^{a}/C^{a}_{max}$ \n')
        
        m=m+1
    
    fig.tight_layout()
    fig.savefig('Figures/1d2c_guess{}_concCa_allcase.tiff'.format(g),dpi=100)
    
#%% conc plots Cv
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from TKfunctions.forward_models.onedim.flow import onedim_twocomp_flow_diff as fmod_f
from TKfunctions.forward_models.onedim.loadcase import load_velocity_1d2c as load_kinetics_v
from TKfunctions.forward_models.onedim.loadcase import load_flow_1d2c as load_kinetics_f
from TKfunctions.utilities import mkdir_p as mkdir
from TKfunctions.forward_models.onedim.fm_utilities import convertTtoL
import seaborn as sns
cases = [6,11,3]
guess = [1,2,3,4,5]
m=0

cmap = mpl.cm.viridis


for g in guess:  
    m = 0 
    plt.rcParams.update({'font.size':30})

    cm = 1/2.54  # centimeters in inches
    two_column = 17.56 #cm
    
    fig_height = (10*3)*cm*2
    fig_width = (10*3.5)*cm*2
    row=3
    col=3
    fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height),sharey=True,sharex=True)
    for c in cases:
        opt_vals = np.load("case{}_guess{}_0snr_0.00005step_2Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,g),allow_pickle=True)
        args = np.load("case{}_guess{}_0snr_0.00005step_2Dt_10000maxiter/args_in_opt.npz".format(c,g),allow_pickle=True)
        
        value_tuple = opt_vals['value_tuple']
        x_fine = args['x_fine']
        xc_fine = args['xc_fine']
        tsamp=args['tsamp']
        i = opt_vals['i']
        t=args['t']
        C_tiss = args['C_tiss']
        C_meas = args['C_meas'][:,:-1]
        Ca = args['Ca_ds'][:,:-1]
        Cv=args['Cv_ds'][:,:-1]
        meas_res = args['meas_res']
        int_res = args['int_res']
        T = tsamp[-1]
        Lx = x_fine[-1]
        dx,dy,dz=meas_res
        dxyz=dx*dy*dz
        
        x_res = np.linspace(0,Lx,int(Lx/dx)+1)
        
        fa1_inv, fv1_inv,fa_inv,fv_inv, F_inv, Da_inv, Dv_inv,va_inv,vaf_inv, vv_inv, v_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = value_tuple
        tissue_params = convertTtoL(value_tuple)
        ua_inv, uv_inv, Kva_inv, Da_inv, Dv_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = tissue_params
        
        fa, fv, Da, Dv, F, va, va_frac, vv, v, Cpa, Cpv, Cna, Cnv, dt, Nt, t = load_kinetics_f(1,2,c)
        ua, uv, Da, Dv, Kva, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_v(1,2,c)
        
        Jpa_inv = Cpa_inv * (fa_inv[0]/(dx))
        Jna_inv = Cna_inv * (abs(fa_inv[-1])/(dx))
        Jpv_inv = Cpv_inv * (abs(fv_inv[0])/(dx))
        Jnv_inv = Cnv_inv * (fv_inv[-1]/(dx))
        
        
        C_inv, Ca_inv, Cv_inv = fmod_f(Lx, dx, tsamp[:i], dt, fa_inv, fv_inv, F_inv, Da_inv, Dv_inv, va_inv, vv_inv,  np.append(Jpa_inv,Jpa_inv[-1]), np.append(Jpv_inv,Jpv_inv[-1]), np.append(Jna_inv,Jna_inv[-1]), np.append(Jnv_inv,Jnv_inv[-1]))
        
        Dt = tsamp[1]-tsamp[0]
        if c==11:
            skip = int(Dt/0.02)
        else:
            skip = int(Dt/0.2)
        Cv_inv= Cv_inv[:,0::skip]
        maxCv = np.max((np.max(Cv),np.max(Cv_inv)))
        if c==11:
            x_plot = x_fine*1
        else:
            x_plot = x_fine/10

        ax[m,0].set(xlabel='t (s)', ylabel='x (cm)')
        ax[m,1].set(xlabel='t (s)')
        ax[m,2].set(xlabel='t (s)')
        
        im1=ax[m,0].pcolor(tsamp,x_plot,Cv_inv,vmin=0,vmax=maxCv,cmap=cmap)
        divider = make_axes_locatable(ax[m,0])
        cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
        fig.add_axes(cax)
        cb=fig.colorbar(im1, cax = cax, orientation = 'vertical')
        cb.ax.tick_params(labelsize=28)
        cb.set_label('$C^{v}~ (mM)$', rotation=90)
        
        im2=ax[m,1].pcolor(tsamp,x_plot,Cv,vmin=0,vmax=maxCv,cmap=cmap)
        divider = make_axes_locatable(ax[m,1])
        cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
        fig.add_axes(cax)
        cb=fig.colorbar(im2, cax = cax, orientation = 'vertical')
        cb.ax.tick_params(labelsize=28)
        cb.set_label('$C^{v} ~(mM)$', rotation=90)
        
        im3=ax[m,2].pcolor(tsamp,x_plot,100*((Cv_inv-Cv)/np.max(Cv)),vmax=np.max(100*(Cv_inv-Cv)/np.max(Cv)),cmap=cmap)
        divider = make_axes_locatable(ax[m,2])
        cax = divider.new_horizontal(size='5%', pad=0.1, pack_start = False)
        fig.add_axes(cax)
        cb = fig.colorbar(im3, cax = cax, orientation = 'vertical')
        cb.ax.tick_params(labelsize=28)
        cb.set_label('$\Delta C^{v} ~(\%)$', rotation=90)
        
        if m==0:
            ax[m,0].title.set_text('$C^{v}_{inv}$ \n')
            ax[m,1].title.set_text('$C^{v}_{meas}$ \n')
            ax[m,2].title.set_text('$\Delta C^{v}/C^{v}_{max}$ \n')
        
        m=m+1
    
    fig.tight_layout()
    fig.savefig('Figures/1d2c_guess{}_concCv_allcase.tiff'.format(g),dpi=100)
