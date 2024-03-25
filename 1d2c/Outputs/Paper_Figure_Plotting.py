#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 10:24:24 2022

@author: pyess
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from TKfunctions.forward_models.onedim.flow import onedim_twocomp_flow_diff as fmod_f
from TKfunctions.forward_models.onedim.loadcase import load_velocity_1d2c as load_kinetics_v
from TKfunctions.forward_models.onedim.loadcase import load_flow_1d2c as load_kinetics_f
from TKfunctions.utilities import mkdir_p as mkdir
from TKfunctions.forward_models.onedim.fm_utilities import convertTtoL
import seaborn as sns
timepoints = np.arange(80,82,2)


cases = [6,11,3]
guess = [1,2,3,4,5]
colors = sns.color_palette("tab10", n_colors=12).as_hex()

colors = colors[4:10]
k=0
plt.rcParams.update({'font.size':34})

cm = 1/2.54  # centimeters in inches
two_column = 17.56 #cm
#%%
plt.rcParams.update({'font.size':34})

cm = 1/2.54  # centimeters in inches
two_column = 17.56 #cm

fig_height = (10*5)*cm*2
fig_width = (10*3)*cm*2
row=5
col=3
fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height))
ax[0,0].set(xlabel='x (cm)', ylabel='$K^{va}$ (1/s)')
ax[1,0].set(xlabel='x (cm)', ylabel='$u^{a}$ (cm/s)')
ax[2,0].set(xlabel='x (cm)', ylabel='$u^{v}$ (cm/s)',ylim=([-10,10]))
ax[3,0].set(xlabel='t (s)', ylabel='$J^{a}_{+}$ (mM/s)')
ax[4,0].set(xlabel='t (s)', ylabel='$J^{a}_{-}$ (mM/s)')
m=-1
for c in cases:  
    k = 0
    m=m+1    
    for g in guess:
        opt_vals = np.load("case{}_guess{}_0snr_0.00005step_2Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,g),allow_pickle=True)
        args = np.load("case{}_guess{}_0snr_0.00005step_2Dt_10000maxiter/args_in_opt.npz".format(c,g),allow_pickle=True)
    
        value_tuple = opt_vals['value_tuple']
        x_fine = args['x_fine']
        xc_fine = args['xc_fine']
        tsamp=args['tsamp']
        i = opt_vals['i']
        t=args['t']
        C_tiss = args['C_tiss']
        C_meas = args['C_meas']
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
        
        fa, fv, Da, Dv, F, va, va_frac, vv, v, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_f(1,2,c)
        ua, uv, Da, Dv, Kva, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_v(1,2,c)
        
        Jpa_inv = Cpa_inv * (fa_inv[0]/(dx))
        Jna_inv = Cna_inv * (abs(fa_inv[-1])/(dx))
        Jpv_inv = Cpv_inv * (abs(fv_inv[0])/(dx))
        Jnv_inv = Cnv_inv * (fv_inv[-1]/(dx))
        
        C_inv, Ca_inv, Cv_inv = fmod_f(Lx, dx, tsamp[:i], dt, fa_inv, fv_inv, F_inv, Da_inv, Dv_inv, va_inv, vv_inv,  np.append(Jpa_inv,Jpa_inv[-1]), np.append(Jpv_inv,Jpv_inv[-1]), np.append(Jna_inv,Jna_inv[-1]), np.append(Jnv_inv,Jnv_inv[-1]))
        Dt = tsamp[1]-tsamp[0]
        skip = int(Dt/0.2)
        C_inv= C_inv[:,0::skip]
        
        
        if g == 1:                    
            if c==11:
                ax[0,m].set(xlim=([0,Lx]),xlabel='x (cm)')
                ax[1,m].set(xlim=([0,Lx]),xlabel='x (cm)',ylim=([-20,20]))
                ax[2,m].set(xlim=([0,Lx]),xlabel='x (cm)',ylim=([-10,10]))
                ax[3,m].set(xlim=([0,T]),xlabel='t (s)')
                ax[4,m].set(xlim=([0,T]),xlabel='t (s)')
                
                ax[0,m].plot(xc_fine,Kva,'k')
                
                ax[1,m].plot(x_fine,ua,'k')
                ax[2,m].plot(x_fine,uv,'k')
                ax[3,m].plot(t,Jpa,'k')
                ax[4,m].plot(t,Jna,'k')
            else:
                
                ax[0,m].set(xlim=([0,Lx/10]),xlabel='x (cm)')
                ax[1,m].set(xlim=([0,Lx/10]),xlabel='x (cm)',ylim=([-4,4]))
                ax[2,m].set(xlim=([0,Lx/10]),xlabel='x (cm)',ylim=([-4,4]))
                ax[3,m].set(xlim=([0,T]),xlabel='t (s)')
                ax[4,m].set(xlim=([0,T]),xlabel='t (s)')
                
                ax[0,m].plot(xc_fine/10,Kva,'k')
                
                ax[1,m].plot(x_fine/10,ua/10,'k')
                ax[2,m].plot(x_fine/10,uv/10,'k')
                ax[3,m].plot(t,Jpa,'k')
                ax[4,m].plot(t,Jna,'k')
        
        if c==11:
            ax[0,m].plot(xc_fine,Kva_inv,colors[k], marker='o', linestyle='',label=g)
            ax[1,m].plot(x_fine,ua_inv,colors[k], marker='o', linestyle='')
            ax[2,m].plot(x_fine,uv_inv,colors[k], marker='o', linestyle='')
            ax[3,m].plot(tsamp[:i-1],Jpa_inv,colors[k], marker='o', linestyle='')
            ax[4,m].plot(tsamp[:i-1], Jna_inv,colors[k], marker='o', linestyle='')
        else:
            ax[0,m].plot(xc_fine/10,Kva_inv,colors[k], marker='o', linestyle='',label=g)
            ax[1,m].plot(x_fine/10,ua_inv/10,colors[k], marker='o', linestyle='')
            ax[2,m].plot(x_fine/10,uv_inv/10,colors[k], marker='o', linestyle='')
            ax[3,m].plot(tsamp[:i-1],Jpa_inv,colors[k], marker='o', linestyle='')
            ax[4,m].plot(tsamp[:i-1], Jna_inv,colors[k], marker='o', linestyle='')
            
        ax[0,m].set(title='Case {}'.format(m+1))
        k+=1
    legend = ax[0,2].legend(title="Guess",
                    loc='upper right', fontsize=24, fancybox=True, ncol=1)

        
fig.tight_layout()
fig.savefig('Figures/paper_1d2c_all_guesses_vel_param_legend_new.png', bbox_inches='tight')

#%%
snr = [5,10,15,20]
repeats= [1,2,3,4,5]

Kva_all = []
Jpa_all  = []
Jna_all  = []
ua_all  = []
uv_all  = []
all_avg_pdiff = []
cost_all = []
panel = 0

for c in cases:
    panel =0
    for s in snr:
        avg_pdiff = []
        cost_r = []
        k=0   
        Kva_all = []
        Jpa_all  = []
        Jna_all  = []
        ua_all  = []
        uv_all  = []
        for r in repeats:
            if c==11:
                opt_vals = np.load("case{}_guess1_r{}_{}snr_0.00005step_2Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,r,s),allow_pickle=True)
                args = np.load("case{}_guess1_r{}_{}snr_0.00005step_2Dt_10000maxiter/args_in_opt.npz".format(c,r,s),allow_pickle=True)
            else:
                opt_vals = np.load("case{}_guess1_{}snr_r{}_0.00005step_2Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,s,r),allow_pickle=True)
                args = np.load("case{}_guess1_{}snr_r{}_0.00005step_2Dt_10000maxiter/args_in_opt.npz".format(c,s,r),allow_pickle=True)
            
            value_tuple = opt_vals['value_tuple']
            x_fine = args['x_fine']
            xc_fine = args['xc_fine']
            tsamp=args['tsamp']
            i = opt_vals['i']
            t=args['t']
            C_tiss = args['C_tiss']
            C_meas = args['C_meas']
            Ca = args['Ca_ds']
            Cv=args['Cv_ds']
            meas_res = args['meas_res']
            int_res = args['int_res']
            T = tsamp[-1]
            Lx = x_fine[-1]
            dx,dy,dz=meas_res
            dxyz=dx*dy*dz
            
            x_res = np.linspace(0,Lx,int(Lx/dx)+1)
            xc_res = np.linspace(dx/2,Lx-dx/2,int(Lx/dx))
            
            fa1_inv, fv1_inv,fa_inv,fv_inv, F_inv, Da_inv, Dv_inv,va_inv,vaf_inv, vv_inv, v_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = value_tuple
            tissue_params = convertTtoL(value_tuple)
            ua_inv, uv_inv, Kva_inv, Da_inv, Dv_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = tissue_params
            
            fa, fv, Da, Dv, F, va, va_frac, vv, v, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_f(1,2,c)
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
            cost = (1/(np.array(C_inv.shape).prod(axis=0)))*abs(np.sqrt(np.sum(np.power(C_meas[:,:-1]-C_inv,2))))
            cost_r += [cost,]


            if r == 1 and s==5:
                plt.rcParams.update({'font.size':34})
                
                cm = 1/2.54  # centimeters in inches
                two_column = 17.56 #cm
                
                fig_height = (10*5)*cm*2
                fig_width = (10*4)*cm*2
                row=5
                col=4
                fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height),sharey='row')
                ax[0,panel].set(xlabel='x (cm)', ylabel='$K^{va}$ (1/s)')
                ax[1,panel].set(xlabel='x (cm)', ylabel='$u^{a}$ (cm/s)')
                ax[2,panel].set(xlabel='x (cm)', ylabel='$u^{v}$ (cm/s)')
                ax[3,panel].set(xlabel='t (s)', ylabel='$J^{a}_{+}$ (mM/s)')
                ax[4,panel].set(xlabel='t (s)', ylabel='$J^{a}_{-}$ (mM/s)')
            ax[0,panel].set(xlabel='x (cm)')
            ax[1,panel].set(xlabel='x (cm)')
            ax[2,panel].set(xlabel='x (cm)')
            ax[3,panel].set(xlabel='t (s)')
            ax[4,panel].set(xlabel='t (s)')


            param = Kva
            param_inv = Kva_inv
            avg_pdiff += [np.mean(abs(param_inv-param)),]
            
            if c==11:
                ax[0,panel].plot(xc_fine,Kva,'k')
                ax[0,panel].plot(xc_fine,Kva_inv,colors[k],linestyle='None',marker='o')
                ax[0,panel].set(title='SNR: {}'.format(s))
                
                ax[1,panel].plot(x_fine,ua,'k')
                ax[1,panel].plot(x_fine,ua_inv,colors[k],linestyle='None',marker='o')
                ax[1,panel].set(ylim=([-25,25]))
                
                ax[2,panel].plot(x_fine,uv,'k')
                ax[2,panel].plot(x_fine,uv_inv,colors[k],linestyle='None',marker='o')
                ax[2,panel].set(ylim=([-10,10]))
                
                ax[3,panel].plot(t,Jpa,'k')
                ax[3,panel].plot(tsamp[:i-1],Jpa_inv,colors[k],linestyle='None',marker='o')
                
                ax[4,panel].plot(t,Jna,'k')
                ax[4,panel].plot(tsamp[:i-1], Jna_inv,colors[k],linestyle='None',marker='o',label=r)
            else:
            
                ax[0,panel].plot(xc_fine/10,Kva,'k')
                ax[0,panel].plot(xc_fine/10,Kva_inv,colors[k],linestyle='None',marker='o')
                ax[0,panel].set(title='SNR: {}'.format(s))
                
                ax[1,panel].plot(x_fine/10,ua/10,'k')
                ax[1,panel].plot(x_fine/10,ua_inv/10,colors[k],linestyle='None',marker='o')
                
                ax[2,panel].plot(x_fine/10,uv/10,'k')
                ax[2,panel].plot(x_fine/10,uv_inv/10,colors[k],linestyle='None',marker='o')
                ax[2,panel].set(ylim=([-3,3]))
                
                ax[3,panel].plot(t,Jpa,'k')
                ax[3,panel].plot(tsamp[:i-1],Jpa_inv,colors[k],linestyle='None',marker='o')
                
                ax[4,panel].plot(t,Jna,'k')
                ax[4,panel].plot(tsamp[:i-1], Jna_inv,colors[k],linestyle='None',marker='o',label=r)
            
            	
            Kva_all += [Kva_inv,]
            Jpa_all += [Jpa_inv,]
            Jna_all += [Jna_inv,]
            ua_all += [ua_inv,]
            uv_all += [uv_inv,]
            k+=1
        
        if c==11:
            all_avg_pdiff += [avg_pdiff,]
            cost_all += [cost_r,]
            Kva_all = np.asarray(Kva_all)
            Jpa_all  = np.asarray(Jpa_all)
            Jna_all  = np.asarray(Jna_all)
            ua_all  = np.asarray(ua_all)
            uv_all  = np.asarray(uv_all)
            
            Kva_m = np.mean(Kva_all,axis=0)
            Jpa_m  = np.mean(Jpa_all,axis=0)
            Jna_m  = np.mean(Jna_all,axis=0)
            ua_m  = np.mean(ua_all,axis=0)
            uv_m  = np.mean(uv_all,axis=0)
            
            Kva_ci = (1.96*np.std(Kva_all,axis=0))/np.sqrt(5)
            Jpa_ci  = (1.96*np.std(Jpa_all,axis=0))/np.sqrt(5)
            Jna_ci  = (1.96*np.std(Jna_all,axis=0))/np.sqrt(5)
            ua_ci  = (1.96*np.std(ua_all,axis=0))/np.sqrt(5)
            uv_ci  = (1.96*np.std(uv_all,axis=0))/np.sqrt(5)
        
                          
            ax[0,panel].fill_between(xc_fine,Kva_m-Kva_ci,Kva_m+Kva_ci)
            ax[1,panel].fill_between(x_fine,ua_m-ua_ci,ua_m+ua_ci)
            ax[2,panel].fill_between(x_fine,uv_m-uv_ci,uv_m+uv_ci)
            ax[3,panel].fill_between(tsamp[:i-1],Jpa_m-Jpa_ci,Jpa_m+Jpa_ci)
            ax[4,panel].fill_between(tsamp[:i-1],Jna_m-Jna_ci,Jna_m+Jna_ci)
            
        else:
            all_avg_pdiff += [avg_pdiff,]
            cost_all += [cost_r,]
            Kva_all = np.asarray(Kva_all)
            Jpa_all  = np.asarray(Jpa_all)
            Jna_all  = np.asarray(Jna_all)
            ua_all  = np.asarray(ua_all)/10
            uv_all  = np.asarray(uv_all)/10
            
            Kva_m = np.mean(Kva_all,axis=0)
            Jpa_m  = np.mean(Jpa_all,axis=0)
            Jna_m  = np.mean(Jna_all,axis=0)
            ua_m  = np.mean(ua_all,axis=0)
            uv_m  = np.mean(uv_all,axis=0)
            
            Kva_ci = (1.96*np.std(Kva_all,axis=0))/np.sqrt(5)
            Jpa_ci  = (1.96*np.std(Jpa_all,axis=0))/np.sqrt(5)
            Jna_ci  = (1.96*np.std(Jna_all,axis=0))/np.sqrt(5)
            ua_ci  = (1.96*np.std(ua_all,axis=0))/np.sqrt(5)
            uv_ci  = (1.96*np.std(uv_all,axis=0))/np.sqrt(5)
        
                          
            ax[0,panel].fill_between(xc_fine/10,Kva_m-Kva_ci,Kva_m+Kva_ci)
            ax[1,panel].fill_between(x_fine/10,ua_m-ua_ci,ua_m+ua_ci)
            ax[2,panel].fill_between(x_fine/10,uv_m-uv_ci,uv_m+uv_ci)
            ax[3,panel].fill_between(tsamp[:i-1],Jpa_m-Jpa_ci,Jpa_m+Jpa_ci)
            ax[4,panel].fill_between(tsamp[:i-1],Jna_m-Jna_ci,Jna_m+Jna_ci)
        panel += 1
    legend = ax[4,3].legend(title="Repeat",loc='upper right', fontsize=24, fancybox=True, ncol=1)
    plt.tight_layout()
    fig.savefig('Figures/paper_1d2c_snr_vel_param_case{}_legend.png'.format(c),bbox_inches='tight')

#%% Layout 4 mixed case separate params
plt.rcParams.update({'font.size':34})

cm = 1/2.54  # centimeters in inches
two_column = 17.56 #cm



snr = [5,10,15,20]
repeats= [1,2,3,4,5]

panel = 0
p_case = [0,2,3,0,0,1,0,0,0,0,11]

Kva_list_all = []
ua_list_all = []
uv_list_all = []
Jpa_list_all = []
Jna_list_all = []
Ja_list_all = []
for c in [6,11,3]:
    snr_params_all = []
    Kva_list = []
    ua_list = []
    uv_list = []
    Jpa_list = []
    Jna_list = []
    Ja_list = []
    for s in snr:
        ua_diff_all = []   
        uv_diff_all = []
        Kva_diff_all =  []  
        Jpa_diff_all =  [] 
        Jna_diff_all =  []
        Ja_diff_all = []
        k=0   
        for r in repeats:
            if c==11:
                opt_vals = np.load("case{}_guess1_r{}_{}snr_0.00005step_2Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,r,s),allow_pickle=True)
                args = np.load("case{}_guess1_r{}_{}snr_0.00005step_2Dt_10000maxiter/args_in_opt.npz".format(c,r,s),allow_pickle=True)
            else:
                opt_vals = np.load("case{}_guess1_{}snr_r{}_0.00005step_2Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,s,r),allow_pickle=True)
                args = np.load("case{}_guess1_{}snr_r{}_0.00005step_2Dt_10000maxiter/args_in_opt.npz".format(c,s,r),allow_pickle=True)
            
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
            xc_res = np.linspace(dx/2,Lx-dx/2,int(Lx/dx))
            
            fa1_inv, fv1_inv,fa_inv,fv_inv, F_inv, Da_inv, Dv_inv,va_inv,vaf_inv, vv_inv, v_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = value_tuple
            tissue_params = convertTtoL(value_tuple)
            ua_inv, uv_inv, Kva_inv, Da_inv, Dv_inv, Cpa_inv, Cpv_inv, Cna_inv, Cnv_inv = tissue_params
            
            fa, fv, Da, Dv, F, va, va_frac, vv, v, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_f(1,2,c)
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
            Jpa = Jpa[0::skip]
            Jna = Jna[0::skip]
            cost = (1/(np.array(C_inv.shape).prod(axis=0)))*abs(np.sqrt(np.sum(np.power(C_meas-C_inv,2))))
            
            conc_diff = 100*((C_inv-C_meas)/np.max(C_meas))
            

            ua_diff = 100*(abs(ua-ua_inv)/np.mean(abs(ua)))
            uv_diff = 100*(abs(uv-uv_inv)/np.mean(abs(uv)))
            Kva_diff = 100*(abs(Kva-Kva_inv)/np.mean(abs(Kva)))
            Jpa_diff = 100*(abs(Jpa[:-2]-Jpa_inv)/np.mean(abs(Jpa)))
            Jna_diff = 100*(abs(Jna[:-2]-Jna_inv)/np.mean(abs(Jna)))
            
            Ja_diff=np.append(Jpa_diff,Jna_diff)
            
            ua_diff_all += [abs(ua_diff)] 
            uv_diff_all += [abs(uv_diff)]
            Kva_diff_all += [abs(Kva_diff)]
            Jpa_diff_all += [abs(Jpa_diff)]
            Jna_diff_all += [abs(Jna_diff)] 
            Ja_diff_all += [abs(Ja_diff)]                
   
        
        ua_diff_all = np.asarray(ua_diff_all).flatten()
        uv_diff_all = np.asarray(uv_diff_all).flatten()
        Kva_diff_all = np.asarray(Kva_diff_all).flatten()
        Jpa_diff_all = np.asarray(Jpa_diff_all).flatten()
        Jna_diff_all = np.asarray(Jna_diff_all).flatten()
        Ja_diff_all = np.asarray(Ja_diff_all).flatten()
        
        if c==6 and r==5 and s==5:
            print('Case, SNR, ua mean, ua std, uv mean, uv std, Kva mean, Kva std, Ja mean, Ja std ')
        print('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}'.format(c,s,np.mean(ua_diff_all),np.std(ua_diff_all),np.mean(uv_diff_all),np.std(uv_diff_all),np.mean(Kva_diff_all),np.std(Kva_diff_all),np.mean((Jpa_diff_all,Jna_diff_all)),np.std((Jpa_diff_all,Jna_diff_all))))
         
        
        Kva_list += [Kva_diff_all,]
        ua_list += [ua_diff_all,]
        uv_list += [uv_diff_all,]
        Jpa_list += [Jpa_diff_all,]
        Jna_list += [Jna_diff_all,]
        Ja_list += [Ja_diff_all,]

        all_param = np.concatenate((ua_diff_all,uv_diff_all,Kva_diff_all,Jpa_diff_all,Jna_diff_all))
        snr_params_all += [all_param,]
    
    Kva_list_all += [Kva_list,]
    ua_list_all += [ua_list,]
    uv_list_all += [uv_list,]
    Jpa_list_all += [Jpa_list,]
    Jna_list_all += [Jna_list,]
    Ja_list_all += [Ja_list,]
        
    fig_height = (10*1)*cm*2
    fig_width = (10*5)*cm*2
    row=1
    col=5
    fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height))
    
    medianprops = dict(linestyle='-', linewidth=3, color='black')
    meanpointprops = dict(marker='o', markersize=8, markeredgecolor='black',markerfacecolor='firebrick')
    bplot1=ax[0].boxplot(Kva_list,vert = True,whis=(5,95),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
    ax[0].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$K^{va}$ (1/s)')
    ax[0].set(ylabel='Absolute Error')
    ax[0].set(xlabel='SNR')
    bplot2=ax[1].boxplot(ua_list,vert = True,whis=(5,95),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
    ax[1].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$u^{a}$ (cm/s)')
    ax[1].set(xlabel='SNR')
    bplot3=ax[2].boxplot(uv_list,vert = True,whis=(5,95),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
    ax[2].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$u^{v}$ (cm/s)')
    ax[2].set(xlabel='SNR')
    bplot4=ax[3].boxplot(Jpa_list,vert = True,whis=(5,95),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
    ax[3].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$J^{a}_{+}$ (mM)')
    ax[3].set(xlabel='SNR')
    if np.sum(Jna_list==0):
        bplot5=ax[4].boxplot([0,0,0,0],vert = True,whis=(5,95),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
    bplot5=ax[4].boxplot(Jna_list,vert = True,whis=(5,95),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
    ax[4].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$J^{a}_{-}$ (mM/s)')
    ax[4].set(xlabel='SNR')
    fig.suptitle('Case {}'.format(p_case[c-1]))
    
    plt.tight_layout()
    
    plt.show()
    
    fig_height = (10*1)*cm*2
    fig_width = (10*1)*cm*2
    fig, ax = plt.subplots(1,1,figsize=(fig_width,fig_height),sharey=True)
    medianprops = dict(linestyle='-', linewidth=3, color='black')
    meanpointprops = dict(marker='o', markersize=8, markeredgecolor='black',markerfacecolor='firebrick')
    bplot=ax.boxplot(snr_params_all,vert = True,whis=(5,95),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
    ax.set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'])
    ax.set(ylabel='Absolute Parameter Error',xlabel='SNR')
    fig.suptitle('Case {}'.format(p_case[c-1]))
    plt.show()    

Kva_list_all = np.transpose(np.concatenate(Kva_list_all, axis=1))
ua_list_all = np.transpose(np.concatenate(ua_list_all, axis=1))
uv_list_all = np.transpose(np.concatenate(uv_list_all, axis=1))
Jpa_list_all = np.transpose(np.concatenate(Jpa_list_all, axis=1))
Jna_list_all = np.transpose(np.concatenate(Jna_list_all, axis=1))
Ja_list_all = np.transpose(np.concatenate(Ja_list_all, axis=1))

bplot_colors = ["#c7acf9", "#c7acf9", "#c7acf9", "#c7acf9"]
fig_height = (10*1)*cm*2
fig_width = (10*4)*cm*2
row=1
col=4
fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height))
medianprops = dict(linestyle='-', linewidth=3, color='black')
meanpointprops = dict(marker='o', markersize=8, markeredgecolor='black',markerfacecolor='firebrick')
bplot1=ax[0].boxplot(Kva_list_all,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot1['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[0].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$K^{va}$ $(1/s)$')
ax[0].set(ylabel=r'$\% \bar{E}_{rel} $')
ax[0].set(xlabel='SNR')
bplot2=ax[1].boxplot(ua_list_all,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot2['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[1].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$u^{a}$ $(cm/s)$')
ax[1].set(xlabel='SNR')
bplot3=ax[2].boxplot(uv_list_all,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot3['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[2].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$u^{v}$ $(cm/s)$')
ax[2].set(xlabel='SNR')
bplot4=ax[3].boxplot(Ja_list_all,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot4['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[3].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$J^{a}$ $(mM/s)$')
ax[3].set(xlabel='SNR')

fig.savefig('Figures/paper_1d2c_SNR_vel_param_relative_error_purple.eps', dpi=500, bbox_inches='tight')

print('Kva SNR relative error means: {}'.format(np.mean(Kva_list_all,axis=0)))
print('ua SNR relative error means: {}'.format(np.mean(ua_list_all,axis=0)))
print('uv SNR relative error means: {}'.format(np.mean(uv_list_all,axis=0)))
print('Ja SNR relative error means: {},{},{},{}'.format(np.mean(Ja_list_all[:,0]),np.mean(Ja_list_all[:,1]),np.mean(Ja_list_all[:,2]),np.mean(Ja_list_all[:,3])))

print('Kva SNR relative error stdev: {}'.format(np.std(Kva_list_all,axis=0)))
print('ua SNR relative error stdev: {}'.format(np.std(ua_list_all,axis=0)))
print('uv SNR relative error stdev: {}'.format(np.std(uv_list_all,axis=0)))
print('Ja SNR relative error stdev: {},{},{},{}'.format(np.std(Ja_list_all[:,0]),np.std(Ja_list_all[:,1]),np.std(Ja_list_all[:,2]),np.std(Ja_list_all[:,3])))


#%% Dt recovery per case
timepoints = [2,4,6,8]

Kva_all = []
Jpa_all  = []
Jna_all  = []
ua_all  = []
uv_all  = []

for c in cases:
    panel = 0
    k=0
    for Dt in timepoints:
        if Dt == 6:
            opt_vals = np.load("case{}_guess1_0snr_0.00005step_{}Dt_10000maxiter/FinalOptimisationValues_78.0s.npz".format(c,Dt),allow_pickle=True)
            args = np.load("case{}_guess1_0snr_0.00005step_{}Dt_10000maxiter/args_in_opt.npz".format(c,Dt),allow_pickle=True)
        else:
            opt_vals = np.load("case{}_guess1_0snr_0.00005step_{}Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,Dt),allow_pickle=True)
            args = np.load("case{}_guess1_0snr_0.00005step_{}Dt_10000maxiter/args_in_opt.npz".format(c,Dt),allow_pickle=True)
    
        value_tuple = opt_vals['value_tuple']
        x_fine = args['x_fine']
        xc_fine = args['xc_fine']
        tsamp=args['tsamp']
        i = opt_vals['i']
        t=args['t']
        C_tiss = args['C_tiss']
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
        
        fa, fv, Da, Dv, F, va, va_frac, vv, v, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_f(1,2,c)
        ua, uv, Da, Dv, Kva, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_v(1,2,c)
        
        Jpa_inv = Cpa_inv * (fa_inv[0]/(dx))
        Jna_inv = Cna_inv * (abs(fa_inv[-1])/(dx))
        Jpv_inv = Cpv_inv * (abs(fv_inv[0])/(dx))
        Jnv_inv = Cnv_inv * (fv_inv[-1]/(dx))

        C_inv, Ca_inv, Cv_inv = fmod_f(Lx, dx, tsamp[:i], dt, fa_inv, fv_inv, F_inv, Da_inv, Dv_inv, va_inv, vv_inv,  np.append(Jpa_inv,Jpa_inv[-1]), np.append(Jpv_inv,Jpv_inv[-1]), np.append(Jna_inv,Jna_inv[-1]), np.append(Jnv_inv,Jnv_inv[-1]))
    
        if Dt==2:
            plt.rcParams.update({'font.size':34})
            
            cm = 1/2.54  # centimeters in inches
            two_column = 17.56 #cm
            
            fig_height = (10*5)*cm*2
            fig_width = (10*4)*cm*2
            row=5
            col=4
            fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height),sharey='row')
            ax[0,panel].set(xlabel='x (cm)', ylabel='$K^{va}$ (1/s)')
            ax[1,panel].set(xlabel='x (cm)', ylabel='$u^{a}$ (cm/s)')
            ax[2,panel].set(xlabel='x (cm)', ylabel='$u^{v}$ (cm/s)')
            ax[3,panel].set(xlabel='t (s)', ylabel='$J^{a}_{+}$ (mM/s)')
            ax[4,panel].set(xlabel='t (s)', ylabel='$J^{a}_{-}$ (mM/s)')
        ax[0,panel].set(xlabel='x (cm)')
        ax[1,panel].set(xlabel='x (cm)')
        ax[2,panel].set(xlabel='x (cm)')
        ax[3,panel].set(xlabel='t (s)')
        ax[4,panel].set(xlabel='t (s)')

        if c ==11:
            ax[0,panel].plot(xc_fine,Kva,'k')
            ax[0,panel].plot(xc_fine,Kva_inv,colors[k],linestyle='None',marker='o')
            ax[0,panel].set(title='Dt: {}s'.format(Dt))
            ax[1,panel].plot(x_fine,ua,'k')
            ax[1,panel].plot(x_fine,ua_inv,colors[k],linestyle='None',marker='o')
            
            ax[2,panel].plot(x_fine,uv,'k')
            ax[2,panel].plot(x_fine,uv_inv,colors[k],linestyle='None',marker='o')
            
            ax[3,panel].plot(t,Jpa,'k')
            ax[3,panel].plot(tsamp[:i-1],Jpa_inv,colors[k],linestyle='None',marker='o')
            
            ax[4,panel].plot(t,Jna,'k')
            ax[4,panel].plot(tsamp[:i-1], Jna_inv,colors[k],linestyle='None',marker='o')
         
        else:
            ax[0,panel].plot(xc_fine/10,Kva,'k')
            ax[0,panel].plot(xc_fine/10,Kva_inv,colors[k],linestyle='None',marker='o')
            ax[0,panel].set(title='Dt: {}s'.format(Dt))
            ax[1,panel].plot(x_fine/10,ua/10,'k')
            ax[1,panel].plot(x_fine/10,ua_inv/10,colors[k],linestyle='None',marker='o')
            
            ax[2,panel].plot(x_fine/10,uv/10,'k')
            ax[2,panel].plot(x_fine/10,uv_inv/10,colors[k],linestyle='None',marker='o')
            
            ax[3,panel].plot(t,Jpa,'k')
            ax[3,panel].plot(tsamp[:i-1],Jpa_inv,colors[k],linestyle='None',marker='o')
            
            ax[4,panel].plot(t,Jna,'k')
            ax[4,panel].plot(tsamp[:i-1], Jna_inv,colors[k],linestyle='None',marker='o')
        panel += 1
    
    plt.tight_layout()       
    fig.savefig('Figures/paper_1d2c_Dtall_vel_param_case{}.png'.format(c),bbox_inches='tight')

#%% timepoints avg param error
cases = [6,11,3]
timepoints = [2,4,6,8]


Kva_list = []
ua_list = []
uv_list = []
Jpa_list = []
Jna_list = []
Ja_list = []

Kva_all = []
Jpa_all  = []
Jna_all  = []
ua_all  = []
uv_all  = []

for Dt in timepoints:    
    ua_diff_all = []  
    uv_diff_all = []
    Kva_diff_all =  [] 
    Jpa_diff_all =  [] 
    Jna_diff_all =  []
    Ja_diff_all = []
    for c in cases:
        Dt=int(Dt)
        if Dt == 6:
            opt_vals = np.load("case{}_guess1_0snr_0.00005step_{}Dt_10000maxiter/FinalOptimisationValues_78.0s.npz".format(c,Dt),allow_pickle=True)
            args = np.load("case{}_guess1_0snr_0.00005step_{}Dt_10000maxiter/args_in_opt.npz".format(c,Dt),allow_pickle=True)
        else:
            opt_vals = np.load("case{}_guess1_0snr_0.00005step_{}Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,Dt),allow_pickle=True)
            args = np.load("case{}_guess1_0snr_0.00005step_{}Dt_10000maxiter/args_in_opt.npz".format(c,Dt),allow_pickle=True)
    
        value_tuple = opt_vals['value_tuple']
        x_fine = args['x_fine']
        xc_fine = args['xc_fine']
        tsamp=args['tsamp']
        i = opt_vals['i']
        t=args['t']
        C_tiss = args['C_tiss']
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
        
        fa, fv, Da, Dv, F, va, va_frac, vv, v, Jpa, Jpv, Jna, Jnv, dt, Nt, t = load_kinetics_f(1,2,c)
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
        Jpa = Jpa[0::skip]
        Jna = Jna[0::skip]
            
        ua_diff = 100*(abs(ua-ua_inv)/np.mean(abs(ua)))
        uv_diff = 100*(abs(uv-uv_inv)/np.mean(abs(uv)))
        Kva_diff = 100*(abs(Kva-Kva_inv)/np.mean(abs(Kva)))
        Jpa_diff = 100*(abs(Jpa[:-2]-Jpa_inv)/np.mean(abs(Jpa)))
        Jna_diff = 100*(abs(Jna[:-2]-Jna_inv)/np.mean(abs(Jna)))
        
        if c==6 and Dt==2:
            print('Case, Dt, ua mean, ua std, uv mean, uv std, Kva mean, Kva std, Ja mean, Ja std ')
        print('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}'.format(c,Dt,np.mean(ua_diff),np.std(ua_diff),np.mean(uv_diff),np.std(uv_diff),np.mean(Kva_diff),np.std(Kva_diff),np.mean((Jpa_diff,Jna_diff)),np.std((Jpa_diff,Jna_diff))))
    
            
        ua_diff_all += [ua_diff,]   
        uv_diff_all += [uv_diff,]
        Kva_diff_all += [Kva_diff,]
        Ja_diff_all += [Jpa_diff,]
        Ja_diff_all += [Jna_diff,]
    
    ua_diff_all = np.concatenate(ua_diff_all)
    uv_diff_all = np.concatenate(uv_diff_all)
    Kva_diff_all = np.concatenate(Kva_diff_all)
    Ja_diff_all = np.concatenate(Ja_diff_all)
        
    Kva_list += [Kva_diff_all,]
    ua_list += [ua_diff_all,]
    uv_list += [uv_diff_all,]
    Jpa_list += [Jpa_diff_all,]
    Jna_list += [Jna_diff_all,]
    Ja_list += [Ja_diff_all,]
  

bplot_colors = ["#78c99c", "#78c99c", "#78c99c", "#78c99c"]
fig_height = (10*1)*cm*2
fig_width = (10*4)*cm*2
row=1
col=4

fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height))
medianprops = dict(linestyle='-', linewidth=3, color='black')
meanpointprops = dict(marker='o', markersize=8, markeredgecolor='black',markerfacecolor='firebrick')
bplot1=ax[0].boxplot(Kva_list,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot1['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[0].set(xticks=[1,2,3,4], xticklabels=['2','4','6','8'],title='$K^{va}$ $(1/s)$')
ax[0].set(ylabel=r'$\% \bar{E}_{rel} $')
ax[0].set(xlabel='Dt (s)')
bplot2=ax[1].boxplot(ua_list,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot2['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[1].set(xticks=[1,2,3,4], xticklabels=['2','4','6','8'],title='$u^{a}$ $(cm/s)$')
ax[1].set(xlabel='Dt (s)')
bplot3=ax[2].boxplot(uv_list,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot3['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[2].set(xticks=[1,2,3,4], xticklabels=['2','4','6','8'],title='$u^{v}$ $(cm/s)$')
ax[2].set(xlabel='Dt (s)')
bplot4=ax[3].boxplot(Ja_list,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot4['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[3].set(xticks=[1,2,3,4], xticklabels=['2','4','6','8'],title='$J^{a}$ $(mM/s)$')
ax[3].set(xlabel='Dt (s)')



fig.savefig('Figures/paper_1d2c_Dtall_rel_diff_vel_param_allcase_green.eps', dpi=500, bbox_inches='tight')

print('Kva Dt relative error means: {}'.format(np.mean(Kva_list,axis=1)))
print('ua Dt relative error means: {}'.format(np.mean(ua_list,axis=1)))
print('uv Dt relative error means: {}'.format(np.mean(uv_list,axis=1)))
print('Ja Dt relative error means: {},{},{},{}'.format(np.mean(Ja_list[0]),np.mean(Ja_list[1]),np.mean(Ja_list[2]),np.mean(Ja_list[3])))



print('Kva Dt relative error stds: {}'.format(np.std(Kva_list,axis=1)))
print('ua Dt relative error stds: {}'.format(np.std(ua_list,axis=1)))
print('uv Dt relative error stds: {}'.format(np.std(uv_list,axis=1)))
print('Ja Dt relative error stds: {},{},{},{}'.format(np.std(Ja_list[0]),np.std(Ja_list[1]),np.std(Ja_list[2]),np.std(Ja_list[3])))

