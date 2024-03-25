#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 17:23:30 2022

@author: pyess
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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
#with open('batch_list.txt') as file:
#    lines = [line[:-5].rstrip() for line in file]
cases = [1,2,3]
guess = [1,2,3,4,5,6]
colors = sns.color_palette("tab10", n_colors=12).as_hex()

colors = colors[4:10]
k=0

#%% Different guess case plot
plt.rcParams.update({'font.size':34})

cm = 1/2.54  # centimeters in inches
two_column = 17.56 #cm

fig_height = (10*5)*cm
fig_width = (10*4)*cm*2
row=2
col=3
fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height),sharey='row')
ax[0,0].set(xlabel='x (cm)', ylabel='$u$ (cm/s)')

m=-1
for c in cases:  
    k = 0
    m=m+1    
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
        
        
        
        
        if g == 1:                    
            ax[0,m].set(xlim=([0,Lx/10]),xlabel='x (cm)')
            ax[1,m].set(xlim=([0,T]),xlabel='t (s)')
            
            ax[0,m].plot(x_fine/10,u/10,'k')
            if c == 1 or c== 2:
                ax[1,m].set(xlabel='t (s)', ylabel='$J_{+}$ (mM/s)',label=g)
                ax[1,m].plot(t,Jp,'k')
            else:
                ax[1,m].set(xlabel='t (s)', ylabel='$J_{-}$ (mM/s)',label=g)
                ax[1,m].plot(t,Jn,'k')
                
            ax[1,m].set(ylim=([0,7]))
        
        ax[0,m].plot(x_fine/10,u_inv/10,colors[k], marker='o', linestyle='',label=g)
        if c == 1 or c== 2:
            ax[1,m].plot(tsamp[:i-1], Jp_inv,colors[k], marker='o', linestyle='')
        else:
            ax[1,m].plot(tsamp[:i-1], Jn_inv,colors[k], marker='o', linestyle='')
        
        ax[0,m].set(title='Case {}'.format(m+1))
        k+=1
    legend = ax[0,2].legend(title="Guess",
                    loc='upper right', fontsize=24, fancybox=True, ncol=2)
        
fig.savefig('Figures/1d1c_all_guesses_vel_param_cases_new.png')


#%% snr plots
snr = [5,10,15,20]
repeats= [1,2,3,4,5]

Jp_all  = []
Jn_all  = []
u_all  = []
all_avg_pdiff = []
cost_all = []

for c in cases:  
    panel = 0
    for s in snr:
        k=0   
        Jp_all  = []
        Jn_all  = []
        u_all  = []
        for r in repeats:
            if c == 1 or c == 2:
                opt_vals = np.load("1d1c_case{}_guess1_{}snr_r{}_0.0001step_2Dt_10000maxiter/FinalOptimisationValues_{}.0s.npz".format(c,s,r,ouput_time),allow_pickle=True)
                args = np.load("1d1c_case{}_guess1_{}snr_r{}_0.0001step_2Dt_10000maxiter/args_in_opt.npz".format(c,s,r),allow_pickle=True)
            if c == 3:
                opt_vals = np.load("1d1c_case{}_guess4_{}snr_r{}_0.0001step_2Dt_10000maxiter/FinalOptimisationValues_{}.0s.npz".format(c,s,r,ouput_time),allow_pickle=True)
                args = np.load("1d1c_case{}_guess4_{}snr_r{}_0.0001step_2Dt_10000maxiter/args_in_opt.npz".format(c,s,r),allow_pickle=True)
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
            xc_res = np.linspace(dx/2,Lx-dx/2,int(Lx/dx))
            
            f_inv, D_inv, v_inv, Cp_inv, Cn_inv = value_tuple
            tissue_params = convertTtoL_oneComp(value_tuple)
            u_inv, D_inv, Cp_inv, Cn_inv = tissue_params
            
            u, D, Jp, Jn, dt, Nt, t = load_kinetics_v(1,1,c)
            
            Jp_inv = Cp_inv * (f_inv/(dx))
            Jn_inv = Cn_inv * (abs(f_inv)/(dx))
            
            
    
            Dt = tsamp[1]-tsamp[0]

            if r ==1 and s == 5:
                plt.rcParams.update({'font.size':34})

                cm = 1/2.54  # centimeters in inches
                two_column = 17.56 #cm

                fig_height = (10*2)*cm*2
                fig_width = (10*4)*cm*2
                row=2
                col=4
                fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height),sharey='row')
                ax[0,panel].set(xlabel='x (cm)', ylabel='$u$ (cm/s)')
                if c==1 or c==2:
                    ax[1,panel].set(xlabel='t (s)', ylabel='$J_{+}$ (mM/s)')
                else:
                    ax[1,panel].set(xlabel='t (s)', ylabel='$J_{-}$ (mM/s)')

            
            ax[0,panel].set(xlabel='x (cm)')
            ax[1,panel].set(xlabel='t (s)')

            ax[0,panel].plot(x_fine/10,u/10,'k')
            ax[0,panel].plot(x_fine/10,u_inv/10,colors[k],linestyle='None',marker='o',label=r)
            ax[0,panel].set(title='SNR: {}'.format(s))
            
            if c ==1 or c==2:
                ax[1,panel].plot(t,Jp,'k')
                ax[1,panel].plot(tsamp[:i-1],Jp_inv,colors[k],linestyle='None',marker='o',label=r)

            else:
                ax[1,panel].plot(t,Jn,'k')
                ax[1,panel].plot(tsamp[:i-1], Jn_inv,colors[k],linestyle='None',marker='o',label=r)

            ax[1,panel].set(ylim=([0,7]))
            
            	
            u_all += [u_inv/10,]
            Jp_all += [Jp_inv,]
            Jn_all += [Jn_inv,]
            k+=1
        u_all = np.asarray(u_all)
        Jp_all  = np.asarray(Jp_all)
        Jn_all  = np.asarray(Jn_all)
        
        u_m = np.mean(u_all,axis=0)
        Jp_m  = np.mean(Jp_all,axis=0)
        Jn_m  = np.mean(Jn_all,axis=0)
        
        u_ci = (1.96*np.std(u_all,axis=0))/np.sqrt(5)
        Jp_ci  = (1.96*np.std(Jp_all,axis=0))/np.sqrt(5)
        Jn_ci  = (1.96*np.std(Jn_all,axis=0))/np.sqrt(5)
    
        ax[0,panel].fill_between(x_fine/10,u_m-u_ci,u_m+u_ci)
        if c ==1 or c==2:
            ax[1,panel].fill_between(tsamp[:i-1],Jp_m-Jp_ci,Jp_m+Jp_ci)
        else:
            ax[1,panel].fill_between(tsamp[:i-1],Jn_m-Jn_ci,Jn_m+Jn_ci)
        legend = ax[1,3].legend(title="Repeat",
                    loc='upper right', fontsize=24, fancybox=True, ncol=1)

        panel += 1
        
    plt.tight_layout()
    fig.savefig('Figures/1d1c_SNR_vel_param_case{}.png'.format(c),bbox_inches='tight')

#%% snr plots relative error plots
cases = [1,2,3]
snr = [5,10,15,20]
repeats= [1,2,3,4,5]

Jp_all  = []
Jn_all  = []
u_all  = []
all_avg_pdiff = []
cost_all = []

u_list_all = []
J_list_all = []

for s in snr:  
    panel = 0
    u_list = []
    J_list = []
    for c in cases:
        k=0   
        Jp_all  = []
        Jn_all  = []
        u_all  = []
        
        u_diff_all =  []  
        J_diff_all =  [] 
        for r in repeats:
            if c == 1 or c == 2:
                opt_vals = np.load("1d1c_case{}_guess1_{}snr_r{}_0.0001step_2Dt_10000maxiter/FinalOptimisationValues_{}.0s.npz".format(c,s,r,ouput_time),allow_pickle=True)
                args = np.load("1d1c_case{}_guess1_{}snr_r{}_0.0001step_2Dt_10000maxiter/args_in_opt.npz".format(c,s,r),allow_pickle=True)
            if c == 3:
                opt_vals = np.load("1d1c_case{}_guess4_{}snr_r{}_0.0001step_2Dt_10000maxiter/FinalOptimisationValues_{}.0s.npz".format(c,s,r,ouput_time),allow_pickle=True)
                args = np.load("1d1c_case{}_guess4_{}snr_r{}_0.0001step_2Dt_10000maxiter/args_in_opt.npz".format(c,s,r),allow_pickle=True)
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
            xc_res = np.linspace(dx/2,Lx-dx/2,int(Lx/dx))
            
            f_inv, D_inv, v_inv, Cp_inv, Cn_inv = value_tuple
            tissue_params = convertTtoL_oneComp(value_tuple)
            u_inv, D_inv, Cp_inv, Cn_inv = tissue_params
            
            u, D, Jp, Jn, dt, Nt, t = load_kinetics_v(1,1,c)
            
            Jp_inv = Cp_inv * (f_inv/(dx))
            Jn_inv = Cn_inv * (abs(f_inv)/(dx))
            
            
    
            Dt = tsamp[1]-tsamp[0]
            u_inv = u_inv/10
            u=u/10
            u_diff = 100*(abs(u-u_inv)/np.mean(abs(u)))
            if c ==1 or c==2:
                J_diff = 100*(abs(Jp[0:-9:10]-Jp_inv)/np.mean(abs(Jp[0:-9:10])))
            else:
                J_diff = 100*(abs(Jn[0:-9:10]-Jn_inv)/np.mean(abs(Jn[0:-9:10])))
            
            u_diff_all += [abs(u_diff),]
            J_diff_all += [abs(J_diff),]
        
        print('u SNR {} case {} relative error means: {}, std: {}'.format(s,c,np.mean(u_diff_all),np.std(u_diff_all)))
        print('J SNR {} case {} relative error means: {}, std: {}'.format(s,c,np.mean(J_diff_all),np.std(J_diff_all)))

        
        u_diff_all = np.concatenate(u_diff_all)
        J_diff_all = np.concatenate(J_diff_all)
        
        u_list += [u_diff_all,]
        J_list += [J_diff_all,]
    
    u_list_all += [np.concatenate(u_list),]
    J_list_all += [np.concatenate(J_list),]


bplot_colors = ["#c7acf9", "#c7acf9", "#c7acf9", "#c7acf9"]
plt.rcParams.update({'font.size':34})
cm = 1/2.54  # centimeters in inches
two_column = 17.56 #cm
fig_height = (10*1)*cm*2
fig_width = (10*2)*cm*2
row=1
col=2
fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height))

medianprops = dict(linestyle='-', linewidth=3, color='black')
meanpointprops = dict(marker='o', markersize=8, markeredgecolor='black',markerfacecolor='firebrick')
bplot1=ax[0].boxplot(u_list_all,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot1['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[0].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$u$ $(cm/s)$')
ax[0].set(ylabel=r'$\% \bar{E}_{rel} $')
ax[0].set(xlabel='SNR')
bplot2=ax[1].boxplot(J_list_all,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot2['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[1].set(xticks=[1,2,3,4], xticklabels=['5','10','15','20'],title='$J$ $(mM/s)$')
ax[1].set(xlabel='SNR')

fig.savefig('Figures/1d1c_SNR_vel_param_relative_error_purple.eps', dpi=500, bbox_inches='tight')

print('u snr relative error means: {} std: {}'.format(np.mean(u_list_all,axis=1),np.std(u_list_all,axis=1)))
print('u snr relative error means: {} std: {}'.format(np.mean(J_list_all,axis=1),np.std(J_list_all,axis=1)))

#%%
cases = [1,2,3]
guess = [1,2,3,4,5]
colors = sns.color_palette("tab10", n_colors=12).as_hex()

colors = colors[4:10]
k=0

timepoints = [2,4,6,8]

Kva_all = []
Jpa_all  = []
Jna_all  = []
ua_all  = []
uv_all  = []

for c in cases:
    panel = 0
    for Dt in timepoints:
        k=0
        for g in [1]:
            if c ==3:
                if Dt == 6 :
                    opt_vals = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/FinalOptimisationValues_78.0s.npz".format(c,3+g,Dt),allow_pickle=True)
                    args = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/args_in_opt.npz".format(c,3+g,Dt),allow_pickle=True)
                else:
                    opt_vals = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,3+g,Dt),allow_pickle=True)
                    args = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/args_in_opt.npz".format(c,3+g,Dt),allow_pickle=True)
            if c == 1 or c ==2:
                if Dt == 6:
                    opt_vals = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/FinalOptimisationValues_78.0s.npz".format(c,g,Dt),allow_pickle=True)
                    args = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/args_in_opt.npz".format(c,g,Dt),allow_pickle=True)
                else:
                    opt_vals = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,g,Dt),allow_pickle=True)
                    args = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/args_in_opt.npz".format(c,g,Dt),allow_pickle=True)
        
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
            
            
            
            if Dt==2 and g ==1:
                plt.rcParams.update({'font.size':34})
                
                cm = 1/2.54  # centimeters in inches
                two_column = 17.56 #cm
                
                fig_height = (10*2)*cm*2
                fig_width = (10*4)*cm*2
                row=2
                col=4
                fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height),sharey='row')
                ax[0,panel].set(xlabel='x (cm)', ylabel='$u$ (cm/s)')
                if c ==1 or c==2:
                    ax[1,panel].set(xlabel='t (s)', ylabel='$J_{+}$ (mM/s)')
                else:
                    ax[1,panel].set(xlabel='t (s)', ylabel='$J_{-}$ (mM/s)')

            ax[0,panel].set(title='Dt: {}s'.format(Dt))
            ax[0,panel].set(xlabel='x (cm)')
            ax[1,panel].set(xlabel='t (s)')


            ax[0,panel].plot(x_fine/10,u/10,'k')
            ax[0,panel].plot(x_fine/10,u_inv/10,colors[-(k+1)],linestyle='None',marker='o')

            if c == 1 or c == 2:
                ax[1,panel].plot(t,Jp,'k')
                ax[1,panel].plot(tsamp[:i-1],Jp_inv,colors[-(k+1)],linestyle='None',marker='o')
            else:
                ax[1,panel].plot(t,Jn,'k')
                ax[1,panel].plot(tsamp[:i-1], Jn_inv,colors[-(k+1)],linestyle='None',marker='o')
                k +=1
panel += 1
    
plt.tight_layout()
fig.savefig('Figures/Dtall_1d1c_vel_param_case{}.png'.format(c),bbox_inches='tight')
fig.savefig('Figures/Dtall_1d1c_vel_param_case{}.tiff'.format(c),bbox_inches='tight')

#%%

cases = [1,2,3]
guess = [1,2,3,4,5]
colors = sns.color_palette("tab10", n_colors=12).as_hex()

colors = colors[4:10]
k=0

timepoints = [2,4,6,8]

Kva_all = []
Jpa_all  = []
Jna_all  = []
ua_all  = []
uv_all  = []

u_list = []
J_list = []

for Dt in timepoints:
    u_diff_all =  []
    J_diff_all =  []
    for c in cases:
        k=0
        Dt = int(Dt)
        for g in [1]:
            if c ==3:
                if Dt == 6 :
                    opt_vals = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/FinalOptimisationValues_78.0s.npz".format(c,3+g,Dt),allow_pickle=True)
                    args = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/args_in_opt.npz".format(c,3+g,Dt),allow_pickle=True)
                else:
                    opt_vals = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,3+g,Dt),allow_pickle=True)
                    args = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/args_in_opt.npz".format(c,3+g,Dt),allow_pickle=True)
            if c == 1 or c ==2:
                if Dt == 6:
                    opt_vals = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/FinalOptimisationValues_78.0s.npz".format(c,g,Dt),allow_pickle=True)
                    args = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/args_in_opt.npz".format(c,g,Dt),allow_pickle=True)
                else:
                    opt_vals = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/FinalOptimisationValues_80.0s.npz".format(c,g,Dt),allow_pickle=True)
                    args = np.load("1d1c_case{}_guess{}_0snr_0.0001step_{}Dt_10000maxiter/args_in_opt.npz".format(c,g,Dt),allow_pickle=True)
        
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
            
            
            
            Dt = tsamp[1]-tsamp[0]
            skip = int(Dt/0.2)
            
            u_inv = u_inv/10
            u=u/10
            u_diff = 100*(abs(u-u_inv)/np.mean(abs(u)))
            if c ==1 or c==2:
                J_diff = 100*(abs(Jp[0:-skip+1:skip]-Jp_inv)/np.mean(abs(Jp[0:-skip+1:skip])))
            else:
                J_diff = 100*(abs(Jn[0:-skip+1:skip]-Jn_inv)/np.mean(abs(Jn[0:-skip+1:skip])))

        u_diff_all += [abs(u_diff),]
        J_diff_all += [abs(J_diff),]
        
        print('u Dt {} case {} relative error means: {}, std: {}'.format(Dt,c,np.mean(u_diff_all),np.std(u_diff_all)))
        print('J Dt {} case {} relative error means: {}, std: {}'.format(Dt,c,np.mean(J_diff_all),np.std(J_diff_all)))


    u_diff_all = np.concatenate(u_diff_all)
    J_diff_all = np.concatenate(J_diff_all)
    
    u_list += [u_diff_all,]
    J_list += [J_diff_all,]

plt.rcParams.update({'font.size':34})
cm = 1/2.54  # centimeters in inches
two_column = 17.56 #cm
bplot_colors = ["#78c99c", "#78c99c", "#78c99c", "#78c99c"]
fig_height = (10*1)*cm*2
fig_width = (10*2)*cm*2
row=1
col=2
fig, ax = plt.subplots(row,col,figsize=(fig_width,fig_height))

medianprops = dict(linestyle='-', linewidth=3, color='black')
meanpointprops = dict(marker='o', markersize=8, markeredgecolor='black',markerfacecolor='firebrick')
bplot1=ax[0].boxplot(u_list,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot1['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[0].set(xticks=[1,2,3,4], xticklabels=['2','4','6','8'],title='$u$ $(cm/s)$')
ax[0].set(ylabel=r'$\% \bar{E}_{rel} $')
ax[0].set(xlabel='Dt (s)')
bplot2=ax[1].boxplot(J_list,vert = True,whis=(10,90),showfliers=False,showmeans=True,meanprops=meanpointprops, medianprops=medianprops, patch_artist=True)
for patch, color in zip(bplot2['boxes'], bplot_colors):
    patch.set_facecolor(color)
ax[1].set(xticks=[1,2,3,4], xticklabels=['2','4','6','8'],title='$J$ $(mM/s)$')
ax[1].set(xlabel='Dt (s)')

fig.savefig('Figures/1d1c_Dt_vel_param_relative_error_green.eps', dpi=500, bbox_inches='tight')

print('u Dt relative error means: {}'.format(np.mean(u_list,axis=1)))
print('J Dt relative error means: {},{},{},{}'.format(np.mean(J_list[0]),np.mean(J_list[1]),np.mean(J_list[2]),np.mean(J_list[3])))
print('u Dt relative error std: {}'.format(np.std(u_list,axis=1)))
print('J Dt relative error std: {},{},{},{}'.format(np.std(J_list[0]),np.std(J_list[1]),np.std(J_list[2]),np.std(J_list[3])))


