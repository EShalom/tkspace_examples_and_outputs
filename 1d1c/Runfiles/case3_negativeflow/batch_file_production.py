#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:37:14 2021

@author: eveshalom
"""
import os.path
from os import path

phs = []
flist = []
pylist = []
case = [3]
guess = [1,2,3,4,5,6]

snr = ['0','5','10', '15', '20']
Dt = [2,4,6,8,10]
maxiter = 10000
step = 0.0001

ph_dir = 'placeholder_value_files'
batch_py_dir = 'batch_python_files'

if not path.exists(ph_dir):
    os.makedirs(ph_dir)
if not path.exists(batch_py_dir):
    os.makedirs(batch_py_dir)

for i in range (1,8):
    phs += ['PH{}'.format(i)]
    
with open('{}/placeholders.txt'.format(ph_dir), 'w') as filehandle:
    for listitem in phs:
        filehandle.write('%s\n' % listitem)

for c in case:
    for g in guess:
        for noise in [0]:
            fname = '1d1c_case{}_guess{}_{}snr_0.0001step_{}Dt_{}maxiter'.format(c,g,noise,Dt[0],maxiter)
            pylist += [fname+'.py']
            flist += [fname+'.txt']
            
            vals = []
            vals += [str(g)]
            vals += [c]
            vals += [str(noise)]
            vals += [str(Dt[0])]
            vals += [str(maxiter)]
            vals += [step]
            vals += [fname]
            
            
            with open('{}/{}.txt'.format(ph_dir,fname), 'w') as filehandle:
                for listitem in vals:
                    filehandle.write('%s\n' % listitem)

    for noise in snr[1:]:
        for i in range(1,6):
            fname = '1d1c_case{}_guess{}_{}snr_r{}_0.0001step_{}Dt_{}maxiter'.format(c,guess[3],noise,i,Dt[0],maxiter)
            pylist += [fname+'.py']
            flist += [fname+'.txt']
            
            vals = []
            vals += [str(guess[3])]
            vals += [c]
            vals += [str(noise)]
            vals += [str(Dt[0])]
            vals += [str(maxiter)]
            vals += [step]
            vals += [fname]
            

            with open('{}/{}.txt'.format(ph_dir,fname), 'w') as filehandle:
                for listitem in vals:
                    filehandle.write('%s\n' % listitem)

for time in Dt[1:]:
        fname = '1d1c_case{}_guess{}_{}snr_0.0001step_{}Dt_{}maxiter'.format(c,guess[3],snr[0],time,maxiter)
        pylist += [fname+'.py']
        flist += [fname+'.txt']
        
        vals = []
        vals += [str(guess[3])]
        vals += [c]
        vals += [str(snr[0])]
        vals += [str(time)]
        vals += [str(maxiter)]
        vals += [step]
        vals += [fname]
        
        
        with open('{}/{}.txt'.format(ph_dir,fname), 'w') as filehandle:
            for listitem in vals:
                filehandle.write('%s\n' % listitem)



with open('{}/batch_list.txt'.format(ph_dir), 'w') as filehandle:
    for listitem in flist:
        filehandle.write('%s\n' % listitem)
        
with open('{}/scriptrun_list.txt'.format(batch_py_dir), 'w') as filehandle:
    for listitem in pylist:
        filehandle.write('%s\n' % listitem)
