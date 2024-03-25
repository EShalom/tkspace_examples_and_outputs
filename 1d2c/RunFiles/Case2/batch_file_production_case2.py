#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 21:09:34 2023

@author: pyess
"""

import os.path
from os import path

phs = []
flist = []
pylist = []

snr = ['0','5','10','15','20']
Dt = ['2','4','6','8','10']
guess = ['1','2', '3', '4', '5']

ph_dir = 'placeholder_value_files'
batch_py_dir = 'batch_python_files'

if not path.exists(ph_dir):
    os.makedirs(ph_dir)
if not path.exists(batch_py_dir):
    os.makedirs(batch_py_dir)

for i in range (1,5):
       phs += ['PH{}'.format(i)]

with open('{}/placeholders.txt'.format(ph_dir), 'w') as filehandle:
       for listitem in phs:
            filehandle.write('%s\n' % listitem)

for g in guess:
    fname = 'case11_guess{}_0snr_0.00005step_1Dt_10000maxiter'.format(g)
    pylist += [fname+'.py']
    flist += [fname+'.txt']
    
    vals = []
    vals += [str(g)]
    vals += [str(snr[0])]
    vals += [str(Dt[0])]
    vals += [fname]
    
    with open('{}/{}.txt'.format(ph_dir,fname), 'w') as filehandle:
        for listitem in vals:
            filehandle.write('%s\n' % listitem)

for g in guess:
    fname = 'case11_guess{}_0snr_0.00005step_2Dt_10000maxiter'.format(g)
    pylist += [fname+'.py']
    flist += [fname+'.txt']
    
    vals = []
    vals += [str(g)]
    vals += [str(snr[0])]
    vals += [str(Dt[1])]
    vals += [fname]
    
    with open('{}/{}.txt'.format(ph_dir,fname), 'w') as filehandle:
        for listitem in vals:
            filehandle.write('%s\n' % listitem)


for noise in snr[1:]:
    for i in range(1,6):
        fname = 'case11_guess1_r{}_{}snr_0.00005step_1Dt_10000maxiter'.format(i,noise)
        pylist += [fname+'.py']
        flist += [fname+'.txt']

        vals = []
        vals += [str(guess[0])]
        vals += [str(noise)]
        vals += [str(Dt[0])]
        vals += [fname]

        with open('{}/{}.txt'.format(ph_dir,fname), 'w') as filehandle:
            for listitem in vals:
                filehandle.write('%s\n' % listitem)

for noise in snr[1:]:
    for i in range(1,6):
        fname = 'case11_guess1_r{}_{}snr_0.00005step_2Dt_10000maxiter'.format(i,noise)
        pylist += [fname+'.py']
        flist += [fname+'.txt']
        
        vals = []
        vals += [str(guess[0])]
        vals += [str(noise)]
        vals += [str(Dt[1])]
        vals += [fname]
        
        with open('{}/{}.txt'.format(ph_dir,fname), 'w') as filehandle:
            for listitem in vals:
                filehandle.write('%s\n' % listitem)

for samp in Dt:
    fname = 'case11_guess1_0snr_0.00005step_{}Dt_10000maxiter'.format(samp)
    pylist += [fname+'.py']
    flist += [fname+'.txt']

    vals = []
    vals += [str(guess[0])]
    vals += [str(snr[0])]
    vals += [str(samp)]
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

