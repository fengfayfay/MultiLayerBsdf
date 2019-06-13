import sys
#sys.path.append('/usr/local/lib/python3.6/site-packages')
#import math
#import matplotlib.pyplot as plt

import io
import os
import glob
from subprocess import call, run
import subprocess
import time

#time.sleep(1000)

#templatePrefix = "area2_gs"
templatePrefix = "gold"
templateFile = 'template/' + templatePrefix + '_temp.pbrt'
#alphas = [9.0, 8.0, 7.0, 5.0, 4.0, 3.0, 2.0, 1.0]
USEMSs = ["true", "false"]
NVPs = ["true", "false"]
MSs = ["false", "gs", "nvp"]
hidden_nvp = 30 
hidden_fresnel = 10
NVPMODELPREFIX = '/Users/fengxie/tensor_results_ms_reflection/bounce_all'
FRESNELMODELPREFIX = '/Users/fengxie/tensor_results_ms_reflection/fresnel'
DEPTH = 3
SAMPLES = 64 


alphas = [9.0, 7.0, 5.0]
#alphas = [9.0]
rotations = [0]
for alpha in alphas:
    alpha_str = '{0:1.2f}'.format(alpha)
    print (alpha_str)
    actual_alpha_str = '{0:1.1f}'.format(alpha * 0.1)
    print(actual_alpha_str)
    file_alpha_str = '{0:1.1f}'.format(alpha * 0.1)
    print(file_alpha_str)
    alpha_dir = 'nvp_' + repr(hidden_nvp) + '_bounce_' + repr(DEPTH) + '_sample_' + repr(SAMPLES) +  '/alpha' + actual_alpha_str 
    print(alpha_dir)
    os.makedirs(alpha_dir, exist_ok = True)
    for rotation in rotations:
        for ms in MSs:
            text = open(templateFile, 'r').read()
            text = text.replace('DEPTH', repr(DEPTH))
            text = text.replace('SAMPLES', repr(SAMPLES))
            text = text.replace('HIDDENNVP', repr(hidden_nvp))
            text = text.replace('HIDDENFRESNEL', repr(hidden_fresnel))
            text = text.replace('ALPHA', file_alpha_str)
            text = text.replace('FRESNELMODELPREFIX', FRESNELMODELPREFIX)
            text = text.replace('NVPMODELPREFIX', NVPMODELPREFIX)
            
            threadCount = 2
            USENVP = "false"
            if ms == "false":
                USEMS = "false"
            else:
                USEMS = "true"
                if ms == "nvp":
                    threadCount = 1
                    USENVP = "true"
                
            text = text.replace('USENVP', USENVP)
            text = text.replace('USEMS', USEMS)
            text = text.replace('MS', ms)
 
            tempFileName = templatePrefix + '_ms_' + ms + '_a_' + actual_alpha_str +  '.pbrt'
            outfile = open(tempFileName, 'w')
            outfile.write(text)
            outfile.close()
            print(tempFileName)
            run(args = ['/Users/fengxie/work/bin/GaussianBSDF/pbrt', '--nthreads', repr(threadCount),  tempFileName])
            tmpfiles = glob.glob('*.png')
            #print(tmpfiles)
            for t in tmpfiles:
                run(args = ['mv', t, alpha_dir])
            run(args = ['mv', tempFileName, alpha_dir])
