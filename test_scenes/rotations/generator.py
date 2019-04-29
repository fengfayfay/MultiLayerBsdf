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
templatePrefix = "area2"
templateFile = 'template/' + templatePrefix + '_temp.pbrt'
#alphas = [9.0, 8.0, 7.0, 5.0, 4.0, 3.0, 2.0, 1.0]
MSs = ["true", "false"]
NVPs = ["true", "false"]
light_intensity = 50

hidden = 30
rotation_step = 20


alphas = [9.0, 8.0, 7.0, 6.0, 5.0]
alphas = [4.0, 3.0, 2.0, 1.0]
for alpha in alphas:
    alpha_str = '{0:1.2f}'.format(alpha)
    print (alpha_str)
    actual_alpha_str = '{0:1.1f}'.format(alpha * 0.1)
    print(actual_alpha_str)
    file_alpha_str = '{0:1.1f}'.format(alpha * 0.1)
    print(file_alpha_str)
    alpha_dir = 'nvp_' + repr(hidden) + '_light_' + repr(light_intensity) + '/alpha' + actual_alpha_str 
    print(alpha_dir)
    os.makedirs(alpha_dir, exist_ok = True)
    for rotation in range(10, 90, rotation_step):
        for ms in MSs:
            text = open(templateFile, 'r').read()
            text = text.replace('HIDDEN', repr(hidden))
            text = text.replace('ROTATION', repr(rotation))
            text = text.replace('ALPHA', file_alpha_str)
            text = text.replace('MS', ms)
            text = text.replace('LIGHTINTENSITY', repr(light_intensity))
            
            tempFileName = templatePrefix + '_ms_' + ms + '_a_' + actual_alpha_str + '_r_' + repr(rotation) + '.pbrt'
            #print(tempFileName)
            outfile = open(tempFileName, 'w')
            outfile.write(text)
            outfile.close()
            threadCount = 1 if ms == "true" else 2
            run(args = ['/Users/fengxie/work/bin/GaussianBSDF/pbrt', '--nthreads', repr(threadCount),  tempFileName])
            tmpfiles = glob.glob('*.png')
            #print(tmpfiles)
            for t in tmpfiles:
                run(args = ['mv', t, alpha_dir])
            run(args = ['mv', tempFileName, alpha_dir])
