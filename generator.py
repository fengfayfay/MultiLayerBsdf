import sys
#sys.path.append('/usr/local/lib/python3.6/site-packages')
#import math
#import matplotlib.pyplot as plt

import io
import os
import glob
from subprocess import call, run
import subprocess
templateFile = 'raytraceHeightField/test_temp.pbrt'
numrays = 30000000

#alphas = [3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
alphas = [1.0, 2.0, 3.0, 4.0, 5.0,  6.0, 7.0, 8.0, 9.0]
#alphas = [1.0, 2.0, 3.0, 4.0]
alphas = [ 6.0, 7.0, 8.0, 4.0, 3.0, 2.0, 1.0]
alphas = [4.0]
for alpha in alphas:
    alpha_str = '{0:1.2f}'.format(alpha)
    print (alpha_str)
    actual_alpha_str = '{0:1.2f}'.format(alpha * 0.1)
    print(actual_alpha_str)
    file_alpha_str = '{0:1.1f}'.format(alpha * 0.1)
    print(file_alpha_str)
    alpha_dir = 'alpha' + actual_alpha_str
    print(alpha_dir)
    run(args = ['mkdir', alpha_dir])
    tempFileName = 'test.'+actual_alpha_str + '.pbrt'
    print(tempFileName)
    text = open(templateFile, 'r').read()
    text = text.replace('alpha', file_alpha_str)
    outfile = open(tempFileName, 'w')
    outfile.write(text)
    outfile.close()
    run(args = ['/Users/fengxie/work/bin/gaussianClean/pbrt', '-alpha', alpha_str, '-numrays', repr(numrays),
        tempFileName])
    tmpfiles = glob.glob('*.p')
    print(tmpfiles)
    for t in tmpfiles:
        run(args = ['mv', t, alpha_dir])
    run(args = ['mv', tempFileName, alpha_dir])
