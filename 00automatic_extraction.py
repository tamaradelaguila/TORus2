# AUTOMATIC EXTRACTION USING PREDRAG'S FUNCTION
"""
Created on Sun Dec 20 23:27:20 2020

@author: Tamara del Águila
"""

import os
import scipy.io
from pathlib import Path

rootpath = os.getcwd() #make sure the current folder is the root folder
sourcepath = os.path.join(rootpath, 'BVdml') #where the dml files are
output_path = os.path.join(rootpath, 'BVdml','output')#where to output the mat files

# GET FILES NAMES FROM FOLDER
listdir = os.listdir(Path(sourcepath))
filelist =[]
for dmlfile in listdir:
    if dmlfile.endswith(".dml"):
        filelist.append (os.path.join(sourcepath, dmlfile))
    
# load python input file and output path
function =os.path.join(rootpath,'functions','external','python2','extract_images3.py')
salida=  output_path


for file in filelist:
    # GET NAME OF FILE
    filename =  file #[-14:-4]+'.dml'
    savename = file[-14:-4]
    origen=  os.path.join(sourcepath,filename)
    # substitute separators for the function
    origen =  origen.replace(os.sep, '/')
    salida = salida.replace(os.sep, '/')
    variables = origen + ' ' + salida
    
    # Apply function to extract array
    runfile(function, args=variables)

    # retry images from dictionary
    data_act= original_data['images']

    # Export to matfile
    savename = savename + '.mat'
    os.chdir(salida)
    scipy.io.savemat(savename, dict(data_act= data_act))