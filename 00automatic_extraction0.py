# AUTOMATIC EXTRACTION USING PREDRAG'S FUNCTION
"""
Created on Sun Dec 20 23:27:20 2020

To make sure it works: 
    - store this script into the experiment folder
    - in a 'BVdml' subfolder, make a 'in' and 'out' subfolders, with the name of the reference; e.g.: 
        - '200616in' : with all the dml to convert into mat-files
        - '200616out': empty folder where the output mat-files will be stored
    - In the script: change the 'reference' parameter to match the folders from previous point
    - Before running the script, make sure that the variable space is empty and that the 'Files' window
    is in the experiment root folder, and not in any subfolder**
    - Run each part of the script in one time
    
    ** Mind that after running the code, the Files-window will open the ouput folder 
    If problem persists: make sure that the function that the code calls is stored where the code looks for it

"""

import os
import scipy.io
from pathlib import Path

rootpath = os.getcwd() #make sure the current folder is the root folder
sourcepath = os.path.join(rootpath, 'data', 'temp_spyder', 'in221020') #where the dml files are
output_path = os.path.join(rootpath, 'data', 'temp_spyder', 'out221020')#where to output the mat files

# GET FILES NAMES FROM FOLDER
listdir = os.listdir(Path(sourcepath))
filelist =[]
for dmlfile in listdir:
    if dmlfile.endswith(".dml"):
        filelist.append (os.path.join(sourcepath, dmlfile))

# load python input file and output path
function =os.path.join(rootpath,'extract_images3.py')
salida=  output_path


# RUN THIS IN A SECOND STEP ..........................................
for file in filelist:    # GET NAME OF FILE
    filename = file  # [-14:-4]+'.dml'
    savename = file[-14:-4]
    origen = os.path.join(sourcepath, filename)
    # substitute separators for the function
    origen = origen.replace(os.sep, '/')
    salida = salida.replace(os.sep, '/')
    variables = origen + ' ' + salida

    # Apply function to extract array
    runfile(function, args=variables)

    # retry images from dictionary
    data_act = original_data['images']

    # Export to matfile
    savename = savename + '.mat'
    os.chdir(salida)
    scipy.io.savemat(savename, dict(data_act=data_act))
