
import os
import scipy.io
from pathlib import Path

sourcefolder = 'in_ultima_test' #inside data/temp_spyder

rootpath = os.getcwd() #make sure the current folder is the folder containing both the Toolbox and the experiment folder (e.g. 'TORus')
sourcepath = os.path.join(rootpath, 'TORus',  'data', 'temp_spyder',sourcefolder ) #where the dml files are
output_path = os.path.join(rootpath,'TORus', 'data','temp_spyder', sourcefolder)#where to output the mat files


# GET FILES NAMES FROM FOLDER
listdir = os.listdir(Path(sourcepath))
filelist =[]
for dmlfile in listdir:
    if dmlfile.endswith(".rsd"):
        filelist.append (os.path.join(sourcepath, dmlfile))

# load python input file and output path
function =os.path.join(rootpath,'VSDI_ourToolbox', 'functions','external','micam_data_extraction','extract_images.py')
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