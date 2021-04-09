# Script to calculate total ion current in each mzML file

from pyopenms import *
from sys import argv
import numpy
import pandas as pd

# argv[1] is the list of mzML files to calc TIC for 
# argv[2] is the name of the norm factor saving list

# parsing list of files supplied
#input = argv[1]
#input_parsed = input.split(" ")

#print(input)
#print(input_parsed)

def calcTIC(exp):
    tic = 0
    for spec in exp:
       if spec.getMSLevel() == 1:
            mz, i = spec.get_peaks()
            tic += sum(i)
    return tic

def load_and_calc(file_name):
    exp =  MSExperiment()
    MzMLFile().load(file_name, exp)
    file_tic = calcTIC(exp)
    exp = None
    return file_tic

all_tics = []
all_files = []

#print(argv[2:])

for unique_file in argv[2:]:
    u_file_tic = load_and_calc(unique_file)
    print(unique_file)
    print(u_file_tic)
    all_tics.append(u_file_tic)
    all_files.append(unique_file)

file_df = pd.DataFrame({'file_name': all_files, 'tic': all_tics})
print(file_df)
file_df.to_csv(argv[1])
