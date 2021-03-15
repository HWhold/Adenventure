# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 08:00:35 2021

@author: X202722
"""


from spyder_kernels.utils.iofuncs import load_dictionary 

filename = 'C:/Users/X202722/Desktop/Prozesssynthese/95_GitEntwicklung/data.spydata'
data_dict = load_dictionary(filename)


import pickle
import tarfile

tar = tarfile.open(filename, "r")
# extract all pickled files to the current working directory
tar.extractall()
extracted_files = tar.getnames()
for f in extracted_files:
    if f.endswith('.pickle'):
         with open(f, 'rb') as fdesc:
             data = pickle.loads(fdesc.read())