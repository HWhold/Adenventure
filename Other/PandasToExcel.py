# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 09:06:36 2021
This script exports lists of dataframes with named sheets

@author: X202722
"""
with pd.ExcelWriter('output.xlsx') as writer:  
    for i in range(len(extraction_res)):
        data = extraction_res[i]
        data.to_excel(writer, sheet_name = str(Ts[i]))