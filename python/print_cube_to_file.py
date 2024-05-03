#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 11:24:53 2024

@author: astahl3
"""
'''
output_file = '/data1/astahl3/aikef_output/cube_runs/G8/aikef_fields/lzfile_temp.txt'
start_value = -87.3214
max_value = 87.5467+0.00001

numNds_X = 672
numNds_Y = 672
numNds_Z = 672
dz = 0.260608197
lines_written = 0

with open(output_file, 'w') as outfile:
    for i in range(1,numNds_X+1):
        for j in range(1,numNds_Y+1):

            current_value = start_value
            while current_value < max_value:                            
                print_value = round(current_value,4)
                outfile.write(f'{print_value}\n')
                current_value += dz
                lines_written += 1
'''
#print(f'Total lines written: {lines_written}')
            
            
test = range(0,10)
print(test[len(test)-1])
            