#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 14:33:29 2021
 This is to change file names so we can easily import them using ClearMap
@author: user
"""
import os
  
os.chdir('/data01/szucca/Projects/SexualImprinting/SWISS_2022/SW24Unfam/Auto')
print(os.getcwd())
COUNT = 0
  
# Function to increment count 
# to make the files sorted.
def increment():
    global COUNT
    COUNT = COUNT + 1

dirFiles=os.listdir()
dirFiles.sort() 
  
for f in dirFiles:
    f_name, f_ext = os.path.splitext(f)
    f_name = "Z" + str(COUNT).zfill(4);
    increment()
  
    new_name = '{} {}'.format(f_name, f_ext)
    os.rename(f, new_name)

