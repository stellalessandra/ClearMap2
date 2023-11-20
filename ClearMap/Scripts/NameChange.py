#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 14:33:29 2021
 This is to change file names so we can easily import them using ClearMap
@author: user
"""
import os
  
os.chdir('/media/astella/Elements/Bovetti/WILD/20230621_W20_cfos/20230629_122136/t000/c000/RES(8486x3850x1040)/00-163/00-163_000000/')
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

