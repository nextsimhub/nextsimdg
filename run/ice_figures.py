#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create figures with output sea ice model. 

@author: ivo
"""

import netCDF4 as nc
import numpy as np 
from matplotlib import pyplot as plt
import xarray as xr
import os, re
import datetime

#Directory containing this file. 
DIR = os.path.dirname(__file__)


class DiagnosticsDB:
    """
    This class reads values of all diagnostic output files. 
    
    path : str 
        Path 
    
    """
    
    def __init__(self, path):
        """ Constructor """
        self.path = path
        self._get_files()
        self._get_variables()
        
    def _get_files(self):
        """
        Get filenames of all diagnostic files in output directory. 
        
        """
        #Find all diagnostic files. 
        pattern = re.compile("diagnostic.(.*).nc")
        files_times = [(file, re.findall(pattern,file)) for file in os.listdir(self.path)]
        
        #Extract times from file names.
        self.files, self.times = [], []
        for file, time in files_times:
            if len(time)>0:
                self.times.append(datetime.datetime.strptime(time[0],
                                                             '%Y-%m-%dT%H:%M:%SZ'))
                self.files.append(os.path.join(self.path, file))
                
        #Sort by time. 
        sorted_indices = np.argsort(self.times)
        self.times = np.array(self.times, dtype=object)[sorted_indices]
        self.files = np.array(self.files, dtype=object)[sorted_indices]
       
    def _get_variables(self):
        """ Get names of variables in files."""
        self.variables = set([])
        for file in self.files:
            with nc.Dataset(file) as data:
                var_in_data = set([str(v) for v in data.groups['data'].variables])
                self.variables = self.variables.union(var_in_data)
                
    def _get_mask(self):
        if len(self.files)>0 and "mask" in self.variables:
            with nc.Dataset(file) as data:
                self.mask = np.array(data.groups["data"]["mask"], dtype=np.bool)
        else:
            self.mask = None
            
    def read(self, variable):
        """ Read the values for a fields 
        
        Parameters
        ----------
        variable : str 
            Name of variable for which data are requested. 
            
        Returns
        -------
        Requested data. 
        
        """
        
        #Check if variable exist. 
        if variable not in self.variables:
            raise ValueError(f"(variable) is not a valid variable. Valid options are (self.variables).")
        
        values = []
        mask = None
        for file in self.files:
            with nc.Dataset(file) as data:
                values1 = np.array(data.groups["data"][variable])
                
                if mask is None and "mask" in self.variables:
                    mask = np.array(data.groups["data"]["mask"])
                elif mask is None:
                    mask = np.ones_like(values)
                    
                print(file,np.shape(values1),np.shape(mask))
                values1 = np.where(mask, values1, np.nan)
                
                values.append(values1)
                
        return np.array(values)   
    