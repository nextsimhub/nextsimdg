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
            
            
    def read(self, variable):
        """ Read the values for a fields 
        
        Parameters
        ----------
        variable : str 
            Name of variable for which data are requested. 
            
        Returns
        -------
        Requested data as np.array. 
        
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
                    
                values1 = np.where(mask, values1, np.nan)
                
                values.append(values1)
                
        return np.array(values)   
    
    
def row_plot(datas):
    """ 
    Plot different datas on a row. 
    
    Parameter
    ---------
    datas : list of tuples
        Each tuple should contain (name, time, field values)
        
    Returns 
    -------
    (figure, axes) handles
    
    """ 
    plt.close('all')
    fig = plt.figure(figsize=(12,6))
    axes = fig.subplots(1,len(datas)).ravel()
    
    #Indices
    IT, INAME, IDATA = 1, 0, 2
    
    for ax, data in zip(axes, datas):
        x, y = np.meshgrid(np.arange(np.size(data[IDATA],0)),
                           np.arange(np.size(data[IDATA],1)))
        p = ax.contourf(x,y,data[IDATA])
        ax.set_title(data[IT])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        
        plt.colorbar(p, ax=ax, label=data[INAME])
        
    return fig, axes
    
#Connect to output files. 
diag = DiagnosticsDB(DIR)


#Plot last time. 
cice = diag.read('cice')
fig, _ = row_plot([('Concentration',diag.times[0],cice[0]), 
              ("Height [m]",diag.times[-1],cice[-1])])
fig.savefig(os.path.join(DIR,"cice_start_end.png"), format="png", dpi=400)

#Plot first time. 
hice = diag.read('hice')
fig, _ = row_plot([('Concentration',diag.times[0],hice[-1]), 
              ("Height [m]",diag.times[-1],hice[-1])])
fig.savefig(os.path.join(DIR,"hice_start_end.png"), format="png", dpi=400)


    
    