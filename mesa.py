
import numpy as np
import os
from astropy.io import ascii

class mesa:

    def init(self):
    ###### parameters ####################
        self.path = '/Users/Jason/code/mesa_progenitors/'  # file path
#        self.type = 'solar'                      # file type
        self.name = 'profileXX.data'                      # file name

    def  __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)

    # sample lab 1 solution
    def read_mesa_profile(self):
        h = ascii.read(self.name,header_start=1,data_start=2,data_end=3)
        f = ascii.read(self.name,header_start=4,data_start=5)
        return h,f

    def read_mesa_star(self):
        h,f = self.read_mesa_profile()
        self.data = f
        self.header = h
