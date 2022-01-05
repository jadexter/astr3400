
import numpy as np
import os

# global variables about where different quantities are in data:

RAD=1
PRES=5
DENS=3
MASS=0
EINT=6
VEL=2
ENTR=7
XH=12
XHE3=13
XHE4=14
XC=15
XO=17

class heger:

    def init(self):
    ###### parameters ####################
        self.path = '/Users/Jason/code/heger_progenitors/'  # file path
        self.type = 'solar'                      # file type
        self.name = 's25.0'                      # file name

    def  __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)

    def read_heger_star(self):
        os.system('gunzip '+self.path+self.type+'/'+self.name+'.gz')
        indx=np.concatenate((np.arange(11)+1,14+np.arange(20)))
        file = open(self.path+self.type+'/'+self.name)
        # read header
        i=0
        header1=file.readline()
        header2=file.readline().split()
        print('header: ',header2)
        # read data
        for line in file:
            ln=np.asarray(line.split())
            ln=ln[indx]
            ln=np.where(ln!='---',ln,'0')
#            print 'i: ', i,len(ln)
#            x=np.array(map(float,ln))
            x=np.array(ln,dtype=float)
            if i==0: data=x
            else: data=np.concatenate((data,x))
            i=i+1

        data=np.transpose(data.reshape((i,len(indx))))
        #print np.shape(data), len(data)

        # store in class:
        self.header=(np.asarray(header2))[indx]
        self.data=data

        os.system('gzip '+self.path+self.type+'/'+self.name)

    def interpolate_heger_star(self,r):
        # interpolate data to locations r
        # for now just piecewise linear from numpy but should try averaging over relevant regions of star too
        pint=np.interp(r,self.data[RAD,:],self.data[PRES,:])
        dint=np.interp(r,self.data[RAD,:],self.data[DENS,:])
        vint=np.interp(r,self.data[RAD,:],self.data[VEL,:])
        mint=np.interp(r,self.data[RAD,:],self.data[MASS,:])
        eint=np.interp(r,self.data[RAD,:],self.data[EINT,:])
        xhint=np.interp(r,self.data[RAD,:],self.data[XH,:])
        xhe4int=np.interp(r,self.data[RAD,:],self.data[XHE4,:])
        xcint=np.interp(r,self.data[RAD,:],self.data[XC,:])
        xoint=np.interp(r,self.data[RAD,:],self.data[XO,:])
        return pint,dint,vint,mint,eint,xhint,xhe4int,xcint,xoint

    def identify_hhe(self):
# find H/He interface location
        a

