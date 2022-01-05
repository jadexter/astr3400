################################################
# simple lagrangian staggered-mesh hydro code 
# Use it to blow up stars
################################################

import numpy as np
import heger
import constants as pc
import matplotlib.pyplot as plt
import copy
import math_ex
plt.ion()

###############################
# define a zone class

class zone:
    d      = np.empty(1)  # density
    p      = np.empty(1)  # pressure
    e      = np.empty(1)  # internal energy
    q      = np.empty(1)  # artificial viscosity
    v      = np.empty(1)  # velocity of upper interface
    r      = np.empty(1)  # radius of upper interface
    dr     = np.empty(1)  # radial size of zone
    vol    = np.empty(1)  # volume of zone
    mass   = np.empty(1)  # mass of zone
    cs     = np.empty(1)  # sound speed
    mcum   = np.empty(1)  # mass coordinate of zone

    def get_zone_values(self, i):
        new=zone()
        new.d=self.d[i]
        new.p=self.p[i]
        new.e=self.e[i]
        new.q=self.q[i]
        new.v=self.v[i]
        new.r=self.r[i]
        new.dr=self.dr[i]
        new.vol=self.vol[i]
        new.mass=self.mass[i]
        new.cs=self.cs[i]
        new.mcum=self.mcum[i]
        return new

    def add_zone_values(self, new):
        self.d=np.append(self.d,new.d)
        self.p=np.append(self.p,new.p)
        self.e=np.append(self.e,new.e)
        self.q=np.append(self.q,new.q)
        self.v=np.append(self.v,new.v)
        self.r=np.append(self.r,new.r)
        self.dr=np.append(self.dr,new.dr)
        self.vol=np.append(self.vol,new.vol)
        self.mass=np.append(self.mass,new.mass)
        self.cs=np.append(self.cs,new.cs)
        self.mcum=np.append(self.mcum,new.mcum)

    def del_zone_values(self, i):
        self.d = np.delete(self.d, i)
        self.p = np.delete(self.p, i)
        self.e = np.delete(self.e, i)
        self.q = np.delete(self.q, i)
        self.v = np.delete(self.v, i)
        self.r = np.delete(self.r, i)
        self.dr = np.delete(self.dr, i)
        self.vol = np.delete(self.vol, i)
        self.mass = np.delete(self.mass, i)
        self.cs = np.delete(self.cs, i)
        self.mcum = np.delete(self.mcum, i)

    def assign_zone(self, indx, zones):
        self.d[indx] = zones.d
        self.p[indx] = zones.p
        self.e[indx] = zones.e
        self.q[indx] = zones.q
        self.v[indx] = zones.v
        self.r[indx] = zones.r
        self.dr[indx] = zones.dr
        self.vol[indx] = zones.vol
        self.mass[indx] = zones.mass
        self.cs[indx] = zones.cs
        self.mcum[indx] = zones.mcum

# artificial viscosity
def get_viscosity(zones,v_inner,C_q):
    v1 = np.roll(zones.v, 1)
    v1[0] = v_inner
    dv = v1 - zones.v
    dv = np.where(dv > 0,dv,0)
    q = C_q*zones.d*dv**2
    return q

# volume of zone
def get_vol(zones,r_inner,n):
    r1 = np.roll(zones.r, 1)
    r1[0] = r_inner
    vv = 4.0*pc.pi/3.0*(zones.r**3 - r1**3)
    if np.any(vv[:n] < 0): print('negative volume!')#,vv,r1,zones.r)
    return vv

# function defining inner piston velocity
def get_v_inner(t,piston_in,piston_stop,v_piston,r_inner,pin):
    if (t < piston_in and pin == 1): return -(r_inner-5e7)/0.45
    if ((t >= piston_in or pin != 1) and t < piston_stop): return v_piston
    else: return 0

#def get_v_inner_out(t,piston_in,piston_stop,v_piston,r_inner):
#    if (t < piston_stop): return v_piston
#    else: return 0

def get_e_inner(t,t_stop,edump):
    if (t < t_stop): return edump
    else: return 0

def calc_ffvel(mcum,r,r0,vscale):
        return -(2*pc.G*mcum)**.5*(1./r-1/r0)**.5/vscale

def calc_sedov(self):
        a

def calc_optical_depth(r,d,kappa=0.2):
# calculates optical depth in radial direction from outside in assuming opacity kappa
        rr = r[::-1]
        dtau = kappa*d[::-1]
        tau = -math_ex.tsum(np.log(rr),rr*dtau,cumu=1)
        return tau[::-1]

# end define  zone class
###############################

## Lagrangian 1D hydro class
class lagrange_hydro_1d:

# Initial condition parameters:
    POWERLAW=0
    HEGER=1

# Boundary condition parameters:
    PISTON=0
    INFLOW=1
    OUTFLOW=2
    FALLBACK=3

# HSE model parameters
    NOHSE=-1
    PTSRC=0
    UDENS=1
    HSENUM=2

## -----------------------------------
## setup initial conditions
## just uniform for now
## -----------------------------------

    def init(self):
    ####### parameters ###################
        self.bcinit = 0         # has initialize_boundary_conditions been run?
        self.stinit = 0         # has setup_initial_conditions been run?
        self.itype  = self.POWERLAW   # Default initial condition
        self.udens  = 4.17e-6   # Density for uniform itype
        self.uslope = 0         # Slope for powerlaw itype
        self.ut0    = 1e4       # initial temperature
        self.hse    = self.NOHSE  # Default type of HSE solution
        self.hpath  = 'heger_progenitors/'   # Path to Heger files
        self.htype  = 'solar'   # Heger file type
        self.hname  = 's25.0'    # Heger file name
        self.hnavg  = 10        # Number of points to average within each zone
        self.hdmin  = 0.4       # maximum allowed density drop at last zone -- otherwise cut it
        self.hpscale = 1.       # arbitrary pressure scale for Heger data
        self.hpinfac = 1e-3     # reduction in pressure in inner zone to mimic collapse
        self.t_stop = 2e4       # time (in seconds) to run to
        self.dt_max = 1000       # maximum time step
        self.dt_min = 0         # minimum time step
        self.dt_lim = 0         # if time step goes below this, break
        self.cfl    = 0.1       # courant condition parameter
        self.gamma  = 4.0/3.0   # adiabatic index
        self.C_q    = 5         # artifical viscosity parameter
        self.pin    = 1         # if 1, move piston in first. otherwise just out
        self.piston_in   = 0.45 # time to stop moving the piston in
        self.v_piston    = 2e9  # velocity of inner piston
        self.piston_stop = 1e3  # time to stop moving the inner piston
        self.piston_eexp = 1e52 # stop piston once eint + ke increases by this amount
        self.hent       = 4     # entropy to put piston at
        self.e_dump     = 0     # energy to dump in innermost zone
        self.n_e_dump   = 1     # number of zones to dump energy in
        self.t_e_dump   = 0.1   # time over which to dump energy
        self.r_inner     = 5e11 # initial inner boundary radius
        self.r_outer     = 1e13 # initial outer boundary radius
        self.v_inner     = 0    # inner boundary (piston) velocity
        self.eorig       = 0    # initial internal energy at t=0
        self.nz          = 200  # number of zones - ghost zones
        self.zones       = zone()    # zone data
        self.save        = zone()    # list of saved accreted zones
        self.checkpt     = -1        # checkpointed copy of hydro object
        self.tsave       = np.array([])   # list of saved times
        self.dtsave      = np.array([])   # list of saved dt values
        self.mdot        = np.array([])   # list of saved mdot values
        self.iphot       = np.array([])   # list of saved photosphere indices
        self.vphot       = np.array([])   # list of saved photosphere velocities
        self.tauphot     = np.array([])   # list of saved photosphere opt depths
        self.rphot       = np.array([])   # list of saved photosphere radii
        self.ns     = 0
        self.it     = 0         # timestep number
        self.t      = 0         # time
        self.dt     = 0         # timestep
        self.rmin   = 1e8       # Break if r_inner < rmin
        self.min_zones = 3      # Break if nz < min_zones
        self.mass_r_inner = 0   # Mass inside of inner boundary
        self.bctype = [self.PISTON,self.OUTFLOW]
        self.vscale = 1
        self.noplot = 0
        self.checkint = 1000
        ########################################

    def  __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)

    def setup_initial_conditions(self):
        self.stinit=1
        if self.rmin >= self.r_inner: self.rmin = self.r_inner/5
        if self.itype==self.POWERLAW:
            self.zones.p = np.zeros(self.nz)
            self.zones.v = np.zeros(self.nz)
            self.zones.d = np.zeros(self.nz)
            self.zones.cs = np.zeros(self.nz)
            self.zones.mcum = np.zeros(self.nz)
            self.zones.q = np.zeros(self.nz)
            self.zones.dr = np.zeros(self.nz)
            # set density, pressure, velocity
            self.zones.r = self.r_inner + self.r_outer*(np.arange(self.nz)+1.0)/(1.0*self.nz)
            self.zones.d = self.udens*(self.r_inner/self.zones.r)**(self.uslope)
            self.zones.p += pc.a*(self.ut0**4.0)

            # set volume, mass, etc...

            # volume of each zone
            self.zones.vol = get_vol(self.zones,self.r_inner,self.nz)

            # mass of each zone (this never changes)
            self.zones.mass = self.zones.d*self.zones.vol

            self.zones.mcum = np.cumsum(self.zones.mass)+self.mass_r_inner

            # zone artificial viscosity
            self.zones.q = get_viscosity(self.zones,self.v_inner,self.C_q)

            # flag for HSE runs
            if self.hse != self.NOHSE:
                #rhs=0.
                pmin=pc.a*(self.ut0**4.0)+pc.G*self.mass_r_inner**2/self.r_inner**5
                #for  i in range(self.nz):
                if self.hse==self.PTSRC:
                    self.zones.p = pmin + pc.G*self.mass_r_inner*self.zones.d*(1.0/self.zones.r-1./self.r_inner)
                elif self.hse==self.UDENS:
                        self.zones.p = pc.G*self.zones.d**2*(self.zones.r**2-self.zones.r**3/self.r_inner+self.r_inner**2-self.r_inner**3/self.zones.r)
                else:
                    # Numerically integrate HSE eqn backwards from pmin at outer edge
                    for i in range(self.nz):
                        ii=self.nz-1-i
                        rhs=pc.G*self.zones.mcum[ii]/self.zones.r[ii]**2
                        if ii==0: dr0=self.zones.r[ii]-self.r_inner
                        else: dr0=self.zones.r[ii]-self.zones.r[ii-1]
                        if ii==(self.nz-1):
                            dr1=self.zones.r[ii]-self.zones.r[ii-1]
                            p0=pmin
                            rho1=self.zones.d[ii]
                        else:
                            rho1=self.zones.d[ii+1]
                            p0=self.zones.p[ii+1]
                        mmean = 0.5*(dr0*self.zones.d[ii]+dr1*rho1)
                        self.zones.p[ii]=p0 + mmean*rhs

            # if HSE fix p:
            if self.hse==self.HSENUM: self.zones.p = self.zones.p - self.zones.p[self.nz-1] + pmin

            #gamma-law equation of state
            self.zones.e = self.zones.p/(self.gamma-1.)/self.zones.d

        elif self.itype==self.HEGER:
            s=heger.heger(path=self.hpath,type=self.htype,name=self.hname)
            s.read_heger_star()
            # set inner boundary at s=4:
            if self.hent > 0: 
                s4, = np.where(s.data[heger.ENTR,:] > self.hent)
                self.r_inner = s.data[heger.RAD,s4[0]]
            
            print('r inner: ',self.r_inner)
            indx = np.arange(len(s.data[heger.RAD,:]))[(s.data[heger.RAD,:] > self.r_inner)&(s.data[heger.RAD,:] < self.r_outer)]
            rvals=s.data[heger.RAD,indx]
            pvals=s.data[heger.PRES,indx]
            dvals=s.data[heger.DENS,indx]
            vvals=s.data[heger.VEL,indx]
            # remove outer zone if density drops suddenly:
#            print(len(dvals), len(dvals[:-2]), len(dvals[:-1]))
            if dvals[-1]/dvals[-2] < self.hdmin:
                dvals=dvals[:-1] 
                rvals=rvals[:-1]
                pvals=pvals[:-1]
                vvals=vvals[:-1]
            
            self.nz = len(rvals)
            # set inner/outer radius:
            self.r_inner=2*rvals[0]-rvals[1]
            self.r_outer=2*rvals[-1]-rvals[-2]
            self.zones.r=rvals
            self.zones.p=pvals*self.hpscale
            self.zones.v=vvals
            self.zones.d=dvals
            # initialize other arrays
            self.zones.cs = np.zeros(self.nz)
            self.zones.mcum = np.zeros(self.nz)
            self.zones.q = np.zeros(self.nz)
            self.zones.dr = np.zeros(self.nz)
                #self.zones.r = self.r_inner + self.r_outer*(i + 1.0)/(1.0*self.nz)
                #if i==0: r0=self.r_inner
                #else: r0=self.zones[i-1].r
                # average hnavg interpolants to find zone averaged quantities (minimum hnavg is 2 for each edge):
                #rvals=r0+(self.zones.r-r0)*np.arange(self.hnavg)/float(self.hnavg-1)
                

                #p,d,v,m,e,xh,xhe,xc,xo=s.interpolate_heger_star(rvals)
                # average:
                #self.zones.p=np.sum(p)/self.hnavg*self.hpscale
                #self.zones.d=np.sum(d)/self.hnavg
                #self.zones.v=np.sum(v)/self.hnavg
                #self.zones.e=np.sum(e)/self.hnavg*self.hpscale
                #gamma-law equation of state
            self.zones.e = self.zones.p/(self.gamma-1)/self.zones.d
                #self.zones.p = (self.gamma-1)*self.zones.d*self.zones.e
#set inner mass to mass coordinate of r_inner:
            self.mass_r_inner=s.data[heger.MASS,indx[0]]
            #self.zones.p*=self.hpinfac

# set volume, mass, etc...

            # volume of each zone
            self.zones.vol = get_vol(self.zones,self.r_inner,self.nz)

            # mass of each zone (this never changes)
            self.zones.mass = self.zones.d*self.zones.vol

            self.zones.mcum = np.cumsum(self.zones.mass)+self.mass_r_inner

            # zone artificial viscosity
            self.zones.q = get_viscosity(self.zones,self.v_inner,self.C_q)

        # fill other cells with fluff:
        #if self.r_outer > self.zones.r[self.nz-1]:

    def initialize_boundary_conditions(self):
        self.bcinit=1
        if self.bctype[1]==self.OUTFLOW:
        # append one extra zone for outer boundary ghost
            self.zones.r=np.append(self.zones.r,self.zones.r[self.nz-1])
            self.zones.v=np.append(self.zones.v,self.zones.v[self.nz-1])
            self.zones.p=np.append(self.zones.p,self.zones.p[self.nz-1])
            self.zones.e=np.append(self.zones.e,self.zones.e[self.nz-1])
            self.zones.d=np.append(self.zones.d,self.zones.d[self.nz-1])
            self.zones.mass=np.append(self.zones.mass,self.zones.mass[self.nz-1])
            self.zones.q=np.append(self.zones.q,self.zones.q[self.nz-1])
            self.zones.mcum=np.append(self.zones.mcum,self.zones.mcum[self.nz-1])
            self.zones.cs=np.append(self.zones.cs,self.zones.cs[self.nz-1])
            self.zones.dr=np.append(self.zones.dr,self.zones.dr[self.nz-1])
            self.zones.vol=np.append(self.zones.vol,self.zones.vol[self.nz-1])

    def apply_boundary_conditions(self):
        a_inner = 0
        #print('bc: ',self.bctype,self.INFLOW,self.OUTFLOW)
        if self.bctype[0]==self.INFLOW:
            # Accelerating inner boundary:
            if self.r_inner >= self.rmin: 
                a_inner = -pc.G*self.mass_r_inner / self.r_inner**2
                self.v_inner=self.zones.v[0]
            elif self.r_inner < self.rmin:
                self.r_inner=self.zones.r[0]
                self.v_inner=self.zones.v[0]
                self.mass_r_inner=self.zones.mcum[0]
                # Save accreted zone data & delete from zone arrays
                self.save.add_zone_values(self.zones.get_zone_values(0))
                self.zones.del_zone_values(0)
                self.tsave=np.append(self.tsave,self.t)
                self.dtsave=np.append(self.dtsave,self.dt)
                iphot,tauphot,rphot,vphot = self.locate_photosphere()
                #print('photosphere: ',iphot,tauphot,rphot,vphot)
                self.iphot   = np.append(self.iphot,iphot)
                self.rphot   = np.append(self.rphot,rphot)
                self.vphot   = np.append(self.vphot,vphot)
                self.tauphot = np.append(self.tauphot,tauphot)
                self.nz -= 1
                self.ns += 1

            #print('ainner: ',a_inner,self.mass_r_inner,self.r_inner**2)
        elif self.bctype[0]==self.PISTON:
            # Outgoing outer boundary:
            self.v_inner = get_v_inner(self.t,self.piston_in,self.piston_stop,self.v_piston,self.r_inner,self.pin)
        elif self.bctype[0]==self.FALLBACK:
            # inflow:
            if self.r_inner >= self.rmin: 
                a_inner = -pc.G*self.mass_r_inner / self.r_inner**2
                self.v_inner=self.zones.v[0]
            elif self.r_inner < self.rmin:
                self.r_inner=self.zones.r[0]
                self.v_inner=self.zones.v[0]
                self.mass_r_inner=self.zones.mcum[0]
                # Save accreted zone data & delete from zone arrays
                self.save.add_zone_values(self.zones.get_zone_values(0))
                self.zones.del_zone_values(0)
                self.tsave=np.append(self.tsave,self.t)
                self.dtsave=np.append(self.dtsave,self.dt)
                iphot,tauphot,rphot,vphot = self.locate_photosphere()
#                print('photosphere: ',iphot,tauphot,rphot,vphot)
                self.iphot   = np.append(self.iphot,iphot)
                self.rphot   = np.append(self.rphot,rphot)
                self.vphot   = np.append(self.vphot,vphot)
                self.tauphot = np.append(self.tauphot,tauphot)
                self.nz -= 1
                self.ns += 1
            # piston if t < piston_stop:
            if (self.t < self.piston_stop or self.zones.v[0] > 0): self.v_inner = get_v_inner(self.t,self.piston_in,self.piston_stop,self.v_piston,self.r_inner,self.pin)
        else: print('ERROR: Unrecognized inner BC type!')
        #self.v_inner += a_inner*self.dt
        #print('v_inner: ',self.v_inner)
        if self.bctype[1]==self.OUTFLOW:
            #print('before outflow bc: ', self.zones.r[-3:])
            self.zones.assign_zone(self.nz,self.zones.get_zone_values(self.nz-1))
            #print('after outflwo bc: ', self.zones.r[-3:])
            
        else:
            print('ERROR: Unrecognized outer BC type!')
    
    def get_time_step(self):
        self.dt = 1e99

        # radial width of each zone
        self.zones.dr = self.zones.r - np.roll(self.zones.r, 1)
        self.zones.dr[0] = self.zones.r[0] - self.r_inner
            
        # sound speed
        self.zones.cs = (self.gamma*self.zones.p/self.zones.d)**0.5
        
        # courant condition
        this_dt = np.min(self.cfl*self.zones.dr[:self.nz]/self.zones.cs[:self.nz])
        if (this_dt < self.dt): self.dt = this_dt

        # check inner boundary
        if (self.v_inner > 0):
            this_dt = self.cfl*(self.zones.r[0] - self.r_inner)/self.v_inner
            if (this_dt < self.dt): self.dt = this_dt

        if (self.dt < self.dt_min): self.dt = self.dt_min
        if (self.dt > self.dt_max): self.dt = self.dt_max
    
    def calc_global_quantities(self):
        
        tot_mass = np.sum(self.zones.mass[:-1]) + self.mass_r_inner
        tot_ke = 0.5*np.sum(self.zones.mass[:-1]*self.zones.v[:-1]**2)
        tot_eint = np.sum(self.zones.e[:-1]*self.zones.mass[:-1])
        tot_egrav = -4.*pc.pi*pc.G*np.sum(self.zones.d[:-1]*self.zones.mcum[:-1]*self.zones.dr[:-1]*self.zones.r[:-1])
        tot_egrav_th = -3./5.*pc.G*tot_mass**2./self.zones.r[self.nz-1]
        tot_mass = tot_mass/pc.m_sun
        tot_e = tot_ke + tot_eint + tot_egrav
        return tot_mass, tot_ke, tot_eint, tot_egrav, tot_e

    def update_vel(self):

#        print('ghost: ',self.zones.p[self.nz-1:], self.zones.d[self.nz-1:], self.zones.q[self.nz-1:])
        # pressure and viscosity gradient
        dp = self.zones.p - np.roll(self.zones.p, -1)
        dq = self.zones.q - np.roll(self.zones.q, -1)

        # mean density * dr
        mmean = 0.5*(self.zones.d*self.zones.dr + np.roll(self.zones.d*self.zones.dr, -1))

        # acceleration from pressure and viscosity gradient
        accel = (dp + dq)/mmean;

        # gravitational acceleration
        g = -pc.G*self.zones.mcum/self.zones.r**2

        accel += g

        # update velocities
        self.zones.v[:self.nz-2] += accel[:self.nz-2]*self.dt
        #self.zones.v[self.nz] = self.zones.v[self.nz-1]

        # let outer zone accel be equal to nearest inner zone
        self.zones.v[self.nz-2:] += accel[self.nz-2]*self.dt
        #print('update v: ', self.zones.v[-4:], accel[-4:], dp[-4:], dq[-4:])

    def update_interface_positions(self):
        # do inner boundary
        self.r_inner += self.v_inner*self.dt
        self.zones.r += self.zones.v*self.dt

    def update_quantities(self):

        # get new volume of each zone
        self.zones.vol = get_vol(self.zones,self.r_inner,self.nz)

        # new density of zone
        new_d = self.zones.mass/self.zones.vol

        # new artificial viscosity
        new_q = get_viscosity(self.zones,self.v_inner,self.C_q)
        q_ave = 0.5*(new_q + self.zones.q)

        # volume (i.e., inverse density) change
        dtau = 1.0/new_d - 1.0/self.zones.d;

        # update zone energy (implicitly)
        fact  = 1 + 0.5*dtau*(self.gamma-1)*new_d
        new_e = (self.zones.e - (0.5*self.zones.p + q_ave)*dtau)/fact

        # check for error...
        if np.any(new_e < 0): print("neg energy!")#,self.zones.e,new_e,dtau,new_d,self.zones.p)
        
        # update other quantities
        self.zones.d = new_d
        self.zones.q = new_q
        self.zones.e = new_e

        # dump energy in
        self.zones.e[0:self.n_e_dump] += get_e_inner(self.t, self.t_e_dump, self.e_dump) / (self.zones.mcum[self.n_e_dump]-self.zones.mcum[0]) * (self.dt / self.t_e_dump) / self.n_e_dump

        # equation of state
        self.zones.p = (self.gamma - 1)*new_e*new_d

##############################################
## Auxillary routines for plotting or analysis
##############################################

    def setup_plot(self):
        # set up velocity plot
        x = self.zones.r[:self.nz]
        y = self.zones.v[:self.nz]
        fig = plt.figure()
        this_plot, = plt.semilogx(x,y,label='fluid velocity')
        two_plot,  = plt.semilogx(x,y,label='escape velocity')
#        this_plot = plt.figure(figsize=(6,5))
#        two_plot = plt.figure(figsize=(6,5))
        if self.bctype[0]==self.INFLOW:
            self.ptype='ff'
            this_plot.axes.set_ylim(-0.2,0.2)
            this_plot.axes.set_ylabel('velocity ($10^9$ cm/s)',fontsize=14)
            this_plot.axes.set_xlabel('radius (cm)',fontsize=14)
            two_plot.axes.set_ylim(-0.2,0.2)
            two_plot.axes.set_ylabel('velocity ($10^9$ cm/s)',fontsize=14)
            two_plot.axes.set_xlabel('radius (cm)',fontsize=14)
            if self.hse==self.NOHSE:
                self.vscale = 1e9
            else:
                print(self.zones.v[0], self.vscale)
        elif self.bctype[0]==self.PISTON or self.bctype[0]==self.FALLBACK:
            self.ptype='explode'
            self.vscale = 1e9
            vval=1.5*self.v_piston/self.vscale
            if self.v_piston==0: vval=(2*self.e_dump/self.zones.mcum[0])**0.5/self.vscale
            this_plot.axes.set_ylim(-vval,vval)
            this_plot.axes.set_ylabel('velocity ($10^9$ cm/s)',fontsize=14)
            this_plot.axes.set_xlabel('radius (cm)',fontsize=14)
            two_plot.axes.set_ylim(-vval,vval)
            two_plot.axes.set_ylabel('velocity ($10^9$ cm/s)',fontsize=14)
            two_plot.axes.set_xlabel('radius (cm)',fontsize=14)
        plt.legend(fontsize=12,loc='lower right')
        return x,this_plot,two_plot,fig

    def update_plot(self,x,this_plot,two_plot,fig):
    # update plot
            y = self.zones.v / self. vscale
    # update x if number of zones changed
            if len(x) > self.nz: x=x[len(x)-self.nz:]
            if self.ptype=='ff':
                if self.hse==self.NOHSE: 
                    y2=calc_ffvel(self.zones.mcum[:self.nz],self.zones.r[:self.nz],x,self.vscale)
                    y2=y2[:self.nz]
            elif self.ptype=='explode':
                E0=1e51
                tot_e = self.zones.mass*(self.zones.e + 0.5*self.zones.v**2)
                y2=2./5.*self.zones.r**(-3./2.)*(E0/self.zones.d)**(.5)/self.vscale
                y2=y2[:self.nz]
            y=y[:self.nz]
            
            this_plot.set_ydata(y)
            this_plot.set_xdata(x)
            if self.hse==self.NOHSE:
                two_plot.set_ydata(y2)
                two_plot.set_xdata(x)
            if (self.it % 100 == 0):
                #print('lenxy: ',len(x),len(y),len(y2))
#                plt.draw()
                fig.canvas.draw()
                fig.canvas.flush_events()
                plt.pause(0.0001)

    def calc_mdot(self):
        mdot=(self.save.mass[2:]/2e33/(self.tsave[1:]-self.tsave[:-1]))
        self.mdot=mdot

    def locate_photosphere(self):
        tau = calc_optical_depth(self.zones.r[:-1],self.zones.d[:-1])
        i = np.argmin(np.abs(tau - 1.))
        return i,tau[i],self.zones.r[i],self.zones.v[i]

    def stretch_star(self,i,f):
# stretch a star by a factor of f in radius outward from index i, keeping rho(i) constant and total M constant
        a

############################################
## -----------------------------------------
############################################

    def run_checkpoint(self,times):
# run & save state at times:
        objlist = []
        for i in range(len(times)):
            self.t_stop = times[i]
            self.run()
            self.calc_mdot()
            objlist.append(copy.deepcopy(self))
        return objlist

    def run(self):
# allow for run restarting:        
        if self.stinit==0: self.setup_initial_conditions()
        if self.bcinit==0: self.initialize_boundary_conditions()

# set up plot
        if self.noplot==0: x,p1,p2,fig=self.setup_plot()

        while (1):

    #-----------------------------
    # apply boundary conditions
    #-----------------------------

            self.apply_boundary_conditions()
            
    #-----------------------------
    # get time step
    #-----------------------------
            
            self.get_time_step()
    # check for code crash
            if self.dt <= self.dt_lim:
                print("ERROR -- Time step of zero!")
                break
        
    #---------------------------------
    # calculate global quantities
    #--------------------------------

            tot_mass,tot_ke,tot_eint,tot_egrav,tot_e=self.calc_global_quantities()
            if (self.it % 100 == 0): 
                print("%6d %8.3e %8.3e" % (self.it,self.t,self.dt))
                print("%8.3e %8.3e %8.3e %8.3e %8.3e %8.3e" % (self.mass_r_inner,tot_mass,tot_ke,tot_eint,tot_egrav,tot_e))

    #-----------------------------
    # calculate accelerations
    #-----------------------------

            self.update_vel()

    #---------------------------------
    # update interface positions
    #--------------------------------

            self.update_interface_positions()

    #---------------------------------
    # update other quantities
    #--------------------------------

            self.update_quantities()

        # printout zone properties, if wanted
        #print " %4d" % i,
        #print " %8.3e" % (self.zones.r),
        #print " %8.3e" % (self.zones.vol),
        #print " %8.3e" % (self.zones.v),
        #print " %8.3e" % (self.zones.d),
        #print " %8.3e" % (self.zones.p),
        #print " %8.3e" % (self.zones.e),

            if (self.it % self.checkint == 0):
                self.checkpt = 0
                self.checkpt = copy.deepcopy(self)

            if self.noplot==0: self.update_plot(x,p1,p2,fig)

    # decide whether we need to change BCs:
            if self.it==0: self.eorig = tot_ke + tot_eint
            if (self.t < self.piston_stop) and (tot_eint+tot_ke-self.eorig) > self.piston_eexp: self.piston_stop = self.t

    # advance time
            self.t  += self.dt
            self.it += 1

#            if (self.it > 2): break
            if (self.t > self.t_stop or self.nz < self.min_zones or np.any(self.zones.d[:-1] < 0) or np.any(self.zones.r[:-1] < 0) or np.any(np.isnan(self.zones.d))): break
           # if self.t > self.piston_stop: break
           # if len(self.tsave) > 0.: break
