# -*- coding: utf-8 -*-
# Gating_with_adatoms.py is a program to visualize the electrostatic potential 
# induced by charged defects on surfaces, as they are placed by the user.
# For the moment, the program deals only with the case of adatoms adsorbed on 
# the vacancy sites of the 2x2 reconstruction of the InAs(111)A semiconductor surface. 
#
# Copyright (c) 2015, Jesús Martínez-Blanco.
# All rights reserved.
#
# Author: Jesús Martínez-Blanco
# Email: jmb.jesus@gmail.com
# Licence: BSD

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import scipy.special as spc
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import wx
#from matplotlib.figure import Figure
#from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

class Charged_Defect:
    def __init__(self,charge=1.,tsd=5.,kappa=8.075,redmass=0.043,infinity=2.):
        '''tsd: tip sample distance (in Angstroms)
           kappa: average static dielectric constant (8.075 for InAs)
           redmass: reduced mass (0.043 in our case)
           charge: idem (+1 in our case)
           infinity: a positive value. I try first with 1, then with 2,
               then with 3... until I don't see significant changes in the resulting potential.
               Actually, a value of 1 is already enough.'''
        self.charge = charge
        self.tsd = tsd 
        self._tsd_AU = self.tsd / 0.529177
        self.kappa=kappa
        self.redmass=redmass
        self.infinity=infinity
        
        self.distance = np.linspace(0,300)
        self.pot_array = self.pot_function(self.distance)
        
    def pot_function(self, R):
        '''Returns the electrostatic potential (in mV) as measured at a distance of R Angstroms.
        Performs the full computation for every value of R.'''
        
        r = np.copy(np.asarray(R)) #- make copy to leave R unmodified
        
        if r.ndim == 0: #- in case r is just a scalar, convert to a 1D array for np.newaxis to work later
            r = r.reshape(1)
        elif r.ndim > 1: #- accepts arrays of maximum 1 dimension
            print('pot_function accepts arrays of maximum 1 dimension.')
            return
        
        s = 2*self.redmass/self.kappa

        potk = np.zeros(10000)
        potk_x = np.linspace(0,self.infinity,len(potk))
 
        r /= 0.529177
        
        potk=(self.charge/self.kappa)*(potk_x/(potk_x+s))*spc.j0(potk_x*r[:,np.newaxis])*np.exp(-potk_x*self._tsd_AU)

        area = trapz(potk,potk_x)

        return area*27.211384*1000
    
    def potential(self, R):
        '''Returns the electrostatic potential (in mV) as measured at a distance of R Angstroms.
        The calculation is done by interpolating the values already computed and stored in self.pot_array.'''
        
        f = interp1d(self.distance,self.pot_array,bounds_error=False,fill_value=0.)
        
        return f(R)
        
#- Create pmap (the potential map)
imx = np.linspace(-100,100,100)
imy = np.linspace(-100,100,100)

imXX, imYY = np.meshgrid(imx, imy)  

pmap = np.zeros(imXX.shape)

#- Create the charged defect (a positively charged indium adatom on InAs(111))
charge=1.
tsd=5.
kappa=8.075
redmass=0.043
infinity=2.

CD=Charged_Defect(charge,tsd,kappa,redmass,infinity)

#- Create the grid of vacancy sites
a = 8.57 #- Size of the unit vector
    
x = np.arange(-12.,12.)
y = np.arange(-7.,7.)
XX, YY = np.meshgrid(x, y)
    
XXs = XX + 0.5
YYs = YY + 0.5
    
vac_X = np.concatenate((XX.flatten(),XXs.flatten()))
vac_Y = np.concatenate((YY.flatten(),YYs.flatten()))
    
vac_X *= a
vac_Y *= a*np.sqrt(3)
    
sizeFILLED = 200.
sizeUNFILLED = 20.
vac_sizes = np.zeros_like(vac_X)
vac_sizes[:] = sizeUNFILLED
vac_colors = np.zeros_like(vac_X,dtype=str)
vac_colors[:] = 'r'

#- Create the figure
fig, inas = plt.subplots(figsize=(10,10))

fig.suptitle('Electrostatic Potential Map for a Nanostructure', fontsize=24)
fig.text(0.25,0.9,'Click a vacancy site (red dots) to place an adatom.\nClick an adatom to remove it.', fontsize=14, color='r')

inas.set_xlim(-100,100)
inas.set_ylim(-100,100)
inas.set_aspect(1.)
inas.set_xlabel('X ($\AA$)')
inas.set_ylabel('Y ($\AA$)')

divider = make_axes_locatable(inas)

profH = divider.append_axes("top", 2., pad=0.2, sharex=inas)
profV = divider.append_axes("right", 2., pad=0.2, sharey=inas)
profH.set_ylim(-1,1)
profV.set_xlim(-1,1)
profH.set_ylabel('Potential @ Y$=0$')
profV.set_xlabel('Potential @ X$=0$')

plt.setp(profH.get_xticklabels() + profV.get_yticklabels(), visible=False)# make some labels invisible

#- Add elements to the figure
im = inas.imshow(pmap, extent=[-100, 100, -100, 100], cmap=plt.cm.gray)
imH, = profH.plot(imx,pmap[:,50])
imV, = profV.plot(pmap[50,:],imy)

inas.vlines(0,-100,100,'w',alpha=0.5)
inas.hlines(0,-100,100,'w',alpha=0.5)
potlabel = inas.text(-50, -130, 'Local Potential: 0.0 mV', fontsize=20, color='b')
vacs = inas.scatter(vac_X,vac_Y, picker=5)
vacs._sizes = vac_sizes
vacs.set_facecolor(vac_colors)

#- Add buttons
def ExportData(event):
    np.savetxt('potential_map.dat', pmap.T, fmt='%.2f', header='charge={0},tsd={1},kappa={2},redmass={3},infinity={4}'.format(charge,tsd,kappa,redmass,infinity))
#
#def ImportData(event):
#    print "clicked:", event    
#        
axExport = plt.axes([0.725, 0.75, 0.175, 0.075])
#axImport = plt.axes([0.75, 0.7, 0.1, 0.075])
#
bExport = Button(axExport, 'Export map to\n\"potential_map.dat\"')
bExport.on_clicked(ExportData)
#
#bImport = Button(axImport, 'Import data')
#bImport.on_clicked(ImportData)

#- Add mouse callbacks
def on_pick(event):
    ind = event.ind[0]
    global pmap
    
    if vac_sizes[ind] == sizeFILLED: #- vacancy has already an adatom
        #pmap -= np.exp((-(imXX-vac_X[ind])**2-(-imYY-vac_Y[ind])**2)/(5**2))
        pmap -= CD.potential(np.sqrt((imXX-vac_X[ind])**2+(-imYY-vac_Y[ind])**2))
        vac_sizes[ind] = sizeUNFILLED
    else: #- vacancy is empty, so put an adatom
        #pmap += np.exp((-(imXX-vac_X[ind])**2-(-imYY-vac_Y[ind])**2)/(5**2))
        pmap += CD.potential(np.sqrt((imXX-vac_X[ind])**2+(-imYY-vac_Y[ind])**2))
        vac_sizes[ind] = sizeFILLED
    
    if vac_sizes.max() < sizeFILLED: #- To avoid zero noise when all atoms are removed
        pmap = np.zeros(imXX.shape)
    
    im.set_data(pmap)
    im.autoscale()
    
    imH.set_data(imx,pmap[50,:])
    imV.set_data(pmap[:,50],-imy)

    profH.set_ylim([pmap[50,:].min(),pmap[50,:].max()])    
    profV.set_xlim([pmap[:,50].min(),pmap[:,50].max()])    
    
    if vac_sizes.max() < sizeFILLED: #- To avoid weird ticks numbers in the profile graphs when all atoms are removed
        profH.set_ylim(-1,1)
        profV.set_xlim(-1,1)
    
    vacs._sizes = vac_sizes
    vac_colors[vac_sizes < sizeFILLED] = 'r'
    vac_colors[vac_sizes > sizeUNFILLED] = 'w'
    vacs.set_facecolor(vac_colors)

    fig.canvas.draw()
    
def on_mouse_move(event):
    try:
        pot = pmap[np.searchsorted(imy,-event.ydata),np.searchsorted(imx,event.xdata)]
        potlabel.set_text('Local Potential: {0:.2f} mV'.format(pot))
    except:
        potlabel.set_text('Local Potential: -')
        
    fig.canvas.draw()

fig.canvas.mpl_connect('pick_event', on_pick)
fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)

#- show graph
plt.show()


