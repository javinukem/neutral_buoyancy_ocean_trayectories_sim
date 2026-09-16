"""
Simple code to read and plot WOCE data+simple neutral trajectory calculation
gsw = Gibbs Seawater toolbox in Python
"""
import secrets
from sys import last_traceback
from tkinter import commondialog
import gsw  
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import wghc_subroutines as wghc


# Directories
data_dir = r'C:\Users\Javier\Documents\Proyecto' #need the r at the front otherwise struggles with \
data_fileT = r'\wghc_params.nc'
dataT = data_dir + data_fileT

# Constants
km=1000
day = 24*3600
OM = 2*np.pi/day
RADIUS = 6371*km
PVU = 1.e-07 #reference value of PV

############################################################################################### Read the file
dataset = Dataset(dataT)
#print('')
#print(dataset)
#print('')
#print(dataset.file_format)
#print(dataset.dimensions.keys())
#print(dataset.variables.keys())
#print('')

x = dataset.variables['LON']  #longitude
y = dataset.variables['LAT']  #latitude
z = dataset.variables['ZAX']  #depth (m)
#print(P) #gives info on P
SIG0 = dataset.variables['SIG0']   #potential density referenced to 0 dbar [depth, latitude,longitude]
GAMN = dataset.variables['GAMMAN'] #approximate neutral surface
PRES = dataset.variables['PRES']   #hydrostatic pressure
TEMP = dataset.variables['TEMP'] #in-situ, in deg C
SALI = dataset.variables['SALINITY'] # in psu
SIG2 = dataset.variables['SIG2'] # potential density referenced to 2000 dbar
BOTD = dataset.variables['BOT_DEP'] #bottom depth (Nz,Ny,Nx)
BOTD = BOTD[0,:,:]  #remove redundant data

#### Set cast to be at 55N. 22.5W
i0 = 675  #longitud (angulo azimutal, desde Greenwich). Si estamos en W, ik = 2(360 - long). Si estamos en E, ik = 2long
j0 = 270 #latitud (angulo polar, desde Ecuador). Si estamos en S, jk = 2*(80 - lat). En N, jk = 2*(80+lat).  LATITUDES EMPIEZAN EN 80S
###TANTO LOS VALORES DE LAT COMO LONG VAN DE 0.5 EN 0.5###


#### Compute thermo variables for that cast
pr0 = PRES[:,j0,i0]  # pressure
sg0 = SIG0[:,j0,i0]  # potential density ref to 0 dbar
gm0 = GAMN[:,j0,i0]  # approximate neutral density
sp = SALI[:,j0,i0]  # in situ (practical) salinity
t = TEMP[:,j0,i0]   # in situ temp
sa0 = gsw.SA_from_SP(sp,pr0,x[i0],y[j0]) #"absolute" salinity
ct0 = gsw.CT_from_t(sa0,t,pr0)  #conservative temp
rho0 = gsw.rho(sa0,ct0,pr0) # in situ density

################################################### Compute neutral path from that cast (=0) at A to point B on a neighbouring cast
#choose point A on cast0
kA = 20 #initial depth index
prA = pr0[kA]
ctA = ct0[kA]
saA = sa0[kA]

j = [j0] #defino una lista para las latitudes
i = [i0] #otra lista para longitudes
pr = [-prA] #otra para presiones

k = 0

while k < 15: #iterate to obtain another 10 different values of our trajectory
    #choose neighbouring cast (=1) to be just north of initial cast (arbitrary)
    j1 = j0+1   #next index of latitude
    i1 = i0     #same index of longitude
    j.append(j1)
    i.append(i1)  #we fill the lists of latitudes/longitudes w/ the new values

    Nabv = 2
    pr1 = PRES[kA-Nabv:kA+Nabv+1,j1,i1] # pressure levels centred on kA (Nabv points above, Nabv points below)
    sg = SIG0[kA-Nabv:kA+Nabv+1,j1,i1]  # potential density ref to 0 dbar
    gm = GAMN[kA-Nabv:kA+Nabv+1,j1,i1]  # approximate neutral density
    sp = SALI[kA-Nabv:kA+Nabv+1,j1,i1]  # in situ (practical) salinity
    t = TEMP[kA-Nabv:kA+Nabv+1,j1,i1]   # in situ temp
    sa1 = gsw.SA_from_SP(sp,pr1,x[i1],y[j1]) #"absolute" salinity
    ct1 = gsw.CT_from_t(sa1,t,pr1)  #conservative temp
    #compute location of B on a neighbouring cast (=1)
    prB, dct, dsa, dzz, flag = wghc.calc_gamman_pressureWGHC(ctA,saA,1,prA,ct1,sa1,2*ct1,pr1,False,0)
    ctB = ctA+dct
    saB = saA+dsa
    print('')
    print('****** Neutral trajectory calculation: lon(A) = ',x[i0],' lat(A) = ',y[j0],' pr(A) = ',prA)
    print('****** Neutral trajectory calculation: lon(B) = ',x[i1],' lat(B) = ',y[j1],' pr(B) = ',prB)
    print('****** Neutral trajectory calculation: ct(B)-ct(A) = ',dct,' (K) and sa(B)-sa(A) = ',dsa,' (g/kg)')
    
    prNew = PRES[:, j1, i1] #we select the new cast
    prDif = abs(prNew - prB) #we define the pressure differences b/w our new cast and the measured value
    prDif_list = list(prDif) #change the result to a list element so we can later get an int index

    kA = prDif_list.index(min(prDif_list))  #we find the index in the new cast for which the pressure difference w/ the measured value is minimum.
    
    prA = pr0[kA] ##So, we can redefine the values for our new cast with the new pressure index
    ctA = ct0[kA]
    saA = sa0[kA]

    pr.append(-prB) #fill the pressure list w/ the newest element
    i0 = i1
    j0 = j1 #we set up the new values of lat/long 

    k += 1


####################################################################################### Plots

### Hydrographic profiles at cast0
fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4)  
ax1.plot(ct0,-pr0)
ax1.set_title('CT (deg C)')
ax1.set_ylabel('Depth (m)')
ax2.plot(sa0,-pr0)
ax2.set_title('SA (g/kg)')
ax3.plot(rho0,-pr0)
ax3.set_title('RHO (kg/m3)')
ax4.plot(sg0,-pr0,'b',gm0,-pr0,'g')
ax4.set_title('Sigma_0 and Gamma_n')

### Bathymetry
plt.figure()
cm=plt.pcolor(x[:],y[:],BOTD, vmin = 0, vmax = 5000, cmap='nipy_spectral')
plt.colorbar(cm)
plt.title('Bottom depth (m)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

### Section of neutral surfaces
plt.figure()
cm=plt.pcolor(y[:],-z[:],GAMN[:,:,i0], vmin = 26, vmax = 28, cmap='nipy_spectral')
plt.colorbar(cm)
plt.title('Neutral surfaces at lon = '+str(x[i0]-360))
plt.xlabel('Latitude')
plt.ylabel('Depth')

### Section of salinity with neutral surfaces
plt.figure()
levg = [27, 27.5, 28, 28.1, 28.2, 28.3, 28.4, 28.5, 28.6]
levg1 = np.arange(27.6,29,.1)
levg2 = np.arange(27.65,29,.1)
cm=plt.pcolor(y[:],-z[:],SALI[:,:,i0], vmin = 34.5, vmax = 35.5, cmap='nipy_spectral')
plt.colorbar(cm)
plt.contour(y[:],-z[:], GAMN[:,:,i0],levels=levg1, colors='w')
plt.contour(y[:],-z[:], GAMN[:,:,i0],levels=levg2, colors='orange')
plt.plot(y[j],pr,'ko') #to check that the neutral path from A to B is close to a surface gamn=constant
plt.title('Salinity and Neutral surfaces at lon = '+str(x[i0]-360))
plt.xlabel('Latitude')
plt.ylabel('Depth')


### Section of salinity with sigma2 surfaces
plt.figure()
levg = [36, 36.75, 37, 37.1]
levg1 = np.arange(36,40,.2)
levg2 = np.arange(36.1,40,.2)
cm=plt.pcolor(y[:],-z[:],SALI[:,:,i0], vmin = 34.5, vmax = 35.5, cmap='nipy_spectral')
plt.colorbar(cm)
plt.contour(y[:],-z[:], SIG2[:,:,i0],levels=levg, colors='w')
plt.title('Salinity and sigma2 surfaces at lon = '+str(x[i0]-360))
plt.xlabel('Latitude')
plt.ylabel('Depth')



plt.show()
quit()
