
from cmath import nan
import gsw
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import wghc_subroutines as wghc
import Functions as fun
import time
#I define a function to change from latitude and longitude to the indices of the lists




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

x = dataset.variables['LON'][:]  #longitude
y = dataset.variables['LAT'][:]  #latitude
z = dataset.variables['ZAX'][:]  #depth (m)

#print(P) #gives info on P
SIG0 = dataset.variables['SIG0']   #potential density referenced to 0 dbar [depth, latitude,longitude]
GAMN = dataset.variables['GAMMAN'] #approximate neutral surface
PRES = dataset.variables['PRES']   #hydrostatic pressure
TEMP = dataset.variables['TEMP'] #in-situ, in deg C
SALI = dataset.variables['SALINITY'] # in psu
SIG2 = dataset.variables['SIG2'] # potential density referenced to 2000 dbar
BOTD = dataset.variables['BOT_DEP'] #bottom depth (Nz,Ny,Nx)
BOTD = BOTD[0,:,:]  #remove redundant data

# 0 <= x <= 359,5
# -80 <= y <= 90

intx = 2
inty = 1
kA = 27 #initial depth index

imax = int(np.round(np.max(x)/intx))
jmax = int(np.round((np.max(y)-np.min(y))/inty))
lon = np.zeros(imax)
lat = np.zeros(jmax)
i = 0
j = 0
while i < imax:
    lon[i] = i*intx
    i += 1
while j < jmax:
    lat[j] = np.min(y)+j*inty
    j += 1
#Start defining the loop
Hghts = np.zeros((len(x),len(y)))
i = 0
j = 0
while i < len(x):
    while j < len(y):
        Hghts[i,j] = gsw.alpha_on_beta (SALI[kA,j,i], TEMP[kA,j,i], PRES[kA,j,i])
        j += 1
    i += 1
    j = 0

inty = 20
intx = 40
ticksx = np.arange(0,len(x),(intx*len(x))/np.max(x))
labelsx = np.arange(np.min(x),np.max(x),intx)
ticksy = np.arange(0,len(y),(inty*len(y))/(np.max(y)-np.min(y)))
labelsy = np.arange(np.min(y),np.max(y),inty)

plot = plt.pcolormesh(np.transpose(Hghts), cmap='plasma')
plt.xticks(ticksx,labelsx)
plt.yticks(ticksy,labelsy)
#ax.set(xlim = (np.min(x),np.max(x)),ylim= (np.min(y),np.max(x)),scaled)
clb = plt.colorbar(plot)
plt.title("Alpha/Beta on presure: {}".format(PRES[kA,1,1]))
clb.set_label("Alpha over Beta")


#fig = plt.figure()



plt.show()
quit()