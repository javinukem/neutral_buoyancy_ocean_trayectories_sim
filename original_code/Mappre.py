""" 
Simple code to read and plot WOCE data+simple neutral trajectory calculation
gsw = Gibbs Seawater toolbox in Python
"""
from cmath import nan
from math import cos, sin
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
intx = 80
inty = 40
Rx = 10
Ry = 27
kA = 20 #initial depth index
graphs = 1 #If 1 plotΔpn, Δsan, ΔcTn, Δp0, Δsa0, ΔcT0

iter = 120
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
Hghtsn, Δpn, Δsan, ΔcTn, Hghts0, Δp0, Δsa0, ΔcT0  = np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)) 
st = 2*np.pi/iter
i = 0
j = 0
while i < imax:
    while j < jmax:

        lat0 = lat[j]
        lon0 = lon[i]
        i0,j0 = fun.latlong(x,y,lat0,lon0,1)
        i1 = np.zeros(iter+1)
        j1 = np.zeros(iter+1)
        i2, j2 = fun.latlong(x,y,lat0+Ry,lon0,0)

        k = 1
        while k<=iter:
            #Get the index in the vector
            i1[k],j1[k] = fun.latlong(x,y,lat0 + Ry*cos(st*k),lon0 + Rx*sin(st*k),1)
            if i1[k]!=i1[k-1] or j1[k]!=j1[k-1]:
                i2.append(int(i1[k]))
                j2.append(int(j1[k]))
            k += 1

        ###################################################################################### Plots
        
        prn, cTn, san = fun.trajectorie(i2,j2,PRES,SIG0,GAMN,SALI,TEMP,x,y,kA,False)
        pr0 , cT0, sa0 = fun.trajectorie(i2,j2,PRES,SIG0,GAMN,SALI,TEMP,x,y,kA,1)
        #True = alpha and beta computed for p = 0

        #If we ran into a problem while computing the trajectorie we stop the program
        if len(prn) == 2 and len(pr0) == 2:
            print(i+1, "out of ", imax,"and", j+1, "out of ", jmax, "*** ERROR BOTH***")
            Hghtsn [i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j], Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = nan,nan,nan,nan,nan,nan,nan,nan
        elif len(prn) == 2:
            print(i+1, "out of ", imax,"and", j+1, "out of ", jmax, "*** ERROR FALSE***")
            Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = nan, nan, nan, nan
            Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = pr0[len(pr0)-1]-pr0[0], abs(np.max(pr0)-np.min(pr0)), abs(np.max(sa0)-np.min(sa0)), abs(np.max(cT0)-np.min(cT0)) 
        elif len(pr0) == 2:
            print(i+1, "out of ", imax,"and", j+1, "out of ", jmax, "*** ERROR TRUE***")
            Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = prn[len(prn)-1]-prn[0], abs(np.max(prn)-np.min(prn)), abs(np.max(san)-np.min(san)), abs(np.max(cTn)-np.min(cTn)) 
            Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = nan, nan,nan,nan
        else:
            print(i+1, "out of ", imax,"and", j+1, "out of ", jmax)
            Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = prn[len(prn)-1]-prn[0], abs(np.max(prn)-np.min(prn)), abs(np.max(san)-np.min(san)), abs(np.max(cTn)-np.min(cTn)) 
            Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = pr0[len(pr0)-1]-pr0[0], abs(np.max(pr0)-np.min(pr0)), abs(np.max(sa0)-np.min(sa0)), abs(np.max(cT0)-np.min(cT0)) 
        j += 1
    i += 1
    j = 0

inty = 20
intx = 40
ticksx = np.arange(0,imax,(intx*imax)/np.max(x))
labelsx = np.arange(np.min(x),np.max(x),intx)
ticksy = np.arange(0,jmax,(inty*jmax)/(np.max(y)-np.min(y)))
labelsy = np.arange(np.min(y),np.max(y),inty)
avgf = np.nanmean(Hghtsn)
avgt = np.nanmean(Hghts0)

plt.figure(1)
plot = plt.pcolormesh(np.transpose(Hghtsn), cmap='bwr', vmin = -10, vmax = 10)
plt.xticks(ticksx,labelsx)
plt.yticks(ticksy,labelsy)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title("Pitch of ellipses with Rx = {}, Ry = {} \n  for initial presure = {}. Average = {} \n Neutral calculation".format(Rx,Ry,PRES[kA,0,0],avgf))
clb = plt.colorbar(plot)
clb.set_label('Pitch of trajectory')


plt.figure(2)
plot = plt.pcolormesh(np.transpose(Hghts0), cmap='bwr', vmin = -10, vmax = 10)
plt.xticks(ticksx,labelsx)
plt.yticks(ticksy,labelsy)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title("Pitch of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Average = {} \n Sigma0 calculation".format(Rx,Ry,PRES[kA,0,0],avgt))
clb = plt.colorbar(plot)
clb.set_label('Pitch of trajectory')
#fig = plt.figure()

if graphs == 1:
    plt.figure(3)
    plot = plt.pcolormesh(np.transpose(Δpn), cmap='Blues', vmin = 0, vmax = 10)
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of pressure of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δp = {} \n Neutral calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(Δpn)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of pressure (normalized)')
    
    plt.figure(4)
    plot = plt.pcolormesh(np.transpose(Δp0), cmap='Blues', vmin = 0, vmax = 10)
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of pressure of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δp  = {} \n Sigma 0 calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(Δp0)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of pressure (normalized)')

    plt.figure(5)
    plot = plt.pcolormesh(np.transpose(Δsan), cmap='Blues', vmin = 0, vmax = 10)
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of salinity of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δsa  = {} \n Neutral calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(Δsan)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of salinity (normalized)')

    plt.figure(6)
    plot = plt.pcolormesh(np.transpose(Δsa0), cmap='Blues', vmin = 0, vmax = 10)
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of salinity of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δsa  = {} \n Sigma 0 calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(Δsa0)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of salinity (normalized)')

    plt.figure(7)
    plot = plt.pcolormesh(np.transpose(ΔcTn), cmap='Blues', vmin = 0, vmax = 10)
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of conserved temperature of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max ΔCT  = {} \n Neutral calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(ΔcTn)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of conserved temperature (normalized)')

    plt.figure(8)
    plot = plt.pcolormesh(np.transpose(ΔcT0), cmap='Blues', vmin = 0, vmax = 10)
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of conserved temperature of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max ΔCT  = {} \n Sigma 0 calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(ΔcT0)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of conserved temperature (normalized)')

    plt.figure(9)
    
    plot = plt.pcolormesh(np.transpose(Hghts0)-np.transpose(Hghtsn), cmap='bwr', vmin = -10, vmax = 10)
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Difference of sigma0 and neutral calculation of ellipses  with Rx = {}, Ry = {} \n for initial presure = {}. ".format(Rx,Ry,PRES[kA,0,0]))
    clb = plt.colorbar(plot)
    clb.set_label('Pitch of trajectory')

plt.show(block = 0)

input()

plt.close('all')



quit()

