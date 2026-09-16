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
SALI = dataset.variables['SALINITY']
#we choose a part of the globe to make it easier to look at specific points
xmin, xmax, ymin, ymax = np.min(x), np.max(x), np.min(y), np.max(y)
#xmin, xmax, ymin, ymax = 360-70, 360-10, 0, 60
x2 = [xmin]
y2 = [ymin]
i= 0

#if xmin != np.min(x) and xmax != np.max(x):
while i < np.size(x):
    if x[i] > xmin and x[i] <= xmax:
        x2.append(x[i])
    i += 1

i = 0
#if ymin != np.min(y) and ymax != np.max(y):
while i < np.size(y):
    if y[i] > ymin and y[i] <= ymax:
        y2.append(y[i])
    i += 1

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
intx = 8
inty = 4
Rx = 10
Ry = 18
kA = 23 #initial depth index
graphs = 0 #If 1 plotΔpn, Δsan, ΔcTn, Δp0, Δsa0, ΔcT0
wfile = 0


iter = 120
imax = int(np.round((np.max(x2)-np.min(x2))/intx))
jmax = int(np.round((np.max(y2)-np.min(y2))/inty))
lon = np.zeros(imax)
lat = np.zeros(jmax)
i = 0
j = 0
while i < imax:
    lon[i] = i*intx + xmin
    i += 1
while j < jmax:
    lat[j] = np.min(y2)+j*inty
    j += 1
#Start defining the loop
Hghtsn, Δpn, Δsan, ΔcTn, Hghts0, Δp0, Δsa0, ΔcT0  = np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)) 
st = 2*np.pi/iter
i = 0
j = 0
if wfile == 1:
    file = open("ERROR_RX{}_RY{}_PR{}kA{}.txt".format(Rx,Ry,PRES[kA,0,0],kA),'w')
while i < imax:
    while j < jmax:

        lat0 = lat[j]
        lon0 = lon[i]
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
        
        prn, cTn, san, flagn, calcflagn = fun.trajectorie(i2,j2,PRES,SIG0,GAMN,SALI,TEMP,x,y,kA,False)
        pr0 , cT0, sa0, flag0, calcflag0 = fun.trajectorie(i2,j2,PRES,SIG0,GAMN,SALI,TEMP,x,y,kA,1)
        #True = alpha and beta computed for p = 0

        #If we ran into a problem while computing the trajectorie we stop the program
        if flagn == 1 and flag0 == 1:
            if wfile == 1:
                str =  "Longitude {} and latitude {} *** ERROR BOTH*** Neutral : {} Sigma0 : {} \n".format(lon0,lat0,calcflagn,calcflag0)
                file.write(str)
            print(i+1, "out of ", imax,"and", j+1, "out of ", jmax, "*** ERROR BOTH*** Neutral : ", calcflagn, " Sigma0 : ", calcflag0)
            Hghtsn [i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j], Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = nan,nan,nan,nan,nan,nan,nan,nan
        elif flagn == 1:
            if wfile == 1:
                str =  "Longitude {} and latitude {} *** ERROR NTRL*** Neutral : {} \n".format(lon0,lat0,calcflagn)
                file.write(str)
            print(i+1, "out of ", imax,"and", j+1, "out of ", jmax, "*** ERROR NTRL*** Neutral : ", calcflagn)
            Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = nan, nan, nan, nan
            Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = pr0[len(pr0)-1]-pr0[0], abs(np.max(pr0)-np.min(pr0)), abs(np.max(sa0)-np.min(sa0)), abs(np.max(cT0)-np.min(cT0)) 
        elif flag0 == 1:
            if wfile == 1:
                str =  "Longitude {} and latitude {} *** ERROR SGM0*** Sigma0 : {} \n".format(lon0,lat0,calcflag0)
                file.write(str)
            print(i+1, "out of ", imax,"and", j+1, "out of ", jmax, "*** ERROR SGM0*** Sigma0 : ", calcflag0)
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
intx = 30
ticksx = np.arange(0,imax,(intx*imax)/(np.max(x2)-np.min(x2)))
labelsx = np.arange(np.min(x2),np.max(x2),intx)
ticksy = np.arange(0,jmax,(inty*jmax)/(np.max(y2)-np.min(y2)))
labelsy = np.arange(np.min(y2),np.max(y2),inty)
avgf = np.nanmean(Hghtsn)
avgt = np.nanmean(Hghts0)

SALI = SALI[0,:,:] > 0

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)
ax1.set_facecolor('k')
plt.contourf(x,y,SALI, colors ='w')
plot = plt.pcolor(lon,lat,np.transpose(Hghtsn), cmap='bwr', vmin = -10, vmax = 10)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title("Pitch of ellipses with Rx = {}, Ry = {} \n  for initial presure = {}. Average = {} \n Neutral calculation".format(Rx,Ry,PRES[kA,0,0],avgf))
clb = plt.colorbar(plot)
clb.set_label('Pitch of trajectory')
ax1.set(xlim=(xmin,xmax))
ax1.set(ylim=(ymin,ymax))

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)
ax2.set_facecolor('k')
plt.contourf(x,y,SALI, colors ='w')
plot = plt.pcolor(lon,lat,np.transpose(Hghts0), cmap='bwr', vmin = -10, vmax = 10)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title("Pitch of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Average = {} \n Sigma0 calculation".format(Rx,Ry,PRES[kA,0,0],avgt))
clb = plt.colorbar(plot)
clb.set_label('Pitch of trajectory')
ax2.set(xlim=(xmin,xmax))
ax2.set(ylim=(ymin,ymax))


if graphs == 1:
    plt.figure(3)
    plot = plt.pcolormesh(np.transpose(Δpn), cmap='Blues')
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of pressure of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δp = {} \n Neutral calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(Δpn)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of pressure ')
    
    plt.figure(4)
    plot = plt.pcolormesh(np.transpose(Δp0), cmap='Blues')
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of pressure of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δp  = {} \n Sigma 0 calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(Δp0)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of pressure ')

    plt.figure(5)
    plot = plt.pcolormesh(np.transpose(Δsan), cmap='Blues')
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of salinity of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δsa  = {} \n Neutral calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(Δsan)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of salinity ')

    plt.figure(6)
    plot = plt.pcolormesh(np.transpose(Δsa0), cmap='Blues')
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of salinity of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δsa  = {} \n Sigma 0 calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(Δsa0)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of salinity ')

    plt.figure(7)
    plot = plt.pcolormesh(np.transpose(ΔcTn), cmap='Blues')
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of conserved temperature of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max ΔCT  = {} \n Neutral calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(ΔcTn)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of conserved temperature ')

    plt.figure(8)
    plot = plt.pcolormesh(np.transpose(ΔcT0), cmap='Blues')
    plt.xticks(ticksx,labelsx)
    plt.yticks(ticksy,labelsy)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title("Range of conserved temperature of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max ΔCT  = {} \n Sigma 0 calculation".format(Rx,Ry,PRES[kA,0,0],np.nanmax(ΔcT0)))
    clb = plt.colorbar(plot)
    clb.set_label('Range of conserved temperature')


plt.show(block = 0)

input()

plt.close('all')

quit()

