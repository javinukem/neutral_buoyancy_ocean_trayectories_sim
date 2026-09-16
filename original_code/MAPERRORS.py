''' 
Simple code to read and plot WOCE data+simple neutral trajectory calculation
gsw = Gibbs Seawater toolbox in Python
'''
from cmath import nan
from math import cos, sin
import matplotlib
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import Functions as fun
from matplotlib import cm
import os
import gsw

#I define a function to change from latitude and longitude to the indices of the lists

intx = 4
inty = 2
Rx = 10
Ry = 25
kAmax = 25 #initial depth index
wfile = 1
Nabv = 7
ellipse = 0
square = 1
iter = 140
graphs = 1
folder = 1

# Directories
data_dir = 'C:/Users/Javier/Documents/Proyecto' #need the r at the front otherwise struggles with \
data_fileT = '/wghc_params.nc'
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

SIG0 = dataset.variables['SIG0']   #potential density referenced to 0 dbar [depth, latitude,longitude]
GAMN = dataset.variables['GAMMAN'] #approximate neutral surface
PRES = dataset.variables['PRES']   #hydrostatic pressure
TEMP = dataset.variables['TEMP'] #in-situ, in deg C
SALI = dataset.variables['SALINITY'] # in psu
SIG2 = dataset.variables['SIG2'] # potential density referenced to 2000 dbar
BOTD = dataset.variables['BOT_DEP'] #bottom depth (Nz,Ny,Nx)
BOTD = BOTD[0,:,:]  #remove redundant data

#CREATE A PATH TO SAVE FIGURES AND TXT


if folder == 1:
    if ellipse == 1:
        pathname = 'E_IX{}IY{}_RX{}RY{}_NABV{}_kA{}P0{}'.format(intx,inty,Rx,Ry,Nabv,kAmax,PRES[kAmax,0,0])
    if square == 1:
        pathname = 'S_IX{}IY{}_RX{}RY{}_NABV{}_kA{}P0{}'.format(intx,inty,Rx,Ry,Nabv,kAmax,PRES[kAmax,0,0])
    if os.path.isdir(os.path.join(data_dir,pathname)):
        i = 1
        pathname = pathname + str(i)
        while os.path.isdir(os.path.join(data_dir,pathname)):
            i += 1
            pathname = pathname[:-1]+ str(i)
        os.mkdir(pathname)
    else:
        os.mkdir(pathname)
else:
    pathname = ''

#we choose a part of the globe to make it easier to look at specific points
xmin, xmax, ymin, ymax = np.min(x), np.max(x), np.min(y), np.max(y)
#xmin, xmax, ymin, ymax = 360-70, 360-10, 0, 60
x2, y2 = [xmin], [ymin]
i= 0

#if xmin != np.min(x) and xmax != np.max(x):
while i < np.size(x):
    if x[i] >= xmin and x[i] <= xmax:
        x2.append(x[i])
    i += 1

i = 0
#if ymin != np.min(y) and ymax != np.max(y):
while i < np.size(y):
    if y[i] >= ymin and y[i] <= ymax:
        y2.append(y[i])
    i += 1

#print(P) #gives info on P


# 0 <= x <= 359,5
# -80 <= y <= 90

if ellipse == 1:
    st = 2*np.pi/iter
if square == 1:
    iter = 4*(Rx+Ry)

imax, jmax = int(np.round((np.max(x2)-np.min(x2))/intx)), int(np.round((np.max(y2)-np.min(y2))/inty))
lon,lat = np.zeros(imax), np.zeros(jmax)
i, j = 0, 0
bn, b0 = ' BELOW NEUTRAL', ' BELOW SIGMA0'
while i < imax:
    lon[i] = i*intx + xmin
    i += 1
while j < jmax:
    lat[j] = np.min(y2)+j*inty
    j += 1
#Start defining the loop
Hghtsn, Δpn, Δsan, ΔcTn, Hghts0, Δp0, Δsa0, ΔcT0  = np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)), np.zeros((imax,jmax)) 
errorsn, errors0, belown, below0  = np.empty((imax,jmax)), np.empty((imax,jmax)), np.empty((imax,jmax)), np.empty((imax,jmax))
errorsn[:], errors0[:], belown[:], below0[:] = np.NaN, np.NaN, np.NaN, np.NaN
i, j = 0, 0
if wfile == 1: 
    if square == 1:
        filefrst = 'ERRORS RX {} // RY {} // INTX {} // INTY {} // PR{} // kA{} // Nabv{} // ITER {} // SQUARE.txt'.format(Rx,Ry,intx,inty,PRES[kAmax,0,0],kAmax,Nabv,iter)
        file = open(os.path.join(data_dir,pathname,'ERRORS.txt'),'w')
        file.write(filefrst + '\n')
    if ellipse ==1:
        filefrst = 'ERRORS RX {} // RY {} // INTX {} // INTY {} // PR{} // kA{} // Nabv{} // ITER {} // ELLIPSE.txt'.format(Rx,Ry,intx,inty,PRES[kAmax,0,0],kAmax,Nabv,iter)
        file = open(os.path.join(data_dir,pathname,'ERRORSE.txt'),'w')
        file.write(filefrst + '\n')
while i < imax:
    while j < jmax:
        lat0, lon0 = lat[j], lon[i]
        k = 1
        if ellipse == 1:
            i1,j1 = np.zeros(iter+1), np.zeros(iter+1)
            i2, j2 = fun.latlong(x,y,lat0+Ry,lon0,0)
            while k<=iter:
                i1[k],j1[k] = fun.latlong(x,y,lat0 + Ry*cos(st*k),lon0 + Rx*sin(st*k),1)
                if i1[k]!=i1[k-1] or j1[k]!=j1[k-1]:
                    i2.append(int(i1[k]))
                    j2.append(int(j1[k]))
                k += 1    
        if square ==1 :
            i1,j1 = np.zeros(iter+1), np.zeros(iter+1)
            i2, j2 = fun.latlong(x,y,lat0+Ry/2,lon0-Rx/2,0)
            while k <= iter:
                if k <= 2*Rx:
                    i1[k],j1[k] = fun.latlong(x,y,lat0 + Ry/2,lon0 - Rx/2 + k*0.5,1)
                elif k <= 2*(Rx+Ry):
                    i1[k],j1[k] = fun.latlong(x,y,lat0 + Ry/2 - (k-2*Rx)*0.5,lon0 + Rx/2 ,1)
                elif k <= 2*(2*Rx+Ry):
                    i1[k],j1[k] = fun.latlong(x,y,lat0 - Ry/2,lon0 + Rx/2 - (k-2*(Rx+Ry))*0.5,1)
                else:
                    i1[k],j1[k] = fun.latlong(x,y,lat0 - Ry/2 + (k-2*(2*Rx+Ry))*0.5,lon0 - Rx/2 ,1)
                if i1[k]!=i1[k-1] or j1[k]!=j1[k-1]:
                    i2.append(int(i1[k]))
                    j2.append(int(j1[k]))
                k += 1    

        ###################################################################################### Plots
        
        prn, cTn, san, flagn, calcflagn = fun.trajectorie(i2,j2,PRES,SIG0,GAMN,SALI,TEMP,x,y,kAmax,False,Nabv)
        pr0 , cT0, sa0, flag0, calcflag0 = fun.trajectorie(i2,j2,PRES,SIG0,GAMN,SALI,TEMP,x,y,kAmax,1,Nabv)
        #True = alpha and beta computed for p = 0


        #If we ran into a problem while computing the trajectorie we stop the program
        if flagn == 1 and flag0 == 1:
            if wfile == 1:
                prn = prn[np.where(prn != 0)]
                pr0 = pr0[np.where(prn != 0)]
                prDif_listn = list(abs(PRES[:, j2[len(prn)], i2[len(prn)]] + prn[-1]))    
                kAn = prDif_listn.index(min(prDif_listn))
                prDif_list0 = list(abs(PRES[:, j2[len(pr0)], i2[len(pr0)]] + pr0[-1]))    
                kA0 = prDif_list0.index(min(prDif_list0))  
                str =  'Longitude {} and latitude {} *** ERROR BOTH*** Neutral : {} Sigma0 : {} '.format(lon0,lat0,calcflagn,calcflag0)
                if fun.check(i2,j2,PRES,BOTD,kAn,Nabv,prn):
                    str = str + bn
                    if calcflagn == 'no_data':
                        belown[i,j] = 1/6
                    elif calcflagn == 'ok':
                        belown[i,j] = 3/6
                    else:
                        belown[i,j] = 5/6
                if fun.check(i2,j2,PRES,BOTD,kA0,Nabv,pr0):
                    str = str + b0
                    if calcflag0 == 'no_data':
                        below0[i,j] = 1/6
                    elif calcflag0 == 'ok':
                        below0[i,j] = 3/6
                    else:
                        below0[i,j] = 5/6
                str = str + '\n'
                file.write(str)
            if calcflagn == 'no_data':
                errorsn[i,j] = 1/6
            elif calcflagn == 'ok':
                errorsn[i,j] = 3/6
            else:
                errorsn[i,j] = 5/6
            if calcflag0 == 'no_data':
                errors0[i,j] = 1/6
            elif calcflag0 == 'ok':
                errors0[i,j] = 3/6
            else:
                errors0[i,j] = 5/6
            
            print(i+1, 'out of ', imax,'and', j+1, 'out of ', jmax, '*** ERROR BOTH*** Neutral : ', calcflagn, ' Sigma0 : ', calcflag0)
            Hghtsn [i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j], Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = nan,nan,nan,nan,nan,nan,nan,nan
        elif flagn == 1:
            if wfile == 1:
                prn = prn[np.where(prn != 0)]
                str =  'Longitude {} and latitude {} *** ERROR NTRL*** Neutral : {} '.format(lon0,lat0,calcflagn)
                prDif_listn = list(abs(PRES[:, j2[len(prn)], i2[len(prn)]] + prn[-1]))    
                kAn = prDif_listn.index(min(prDif_listn))
                if fun.check(i2,j2,PRES,BOTD,kAn,Nabv,prn):
                    str = str + bn
                    if calcflagn == 'no_data':
                        belown[i,j] = 1/6
                    elif calcflagn == 'ok':
                        belown[i,j] = 3/6
                    else:
                        belown[i,j] = 5/6
                str = str + '\n'
                file.write(str)
            print(i+1, 'out of ', imax,'and', j+1, 'out of ', jmax, '*** ERROR NTRL*** Neutral : ', calcflagn)
            if calcflagn == 'no_data':
                errorsn[i,j] = 1/6
            elif calcflagn == 'ok':
                errorsn[i,j] = 3/6
            else:
                errorsn[i,j] = 5/6
            Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = nan, nan, nan, nan
            Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = pr0[-1]-pr0[0], abs(np.max(pr0)-np.min(pr0)), abs(np.max(sa0)-np.min(sa0)), abs(np.max(cT0)-np.min(cT0)) 
        elif flag0 == 1:
            if wfile == 1:
                pr0 = pr0[np.where(pr0 != 0)]
                print(i+1, 'out of ', imax,'and', j+1, 'out of ', jmax, '*** ERROR SGM0*** Sigma0 : ', calcflag0)
                prDif_list0 = list(abs(PRES[:, j2[len(pr0)], i2[len(pr0)]] + pr0[-1]))    
                kA0 = prDif_list0.index(min(prDif_list0))  
                str =  'Longitude {} and latitude {} *** ERROR SGM0 *** Sigma0 : {} '.format(lon0,lat0,calcflagn)
                if fun.check(i2,j2,PRES,BOTD,kA0,Nabv,pr0):
                    str = str + b0
                    if calcflag0 == 'no_data':
                        below0[i,j] = 1/6
                    elif calcflag0 == 'ok':
                        below0[i,j] = 3/6
                    else:
                        below0[i,j] = 5/6
                str = str + '\n'
                file.write(str)
            if calcflag0 == 'no_data':
                errors0[i,j] = 1/6
            elif calcflag0 == 'ok':
                errors0[i,j] = 3/6
            else:
                errors0[i,j] = 5/6
            Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = prn[-1]-prn[0], abs(np.max(prn)-np.min(prn)), abs(np.max(san)-np.min(san)), abs(np.max(cTn)-np.min(cTn)) 
            Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = nan, nan,nan,nan
        else:
            print(i+1, 'out of ', imax,'and', j+1, 'out of ', jmax)
            Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = prn[-1]-prn[0], abs(np.max(prn)-np.min(prn)), abs(np.max(san)-np.min(san)), abs(np.max(cTn)-np.min(cTn)) 
            Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = pr0[-1]-pr0[0], abs(np.max(pr0)-np.min(pr0)), abs(np.max(sa0)-np.min(sa0)), abs(np.max(cT0)-np.min(cT0)) 
        j += 1
    i += 1
    j = 0

if wfile == 1:
    numerrorsn, numerrors0 = np.count_nonzero(~np.isnan(errorsn)), np.count_nonzero(~np.isnan(errors0))
    nodatan, okn, critn = len((errorsn[np.where(errorsn == 1/6)])), len((errorsn[np.where(errorsn == 3/6)])), len((errorsn[np.where(errorsn == 5/6)]))
    nodatanb, oknb, critnb = len((belown[np.where(belown == 1/6)])), len((belown[np.where(belown == 3/6)])), len((belown[np.where(belown == 5/6)]))
    nodata0, ok0, crit0 = len((errors0[np.where(errors0 == 1/6)])), len((errors0[np.where(errors0 == 3/6)])), len((errors0[np.where(errors0 == 5/6)]))
    nodata0b, ok0b, crit0b = len((below0[np.where(below0 == 1/6)])), len((below0[np.where(below0 == 3/6)])), len((below0[np.where(below0 == 5/6)]))
    file.write('***********NEUTRAL***********\n')
    file.write('TOTAL : NO DATA  {} OK  {} CRITNEGATIVE : {}  \n'.format(nodatan,okn,critn))
    file.write('BELOW : NO DATA  {} OK  {} CRITNEGATIVE : {}  \n'.format(nodatanb,oknb,critnb))
    file.write('***********SIGMA0***********\n')
    file.write('TOTAL : NO DATA  {} OK  {} CRITNEGATIVE : {}  \n'.format(nodata0,ok0,crit0))
    file.write('BELOW : NO DATA  {} OK  {} CRITNEGATIVE : {}  \n'.format(nodata0b,ok0b,crit0b))
    file.write('NUMBER OF NEUTRAL ERRORS {} NUMBER OF SIGMA0 ERRORS {} OUT OF {}'.format(numerrorsn,numerrors0,(imax+1)*(jmax+1)))

intx, inty = 30, 20
ticksx , ticksy = np.arange(0,imax,(intx*imax)/(np.max(x2)-np.min(x2))),  np.arange(0,jmax,(inty*jmax)/(np.max(y2)-np.min(y2)))
labelsx, labelsy = np.arange(np.min(x2),np.max(x2),intx), np.arange(np.min(y2),np.max(y2),inty)
avgf, avgt, avgs = np.nanmean(Hghtsn), np.nanmean(Hghts0), np.nanmean(-Hghts0+Hghtsn)
avgf, avgt, avgs = "{:.1f}".format(avgf),"{:.1f}".format(avgt),"{:.1f}".format(avgs)

if graphs == 1:
    i,j = 0,0
    Hghts = np.zeros((len(x),len(y)))
    while i < len(x):
        while j < len(y):
            print(i+1, 'out of', len(x), ' ', j+1 , ' out of ', len(y))
            Hghts[i,j] = gsw.alpha_on_beta (SALI[kAmax,j,i], TEMP[kAmax,j,i], PRES[kAmax,j,i])
            j += 1
        i += 1
        j = 0

SALI = SALI[0,:,:] > 0

if square == 1:
    shape = 'square'
if ellipse == 1:
    shape = 'elliptical'

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)
ax1.set_facecolor('k')
plt.contourf(x,y,SALI, colors ='w')
plot = plt.pcolor(lon,lat,np.transpose(Hghtsn), cmap='bwr', vmin = -10, vmax = 10)
plt.xlabel('Longitude / º')
plt.ylabel('Latitude / º')
fig1.suptitle('Pitch of {} neutral buoyancy trajectories with Rx = {} º, Ry = {} º \n  for initial height = {} m. Average = {} m \n Neutral calculation'.format(shape,Rx,Ry,PRES[kAmax,0,0],avgf))
clb = plt.colorbar(plot)
clb.set_label('Pitch of trajectory / m (Normalized)')
ax1.set(xlim=(xmin,xmax))
ax1.set(ylim=(ymin,ymax))

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)
ax2.set_facecolor('k')
plt.contourf(x,y,SALI, colors ='w')
plot = plt.pcolor(lon,lat,np.transpose(Hghts0), cmap='bwr', vmin = -10, vmax = 10)
plt.xlabel('Longitude / º')
plt.ylabel('Latitude / º')
fig2.suptitle('Pitch of {} neutral buoyancy trajectories with Rx = {} º, Ry = {} º\n for initial height = {} m. Average = {} m\n Sigma0 calculation'.format(shape,Rx,Ry,PRES[kAmax,0,0],avgt))
clb = plt.colorbar(plot)
clb.set_label('Pitch of trajectory / m (Normalized)')
ax2.set(xlim=(xmin,xmax))
ax2.set(ylim=(ymin,ymax))

fig3 = plt.figure(3)
ax3 = fig3.add_subplot(1,1,1)
ax3.set_facecolor('k')
clmaperror = cm.get_cmap('brg', 3)
plt.contourf(x,y,SALI, colors ='w')
plot = plt.pcolor(lon,lat,np.transpose(errorsn), cmap = clmaperror,vmin = 0, vmax = 1)
plt.xlabel('Longitude / º')
plt.ylabel('Latitude / º')
fig3.suptitle('Error codes with Rx = {} º, Ry = {} º \n for initial height = {}\n Neutral calculation'.format(Rx,Ry,PRES[kAmax,0,0]))
clb = plt.colorbar(plot)
clb.set_ticks([1/6,3/6,5/6])
clb.set_ticklabels(['No data', 'N2 neg', 'Crit neg'])
ax3.set(xlim=(xmin,xmax))
ax3.set(ylim=(ymin,ymax))

fig4 = plt.figure(4)
ax4 = fig4.add_subplot(1,1,1)
ax4.set_facecolor('k')
plt.contourf(x,y,SALI, colors ='w')
plot = plt.pcolor(lon,lat,np.transpose(errors0), cmap = clmaperror,vmin = 0, vmax = 1)
plt.xlabel('Longitude / º')
plt.ylabel('Latitude / º')
fig4.suptitle('Error codes with Rx = {} º, Ry = {} º\n for initial height = {}\n Sigma0 calculation'.format(Rx,Ry,PRES[kAmax,0,0]))
clb = plt.colorbar(plot)
clb.set_ticks([1/6,3/6,5/6])
clb.set_ticklabels(['No data', 'N2 neg', 'Crit neg'])
ax4.set(xlim=(xmin,xmax))
ax4.set(ylim=(ymin,ymax))



fig5 = plt.figure(5)
ax5 = fig5.add_subplot(1,1,1)
ax5.set_facecolor('k')
plt.contourf(x,y,SALI, colors ='w')
plot = plt.pcolor(lon,lat,np.transpose(-Hghts0+Hghtsn), cmap='bwr', vmin = -10, vmax = 10)
plt.xlabel('Longitude / º')
plt.ylabel('Latitude / º')
fig5.suptitle('Substraction of neutral minus sigma0 calculation for \n {} neutral buoyancy trajectories. Rx = {} º, Ry = {} º\n for initial height = {} m. Average = {} m'.format(shape,Rx,Ry,PRES[kAmax,0,0],avgs))
clb5 = plt.colorbar(plot)
clb5.set_label('Pitch of trajectory / m (Normalized)')
ax5.set(xlim=(xmin,xmax))
ax5.set(ylim=(ymin,ymax))


fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig5.tight_layout()

if folder == 1: 

    fig1.savefig(pathname + '/mapn.png')
    fig2.savefig(pathname + '/map0.png')
    fig3.savefig(pathname + '/errn.png')
    fig4.savefig(pathname + '/err0.png')
    fig5.savefig(pathname + '/mapsub.png')

if graphs == 1:

    fig6 = plt.figure(6)
    ax6 = fig6.add_subplot(1,1,1)
    ax6.set_facecolor('k')
    ax6.set(xlim=(xmin,xmax))
    ax6.set(ylim=(ymin,ymax))
    plt.contourf(x,y,SALI, colors ='w')
    plot = plt.pcolor(lon,lat,np.transpose(Δpn), cmap='Blues')
    plt.xlabel('Longitude / º')
    plt.ylabel('Latitude / º')
    fig6.suptitle("Range of heights of {} neutral trajectories with Rx = {} º, Ry = {} º \n for initial height = {} m. Max Δz = {} m\n Neutral calculation".format(shape,Rx,Ry,PRES[kAmax,0,0],"{:.1f}".format(np.nanmax(Δpn))))
    clb = plt.colorbar(plot)
    clb.set_label('Δz / m')
    
    fig7 = plt.figure(7)
    ax7 = fig7.add_subplot(1,1,1)
    ax7.set_facecolor('k')
    ax7.set(xlim=(xmin,xmax))
    ax7.set(ylim=(ymin,ymax))
    plt.contourf(x,y,SALI, colors ='w')    
    plot = plt.pcolor(lon,lat,np.transpose(Δp0), cmap='Blues')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig7.suptitle("Range of pressure of {} neutral trajectories with Rx = {}, Ry = {} \n for initial presure = {}. Max Δp  = {} \n Sigma 0 calculation".format(shape,Rx,Ry,PRES[kAmax,0,0],"{:.1f}".format(np.nanmax(Δp0))))
    clb = plt.colorbar(plot)
    clb.set_label('Range of pressure ')

    fig8 = plt.figure(8)
    ax8 = fig8.add_subplot(1,1,1)
    ax8.set_facecolor('k')
    ax8.set(xlim=(xmin,xmax))
    ax8.set(ylim=(ymin,ymax))
    plt.contourf(x,y,SALI, colors ='w')
    plot = plt.pcolor(lon,lat,np.transpose(Δsan), cmap='Blues')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig8.suptitle("Range of salinity of {} neutral trajectories with Rx = {}, Ry = {} \n for initial presure = {}. Max Δsa  = {} \n Neutral calculation".format(shape,Rx,Ry,PRES[kAmax,0,0],"{:.1f}".format(np.nanmax(Δsan))))
    clb = plt.colorbar(plot)
    clb.set_label('Range of salinity ')

    fig9 = plt.figure(9)
    ax9 = fig9.add_subplot(1,1,1)
    ax9.set_facecolor('k')
    ax9.set(xlim=(xmin,xmax))
    ax9.set(ylim=(ymin,ymax))
    plt.contourf(x,y,SALI, colors ='w')
    plot = plt.pcolor(lon,lat,np.transpose(Δsa0), cmap='Blues')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig9.suptitle("Range of salinity of {} neutral trajectories with Rx = {}, Ry = {} \n for initial presure = {}. Max Δsa  = {} \n Sigma 0 calculation".format(shape,Rx,Ry,PRES[kAmax,0,0],"{:.1f}".format(np.nanmax(Δsa0))))
    clb = plt.colorbar(plot)
    clb.set_label('Range of salinity ')

    fig10 = plt.figure(10)
    ax10 = fig10.add_subplot(1,1,1)
    ax10.set_facecolor('k')
    ax10.set(xlim=(xmin,xmax))
    ax10.set(ylim=(ymin,ymax))
    plt.contourf(x,y,SALI, colors ='w')
    plot = plt.pcolor(lon,lat,np.transpose(ΔcTn), cmap='Blues')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig10.suptitle("Range of conserved temperature of {} neutral trajectories with Rx = {}, Ry = {} \n for initial presure = {}. Max ΔCT  = {} K \n Neutral calculation".format(shape,Rx,Ry,PRES[kAmax,0,0],"{:.1f}".format(np.nanmax(ΔcTn))))
    clb = plt.colorbar(plot)
    clb.set_label('Range of conserved temperature / K')

    fig11 = plt.figure(11)
    ax11 = fig11.add_subplot(1,1,1)
    ax11.set_facecolor('k')
    ax11.set(xlim=(xmin,xmax))
    ax11.set(ylim=(ymin,ymax))
    plt.contourf(x,y,SALI, colors ='w')
    plot = plt.pcolor(lon,lat,np.transpose(ΔcT0), cmap='Blues')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig11.suptitle("Range of conserved temperature of {} neutral trajectories with Rx = {}, Ry = {} \n for initial presure = {}. Max ΔCT  = {} \n Sigma 0 calculation".format(shape,Rx,Ry,PRES[kAmax,0,0],"{:.1f}".format(np.nanmax(ΔcT0))))
    clb = plt.colorbar(plot)
    clb.set_label('Range of conserved temperature')

    fig12 = plt.figure(12)
    ax12 = fig12.add_subplot(1,1,1)
    ax12.set(xlim=(xmin,xmax))
    ax12.set(ylim=(ymin,ymax))
    clb = plt.colorbar(plot)
    fig12.suptitle("Alpha/Beta on presure: {}".format(PRES[kAmax,1,1]))
    clb.set_label("Alpha over Beta")
    plot = plt.pcolor(x,y,np.transpose(Hghts), cmap='plasma')

    if folder == 1: 
        fig6.tight_layout()
        fig7.tight_layout()
        fig8.tight_layout()
        fig9.tight_layout()
        fig10.tight_layout()
        fig11.tight_layout()
        fig12.tight_layout()

        fig6.savefig(pathname + '/Δpn.png')
        fig7.savefig(pathname + '/Δp0.png')
        fig8.savefig(pathname + '/ΔSan.png')
        fig9.savefig(pathname + '/ΔSa0.png')
        fig10.savefig(pathname + '/ΔcTn.png')
        fig11.savefig(pathname + '/ΔcT0.png')
        fig12.savefig(pathname + '/alphabeta.png')

figblwn = plt.figure()
axblwn = figblwn.add_subplot(1,1,1)
axblwn.set_facecolor('k')
plot = plt.pcolor(lon,lat,np.transpose(belown), cmap = clmaperror,vmin = 0, vmax = 1)
plt.xlabel('Longitude / º')
plt.ylabel('Latitude / º')
figblwn.suptitle('Error codes BELOW with Rx = {}, Ry = {} \n for initial presure = {}\n Neutral calculation'.format(Rx,Ry,PRES[kAmax,0,0]))
clb.set_ticks([1/6,3/6,5/6])
clb.set_ticklabels(['No data', 'N2 neg', 'Crit neg'])
axblwn.set(xlim=(xmin,xmax))
axblwn.set(ylim=(ymin,ymax))

figblwns = plt.figure()
axblwns = figblwns.add_subplot(1,1,1)
axblwns.set_facecolor('k')
plot = plt.pcolor(lon,lat,np.transpose(belown - errorsn), cmap = clmaperror,vmin = 0, vmax = 1)
plt.xlabel('Longitude / º')
plt.ylabel('Latitude / º')
figblwns.suptitle('Error SUBSTRACTION BELOW with Rx = {}, Ry = {} \n for initial presure = {}\n Neutral calculation'.format(Rx,Ry,PRES[kAmax,0,0]))
clb.set_ticks([1/6,3/6,5/6])
clb.set_ticklabels(['No data', 'N2 neg', 'Crit neg'])
axblwns.set(xlim=(xmin,xmax))
axblwns.set(ylim=(ymin,ymax))

#MAP WITHOUT NORMALICING
subtractm = -Hghts0+Hghtsn
#we center the colormap in 0
shifted_cmap = fun.shiftedColorMap(matplotlib.cm.bwr, midpoint = 1 - np.nanmax(subtractm)/(np.nanmax(subtractm)+ abs(np.nanmin(subtractm))))
fignn = plt.figure()
axnn = fignn.add_subplot(1,1,1)
axnn.set_facecolor('k')
plt.contourf(x,y,SALI, colors ='w')
plot = plt.pcolor(lon,lat,np.transpose(-Hghts0+Hghtsn), cmap=shifted_cmap)
plt.xlabel('Longitude / º')
plt.ylabel('Latitude / º')
fignn.suptitle('Substraction of neutral minus sigma0 calculation for \n {} neutral buoyancy trajectories. Rx = {} º, Ry = {} º\n for initial height = {} m. Average = {} m'.format(shape,Rx,Ry,PRES[kAmax,0,0],avgs))
clb = plt.colorbar(plot)
clb.set_label('Pitch of trajectory / m ')
axnn.set(xlim=(xmin,xmax))
axnn.set(ylim=(ymin,ymax))

if folder == 1: 
    fignn.tight_layout()
    fignn.savefig(pathname + '/notnormalizedsub.png')

plt.show(block = 0)

input()

plt.close('all')

quit()
