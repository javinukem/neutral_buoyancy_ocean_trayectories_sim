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
from matplotlib import animation
import numpy.ma as ma
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
BOTD = dataset.variables['BOT_DEP'][:] #bottom depth (Nz,Ny,Nx)
BOTD = BOTD[0,:,:]  #remove redundant data


#Start defining the loop

Nabv = 3
iter = 140
lon0 = 120
lat0 = 5
Rx = 15
Ry = 25
kA = 20
ellipse = 0
square = 1

i = 0
graphs = 0
scale = 1.2
scalesmall = 00.1

#define intial values

i0,j0 = fun.latlong(x,y,lat0,lon0,1)
lat = np.zeros(iter+1)
lon = np.zeros(iter+1)
i1 = np.zeros(iter+1)
j1 = np.zeros(iter+1)
lat1 = [lat0+Ry]
lon1 = [lon0]
i2, j2 = fun.latlong(x,y,lat0+Ry,lon0,0)

#compute the points of the ellipse
if ellipse == 1:
    lat = np.zeros(iter+1)
    lon = np.zeros(iter+1)
    i1 = np.zeros(iter+1)
    j1 = np.zeros(iter+1)
    lat1 = [lat0+Ry]
    lon1 = [lon0]
    i2, j2 = fun.latlong(x,y,lat0+Ry,lon0,0)
    st = 2*3.14/iter
    while i<=iter:
        lat[i] = lat0 + Ry*cos(st*i)
        lon[i] = lon0 + Rx*sin(st*i)
        #Get the index in the vector
        i1[i],j1[i] = fun.latlong(x,y,lat[i],lon[i],1)
        if i != 0 :
            #Check we are not repeating the same cast
            if i1[i]!=i1[i-1] or j1[i]!=j1[i-1]:
                lat1.append(lat[i]) 
                lon1.append(lon[i]) 
                i2.append(int(i1[i]))
                j2.append(int(j1[i]))
        i += 1    

if square ==1 :

    iter = int(4*(Rx+Ry))
    lat = np.zeros(iter+1)
    lon = np.zeros(iter+1)
    i1 = np.zeros(iter+1)
    j1 = np.zeros(iter+1)
    lat1 = [lat0+Ry/2]
    lon1 = [lon0-Rx/2]
    i2, j2 = fun.latlong(x,y,lat0+Ry/2,lon0-Rx/2,0)
    while i <= iter:
        if i <= 2*Rx:
            lat[i] = lat0 + Ry/2
            lon[i] = lon0 - Rx/2 + i*0.5
        elif i <= 2*(Rx+Ry):
            lat[i] = lat0 + Ry/2 - (i-2*Rx)*0.5
            lon[i] = lon0 + Rx/2 
        elif i <= 2*(2*Rx+Ry):
            lat[i] = lat0 - Ry/2
            lon[i] = lon0 + Rx/2 - (i-2*(Rx+Ry))*0.5
        else:
            lat[i] = lat0 - Ry/2 + (i-2*(2*Rx+Ry))*0.5
            lon[i] = lon0 - Rx/2 

        i1[i],j1[i] = fun.latlong(x,y,lat[i],lon[i],1)

        if i != 0 :
            #Check we are not repeating the same cast
            if i1[i]!=i1[i-1] or j1[i]!=j1[i-1]:
                lat1.append(lat[i]) 
                lon1.append(lon[i]) 
                i2.append(int(i1[i]))
                j2.append(int(j1[i]))
        i += 1    



figxd = plt.figure()
plt.plot(lon,lat)
#plt.plot(lon1,lat1,'ko')
plt.plot(x[i2],y[j2],'ro')
plt.ylabel('Latitude / º')
plt.xlabel('Longitude / º /g/Kg')

#computes the neutral and sigma 0 trajectories
prn, cTn, sAn, flagn, calcflagn = fun.trajectorie(i2,j2,PRES,SIG0,GAMN,SALI,TEMP,x,y,kA,False,Nabv)

if flagn == 1:

    print(calcflagn)
    prn = prn[np.where(prn != 0)]

    i2red = i2[0:len(prn)]
    j2red = j2[0:len(prn)]
    plt.plot(x[i2[len(prn)]],y[j2[len(prn)]],'g*')

    #plt.plot(y[j2[len(prn)+1]],x[i2[len(prn)+1]],'ro')


    prDif_list = list(abs(PRES[:, j2[len(prn)], i2[len(prn)]] + prn[-1]))    
    kA = prDif_list.index(min(prDif_list))  
    if len(prn) == 1:
        if np.any(PRES[kA-Nabv:kA+Nabv+1,j2[len(prn)-1],i2[len(prn)-1]] > BOTD[j2[len(prn)-1],i2[len(prn)-1]]):
            print('xD LEN1')
        if np.any(np.isnan(PRES[kA-Nabv:kA+Nabv+1,j2[len(prn)-1],i2[len(prn)-1]])):
            print('xDDD LEN1')
        if np.any(np.isnan(TEMP[kA-Nabv:kA+Nabv+1,j2[len(prn)-1],i2[len(prn)-1]])):
            print('xDDDDD LEN1')
        if np.any(ma.is_masked(TEMP[kA-Nabv:kA+Nabv+1,j2[len(prn)-1],i2[len(prn)-1]])):
            print('xDD maskeaooooD LEN1')
        elif BOTD[j2[len(prn)-1],i2[len(prn)-1]] < 0:
            print('sale por arriba Uwu')
    if np.any(PRES[kA-Nabv:kA+Nabv+1,j2[len(prn)],i2[len(prn)]] > BOTD[j2[len(prn)],i2[len(prn)]]) :
        print('xD')
    if np.any(np.isnan(PRES[kA-Nabv:kA+Nabv+1,j2[len(prn)],i2[len(prn)]])):
        print('xDDD')
    if np.any(np.isnan(TEMP[kA-Nabv:kA+Nabv+1,j2[len(prn)],i2[len(prn)]])):
        print('xDDDDD')
    if np.any(ma.is_masked(TEMP[kA-Nabv:kA+Nabv+1,j2[len(prn)],i2[len(prn)]])):
        print('xDD maskeaooooD')
    elif BOTD[j2[len(prn)],i2[len(prn)]] < 0:
        print('sale por arriba Uwu')
    if PRES[kA-Nabv:kA+Nabv+1,j2[len(prn)],i2[len(prn)]].size == 0:
        print('EMPTY')
    i3max, j3max = fun.latlong(x,y,y[j2red[-1]]+scalesmall*Ry,x[i2red[-1]]+scalesmall*Ry,1)
    i3min,j3min = fun.latlong(x,y,y[j2red[-1]]-scalesmall*Ry,x[i2red[-1]]-scalesmall*Ry,1)
    i3 = np.arange(i3min, i3max)
    j3 = np.arange(j3min,j3max)
    X, Y = np.meshgrid(x[i3], y[j3])

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    cb = ax.scatter(x[i2red],y[j2red],prn, c = cTn[0:len(prn)])
    #ax.plot_surface(X,Y,-BOTD[np.ix_(j3,i3)],alpha = 0.5)
    fig.colorbar(cb, shrink=0.5, aspect=5)
    ax.plot(x[i2[0]],y[j2[0]],prn[0],markerfacecolor='y', markeredgecolor='k', marker='*', markersize=10, alpha=1)
    ax.plot(x[i2red[-1]],y[j2red[-1]],prn[-1],markerfacecolor='y', markeredgecolor='k', marker='v', markersize=10, alpha=1)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    ax.set(xlim=(lon0-scale*Ry,lon0+scale*Ry))
    ax.set(ylim=(lat0-scale*Ry,lat0+scale*Ry))
    plt.title("Error at longitude: {} and latitude: {} ".format(x[i2red[-1]],y[j2red[-1]]) )

    fig2 = plt.figure()
    ax = plt.axes(projection='3d')
    cb = ax.scatter(x[i2red],y[j2red],prn, c = cTn[0:len(prn)])
    ax.plot_surface(X,Y,-BOTD[np.ix_(j3,i3)],alpha = 0.5)
    #clb = fig2.colorbar(cb, shrink=0.5, aspect=5)
    ax.plot(x[i2[0]],y[j2[0]],prn[0],markerfacecolor='y', markeredgecolor='k', marker='*', markersize=10, alpha=1)
    ax.plot(x[i2red[-1]],y[j2red[-1]],prn[-1],markerfacecolor='y', markeredgecolor='k', marker='v', markersize=10, alpha=1)
    ax.plot(x[i2[len(prn)]],y[j2[len(prn)]],prn[-1],markerfacecolor='r', markeredgecolor='k', marker='^', markersize=10, alpha=1)
    plt.xlabel('Longitude º')
    plt.ylabel('Latitude º')
    ax.set_zlabel('Height / m')
    ax.set(xlim=(x[i2red[-1]]-scalesmall*Ry,x[i2red[-1]]+scalesmall*Ry))
    ax.set(ylim=(y[j2red[-1]]-scalesmall*Ry,y[j2red[-1]]+scalesmall*Ry))
    plt.title("Error at longitude: {} and latitude: {} . Zoom in".format(x[i2red[-1]],y[j2red[-1]]) )
    #clb.set_label('Conserved Temperature / K') 

    fig3 = plt.figure()
    ax3 = plt.axes()
    cb = ax3.scatter(x[i2red],prn, c = cTn[0:len(prn)])
    ax3.plot(x[i3],-BOTD[j3[-1],i3])
    #clb = fig3.colorbar(cb, shrink=0.5, aspect=5)
    ax3.plot(x[i2[0]],prn[0],markerfacecolor='y', markeredgecolor='k', marker='*', markersize=10, alpha=1)
    ax3.plot(x[i2red[-1]],prn[-1],markerfacecolor='y', markeredgecolor='k', marker='v', markersize=10, alpha=1)
    plt.xlabel('Longitude')
    plt.ylabel('height')
    plt.title("Error at longitude: {} and latitude: {} . Zoom in".format(x[i2red[-1]],y[j2red[-1]]) )
    #clb.set_label('Conserved Temperature / K') 

    fig4 = plt.figure()
    ax4 = plt.axes()
    cb = ax4.scatter(y[j2red],prn, c = cTn[0:len(prn)])
    ax4.plot(y[j3],-BOTD[j3,i3[-1]])
    clb = fig4.colorbar(cb, shrink=0.5, aspect=5)
    ax4.plot(y[j2[0]],prn[0],markerfacecolor='y', markeredgecolor='k', marker='*', markersize=10, alpha=1)
    ax4.plot(y[j2red[-1]],prn[-1],markerfacecolor='y', markeredgecolor='k', marker='v', markersize=10, alpha=1)
    plt.xlabel('Latitude')
    plt.ylabel('height')
    plt.title("Error at longitude: {} and latitude: {} . Zoom in".format(x[i2red[-1]],y[j2red[-1]]) )
    clb.set_label('Conserved Temperature / K') 

    #print("I ran into a problem with bool = false in longitude : ", x[i2[prn[0]]], "and latitude", y[j2[prn[0]]])
    
    plt.show()
    quit()

pr0, cT0, sa0, flag0, calcflag0 = fun.trajectorie(i2,j2,PRES,SIG0,GAMN,SALI,TEMP,x,y,kA,True,Nabv) #True = alpha and beta computed for p = 0


#If we ran into a problem while computing the trajectorie we stop the program

if flag0 == 1:
    plt.figure()
    plt.pcolor(x[:],y[:],SALI[pr0[1],:,:])
    plt.plot(x[i2[pr0[0]]],y[j2[pr0[0]]],'ko')
    print("I ran into a problem with bool = true in longitude : ", x[i2[pr0[0]]], "and latitude", y[j2[pr0[0]]])
    plt.show()
    quit()

prdif = prn[0]-prn[len(prn)-1]
print("the difference in height is (for bool = false ):",prdif)
prdif = pr0[0]-pr0[len(pr0)-1]
print("the difference in height is (for bool = true ):",prdif)

#PLOTS
i3max, j3max = fun.latlong(x,y,lat0+scale*Ry,lon0+scale*Ry,1)
i3min,j3min = fun.latlong(x,y,lat0-scale*Ry,lon0-scale*Ry,1)
i3 = np.arange(i3min, i3max)
j3 = np.arange(j3min,j3max)
X, Y = np.meshgrid(x[i3], y[j3])

    
    
fig2 = plt.figure()
ax = fig2.add_subplot(111,projection='3d')
cb = ax.scatter(lon1,lat1,prn, c = cTn)
ax.plot(lon1[0],lat1[0],prn[0],markerfacecolor='y', markeredgecolor='k', marker='*', markersize=10, alpha=1)
ax.plot(lon1[-1],lat1[-1],prn[-1],markerfacecolor='y', markeredgecolor='k', marker='v', markersize=10, alpha=1)
plt.xlabel('Longitude º')
plt.ylabel('Latitude º')
ax.set(xlim=(lon0-scale*Ry,lon0+scale*Ry))
ax.set(ylim=(lat0-scale*Ry,lat0+scale*Ry))
ax.set_zlabel('Height / m')
clb = fig2.colorbar(cb, shrink=0.5, aspect=5)
clb.set_label('Conserved Temperature / K')
'''
def init():
    ax.view_init(elev=90., azim=0)
    return [ax]

def animate(i):
    ax.view_init(elev=90., azim=i)
    return [ax]

anim = animation.FuncAnimation(fig2, animate, init_func=init,frames=360,interval = 30)
anim.save('animationellipabove.gif')
'''
if graphs == 1:
        
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter(cTn,sAn,prn)
    ax.scatter(cT0,sa0,pr0, color='r', linewidth = 3, linestyle = '--')
    ax.plot(cTn[0],sAn[0],prn[0],markerfacecolor='y', markeredgecolor='k', marker='*', markersize=10, alpha=1)
    ax.plot(cT0[0],sa0[0],pr0[0],markerfacecolor='g', markeredgecolor='k', marker='*', markersize=10, alpha=1)
    plt.xlabel('Conservative Temperature')
    plt.ylabel('Absolute Salinity')
    plt.title('Trajectory in thermodynamic space')

    plt.figure()
    plt.plot(cTn, -prn)
    plt.xlabel('Temperature')
    plt.ylabel('Pressure')
    plt.title('Alpha/Beta from Bool = false for a point in latitude '+str(lat0)+', longitude '+str(lon0)+' and depth '+str(pr0[0]))

    alpha_beta_f = gsw.alpha_on_beta (sAn, cTn, -prn)
    alpha_beta_t = gsw.alpha_on_beta(sa0, cT0, np.zeros(len(pr0)))

    plt.figure()
    plt.plot(alpha_beta_f, prn)
    plt.xlabel('Alpha/Beta')
    plt.ylabel('Pressure')
    plt.title('Alpha/Beta from Bool = false for a point in latitude '+str(lat0)+', longitude '+str(lon0)+' and depth '+str(pr0[0]))

    plt.figure()
    plt.plot(alpha_beta_t, pr0)
    plt.xlabel('Alpha/Beta')
    plt.ylabel('Pressure')
    plt.title('Alpha/Beta from Bool = true for a point in latitude '+str(lat0)+', longitude '+str(lon0)+' and depth '+str(pr0[0]))

    i = 0
    ΔSan = np.zeros(len(sAn)-1)
    ΔcT = np.zeros(len(sAn)-1)
    checkn = np.zeros(len(sAn)-1)

    while i < len(sAn)-1:
        ΔSan[i] = sAn[i+1]-sAn[i]
        ΔcT[i] = cTn[i+1]-cTn[i]
        checkn[i] = alpha_beta_f[i]*(ΔcT[i]/ΔSan[i])
        i += 1
        
    plt.figure()

    plt.plot(checkn,prn[0:-1])
    plt.title('neutral check')
    
    i = 0
    ΔSa0 = np.zeros(len(sa0)-1)
    ΔcT0  = np.zeros(len(sa0)-1)
    check0 = np.zeros(len(sa0)-1)
    while i < len(sa0)-1:
        ΔSa0[i] = sa0[i+1]-sa0[i]
        ΔcT0[i] = cT0[i+1]-cT0[i]
        check0[i] = alpha_beta_t[i]*(ΔcT0[i]/ΔSa0[i])
        i += 1




    plt.figure()
    plt.plot(check0,pr0[0:-1])
    plt.title('sigma 0 check')


    plt.show(block = 0)

    input()

    plt.close('all')

plt.show()

quit()