"""
Simple code to read and plot WOCE data+simple neutral trajectory calculation
gsw = Gibbs Seawater toolbox in Python
"""
from math import cos, sin
import gsw
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import Functions as fun

# Directories
data_dir = r'C:\Users\Juan\Desktop\Programas Python' #need the r at the front otherwise struggles with \
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


#Start defining the loop

iter = 200
st = 2*3.14/iter
lat0 = -30
lon0 = 240
i0,j0 = fun.latlong(lat0,lon0,1)
R = 8
i = 0
lat = np.zeros(iter+1)
lon = np.zeros(iter+1)
i1 = np.zeros(iter+1)
j1 = np.zeros(iter+1)
lat1 = [lat0+R]
lon1 = [lon0]
i2, j2 = fun.latlong(lat0+R,lon0,0)

kA = 15 #initial depth index

while i<=iter:
    #Parametric equation of a circle
    lat[i] = lat0 + R*cos(st*i)
    lon[i] = lon0 + R*sin(st*i)
    #Get the index in the vector
    i1[i],j1[i] = fun.latlong(lat[i],lon[i],1)
    if i != 0 :
        #Check we are not repeating the same cast
        if i1[i]!=i1[i-1] or j1[i]!=j1[i-1]:
            lat1.append(lat[i])
            lon1.append(lon[i])
            i2.append(int(i1[i]))
            j2.append(int(j1[i]))
    i += 1


####################################################################################### Plots

prfalse, cT, sA = fun.trajectorie(i2,j2,i0,j0,PRES,SIG0,GAMN,SALI,TEMP,x,y,kA,False)
prtrue, cT, sA = fun.trajectorie(i2,j2,i0,j0,PRES,SIG0,GAMN,SALI,TEMP,x,y,kA,True) #True = alpha and beta computed for p = 0

prf_aux = prfalse
prt_aux = prtrue

#define auxiliar values of pressures list, so they are postive and we can calculate alpha/beta with pressure > 0

alpha_beta_f = gsw.alpha_on_beta (sA, cT, prf_aux)
alpha_beta_t = gsw.alpha_on_beta(sA, cT, prt_aux)

#If we ran into a problem while computing the trajectorie we stop the program
if len(prfalse) == 2:
    plt.figure()
    plt.pcolor(x[:],y[:],SALI[prfalse[1],:,:])
    plt.plot(x[prfalse[0]],y[prfalse[0]],'ko')
    plt.show()
    quit()


plt.figure()
plt.pcolor(x[:],y[:],SALI[0,:,:])


i3max, j3max = fun.latlong(lat0+1.5*R,lon0+1.5*R,1)
i3min,j3min = fun.latlong(lat0-1.5*R,lon0-1.5*R,1)

i3 = np.arange(i3min, i3max)
j3 = np.arange(j3min,j3max)
"""""
i = 0
j = 0
k = 0
l = 0
t = 0
level = 27
superficie = np.zeros((len(i3),len(j3)))
while i < len(i3):
    while j < len(j3):
        while k < 44:
            if GAMN[k,j3[j],i3[i]] <= level+0.1 and GAMN[k,j3[j],i3[i]] >= level-0.1:
                if superficie[i,j]  != None:
                    print("trouble")
                superficie[i,j] = -z[k]
            k += 1
        k = 0
        j += 1
    j = 0
    i += 1
print(*superficie[:,:])
"""

#colorsT = (cT - np.amin(cT))/np.amax(cT)
#viridis = plt.cm.jet(colorsT)

X, Y = np.meshgrid(x[i3], y[j3])

plt.figure()
ax = plt.axes(projection='3d')
ax.plot(x[i2],y[j2],prfalse)
ax.plot(x[i2],y[j2],prtrue, color='r', linewidth = 3, linestyle = '--')
ax.plot(x[i2[0]],y[j2[0]],prfalse[0],markerfacecolor='y', markeredgecolor='k', marker='*', markersize=10, alpha=1)

#ax.contour3D(X,Y,GAMN[kA,i3,j3],100)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
#ax.plot_surface(X,Y,superficie)
prdif = prfalse[0]-prfalse[len(prfalse)-1]
print("the difference in height is (for bool = false ):",prdif)
prdif = prtrue[0]-prtrue[len(prtrue)-1]
print("the difference in height is (for bool = true ):",prdif)

ax.set(xlim=(lon0-1.5*R,lon0+1.5*R))
ax.set(ylim=(lat0-1.5*R,lat0+1.5*R))


"""""
ax.contourf3D(x[i3],y[j3],-z[kA]+GAMN[kA,j3,i3],100,alpha = 0.2)
ax.contourf3D(x[i3],y[j3],-z[kA+1]+GAMN[kA,j3,i3],100,alpha = 0.2)
ax.contourf3D(x[i3],y[j3],-z[kA-1]+GAMN[kA,j3,i3],100,alpha = 0.2)
ax.contourf3D(x[i3],y[j3],-z[kA+2]+GAMN[kA,j3,i3],100,alpha = 0.2)
ax.contourf3D(x[i3],y[j3],-z[kA+3]+GAMN[kA,j3,i3],100,alpha = 0.2)

"""
plt.figure()
ax = plt.axes(projection='3d')
ax.plot(cT,sA,prfalse)
ax.plot(cT,sA,prtrue, color='r', linewidth = 3, linestyle = '--')
ax.plot(cT[0],sA[0],prfalse[0],markerfacecolor='y', markeredgecolor='k', marker='*', markersize=10, alpha=1)
#ax.contour3D(X,Y,GAMN[kA,i3,j3],100)
plt.xlabel('Conservative Temperature')
plt.ylabel('Absolute Salinity')
plt.title('Trajectory in thermodynamic space for a point in latitude '+str(lat0)+', longitude '+str(lon0)+' and depth '+str(prtrue[0]))
#ax.plot_surface(X,Y,superficie)
prdif = prfalse[0]-prfalse[len(prfalse)-1]
print("the difference in height is (for bool = false ):",prdif)
prdif = prtrue[0]-prtrue[len(prtrue)-1]
print("the difference in height is (for bool = true ):",prdif)

plt.figure()
plt.plot(alpha_beta_f, prfalse)
plt.xlabel('Alpha/Beta')
plt.ylabel('Pressure')
plt.title('Alpha/Beta from Bool = false for a point in latitude '+str(lat0)+', longitude '+str(lon0)+' and depth '+str(prtrue[0]))


plt.figure()
plt.plot(alpha_beta_t, prtrue)
plt.xlabel('Alpha/Beta')
plt.ylabel('Pressure')
plt.title('Alpha/Beta from Bool = true for a point in latitude '+str(lat0)+', longitude '+str(lon0)+' and depth '+str(prtrue[0]))
plt.show()
quit()
