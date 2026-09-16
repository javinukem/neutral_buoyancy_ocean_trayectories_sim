import math
import numpy as np
from cmath import nan
from math import cos, sin
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import Functions as fun
from matplotlib import cm
import os
import gsw
import wghc_subroutines as wghc
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def latlong(x,y,latitud, longitud, bool):
    # Returns the indices in the x and y arrays which correspond to the latitud and longitud
    # The bool argument controls wheter the indices are returned as a set or as a set of lists each containing the index

    if longitud < 0 or longitud > 360:
        longitud = longitud - np.sign(longitud)*360
    if latitud < -90 or latitud > 90:
        latitud = latitud - np.sign(latitud) * 180
    xdif = abs(x - longitud) #we define the pressure differences b/w our new cast and the measured value
    xdif_list = list(xdif) #change the result to a list element so we can later get an int index

    i0 = xdif_list.index(min(xdif_list))  #we find the index in the new cast for which the pressure difference w/ the measured value is minimum.

    ydif = abs(y - latitud) #we define the pressure differences b/w our new cast and the measured value
    ydif_list = list(ydif) #change the result to a list element so we can later get an int index

    j0 = ydif_list.index(min(ydif_list))

    #j0 = int(np.round(latitud*2+160))
    #i0 = int(np.round(longitud*2))

    if bool == 0:
        return([i0] ,[j0])
    else:
        return(i0 ,j0)

def trajectorie(i2,j2,PRES,SIG0,GAMN,SALI,TEMP,x,y,kA,Bool_sig0,Nabv):
    flag = 0
    i = 0
    int(kA)
    pr0 = PRES[:,j2[0],i2[0]]  # pressure
    #sg0 = SIG0[:,j2[0],i2[0]]  # potential density ref to 0 dbar
    #gm0 = GAMN[:,j2[0],i2[0]]  # approximate neutral density
    sp = SALI[:,j2[0],i2[0]]  # in situ (practical) salinity
    t = TEMP[:,j2[0],i2[0]]   # in situ temp
    sa0 = gsw.SA_from_SP(sp,pr0,x[i2[0]],y[j2[0]]) #"absolute" salinity
    ct0 = gsw.CT_from_t(sa0,t,pr0)  #conservative temp
    #rho0 = gsw.rho(sa0,ct0,pr0) # in situ density
    #prA = pr0[kA]
    #ctA = ct0[kA]
    ctA = np.zeros(len(i2))
    saA = np.zeros(len(i2))
    prA = np.zeros(len(i2))
    ctA[0] = ct0[kA]
    prA[0] = -pr0[kA]
    saA[0] = sa0[kA]
    i = 1

    while i <= len(i2)-1: #iterate to obtain another 10 different values of our trajectory

        if kA < Nabv:
            pr1 = PRES[0:2*Nabv,j2[i],i2[i]]
            #sg = SIG0[0:5,j2[i],i2[i]]  # potential density ref to 0 dbar
            #gm = GAMN[0:5,j2[i],i2[i]]  # approximate neutral density
            sp = SALI[0:2*Nabv,j2[i],i2[i]]  # in situ (practical) salinity
            t = TEMP[0:2*Nabv,j2[i],i2[i]]   # in situ temp
        else:
            pr1 = PRES[kA-Nabv:kA+Nabv+1,j2[i],i2[i]] # pressure levels centred on kA (Nabv points above, Nabv points below)
            #sg = SIG0[kA-Nabv:kA+Nabv+1,j2[i],i2[i]]  # potential density ref to 0 dbar
            #gm = GAMN[kA-Nabv:kA+Nabv+1,j2[i],i2[i]]  # approximate neutral density
            sp = SALI[kA-Nabv:kA+Nabv+1,j2[i],i2[i]]  # in situ (practical) salinity
            t = TEMP[kA-Nabv:kA+Nabv+1,j2[i],i2[i]]  # in situ temp
        sa1 = gsw.SA_from_SP(sp,pr1,x[i2[i]],y[j2[i]]) #"absolute" salinity
        ct1 = gsw.CT_from_t(sa1,t,pr1)   #conservative temp
        #compute location of B on a neighbouring cast (=1)


        prB, dct, dsa, plcholder, calcflag = wghc.calc_gamman_pressureWGHC(ctA[i-1],saA[i-1],1,-prA[i-1],ct1,sa1,2*ct1,pr1,Bool_sig0,0)

        if math.isnan(prB) == True:
             #prB, dct, dsa, plcholder, calcflag = wghc.calc_gamman_pressureWGHC(ctA[i-1],saA[i-1],1,-prA[i-1],ct1,sa1,2*ct1,pr1,Bool_sig0,1)
             flag = 1
             #print("Out of range in latitude:",i2[i],"and longtiude:", j2[i] ,"and height: ",prA[i-1] )
             return(prA, ctA, saA, flag ,calcflag)


        #now we redefine the thermodynamic values for our new cast

        ctA[i] = ctA[i-1] + dct
        saA[i] = saA[i-1] + dsa
        prA[i] = -prB
        prDif_list = list(abs(PRES[:, j2[i], i2[i]] - prB)) #change the result to a list element so we can later get an int index

        kA = prDif_list.index(min(prDif_list))  #we find the index in the new cast for which the pressure difference w/ the measured value is minimum.




        #print(ctA[i]-ct0[kA])
        #saA = sa0[kA]
        #ctA = ct0[kA]
        #prA = pr0[kA] ##So, we can redefine the values for our new cast with the new pressure index
        i += 1
    return(prA,ctA,saA, flag, calcflag)

def check(i2,j2,PRES,BOTD,kA,Nabv,pr):
    flag = 0

    if len(pr) == 1:
        if np.any(PRES[kA-Nabv:kA+Nabv+1,j2[len(pr)-1],i2[len(pr)-1]] > BOTD[j2[len(pr)-1],i2[len(pr)-1]]):
            flag = 1
    elif np.any(PRES[kA-Nabv:kA+Nabv+1,j2[len(pr)],i2[len(pr)]] > BOTD[j2[len(pr)],i2[len(pr)]]) :
        flag = 1
    elif PRES[kA-Nabv:kA+Nabv+1,j2[len(pr)],i2[len(pr)]].size == 0:
        flag = 1

    return(flag)

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def map(intx:int,inty:int,Rx: float,Ry:float,kAmax: int,wfile: bool,Nabv: int,ellipse: bool,square: bool,iter: int,graphs: bool):
    int(kAmax)
    int(Nabv)

    # Directories
    data_dir = 'C:/Users/Javier/Documents/Proyecto' #need the r at the front otherwise errorcodeuggles with \
    data_fileT = '/wghc_params.nc'
    dataT = data_dir + data_fileT

    # Constants
    km=1000
    day = 24*3600


    ############################################################################################### Read the file
    dataset = Dataset(dataT)

    x : np.ndarray = dataset.variables['LON'][:]  #longitude
    y : np.ndarray = dataset.variables['LAT'][:]  #latitude

    SIG0 = dataset.variables['SIG0']   #potential density referenced to 0 dbar [depth, latitude,longitude]
    GAMN = dataset.variables['GAMMAN'] #approximate neutral surface
    PRES = dataset.variables['PRES']   #hydrostatic pressure
    TEMP = dataset.variables['TEMP'] #in-situ, in deg C
    SALI = dataset.variables['SALINITY'] # in psu
    BOTD = dataset.variables['BOT_DEP'] #bottom depth (Nz,Ny,Nx)
    BOTD = BOTD[0,:,:]  #remove redundant data

    #CREATE A PATH TO SAVE FIGURES AND TXT

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
                    errorcode =  'Longitude {} and latitude {} *** ERROR BOTH*** Neutral : {} Sigma0 : {} '.format(lon0,lat0,calcflagn,calcflag0)
                    if np.any(PRES[kAn-Nabv:kAn+Nabv+1,j2[len(prn)],i2[len(prn)]] > BOTD[j2[len(prn)],i2[len(prn)]]):
                        errorcode = errorcode + bn
                        if calcflagn == 'no_data':
                            belown[i,j] = 1/6
                        elif calcflagn == 'ok':
                            belown[i,j] = 3/6
                        else:
                            belown[i,j] = 5/6
                    if np.any(PRES[kA0-Nabv:kA0+Nabv+1,j2[len(pr0)],i2[len(pr0)]] > BOTD[j2[len(pr0)],i2[len(pr0)]]):
                        errorcode = errorcode + b0
                        if calcflag0 == 'no_data':
                            below0[i,j] = 1/6
                        elif calcflag0 == 'ok':
                            below0[i,j] = 3/6
                        else:
                            below0[i,j] = 5/6
                    errorcode = errorcode + '\n'
                    file.write(errorcode)
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
                    errorcode =  'Longitude {} and latitude {} *** ERROR NTRL*** Neutral : {} '.format(lon0,lat0,calcflagn)
                    prDif_listn = list(abs(PRES[:, j2[len(prn)], i2[len(prn)]] + prn[-1]))
                    kAn = prDif_listn.index(min(prDif_listn))
                    if np.any(PRES[kAn-Nabv:kAn+Nabv+1,j2[len(prn)],i2[len(prn)]] > BOTD[j2[len(prn)],i2[len(prn)]]):
                        errorcode = errorcode + bn
                        if calcflagn == 'no_data':
                            belown[i,j] = 1/6
                        elif calcflagn == 'ok':
                            belown[i,j] = 3/6
                        else:
                            belown[i,j] = 5/6
                    errorcode = errorcode + '\n'
                    file.write(errorcode)
                print(i+1, 'out of ', imax,'and', j+1, 'out of ', jmax, '*** ERROR NTRL*** Neutral : ', calcflagn)
                if calcflagn == 'no_data':
                    errorsn[i,j] = 1/6
                elif calcflagn == 'ok':
                    errorsn[i,j] = 3/6
                else:
                    errorsn[i,j] = 5/6
                Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = nan, nan, nan, nan
                Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = pr0[len(pr0)-1]-pr0[0], abs(np.max(pr0)-np.min(pr0)), abs(np.max(sa0)-np.min(sa0)), abs(np.max(cT0)-np.min(cT0))
            elif flag0 == 1:
                if wfile == 1:
                    pr0 = pr0[np.where(pr0 != 0)]
                    print(i+1, 'out of ', imax,'and', j+1, 'out of ', jmax, '*** ERROR SGM0*** Sigma0 : ', calcflag0)
                    prDif_list0 = list(abs(PRES[:, j2[len(pr0)], i2[len(pr0)]] + pr0[-1]))
                    kA0 = prDif_list0.index(min(prDif_list0))
                    errorcode =  'Longitude {} and latitude {} *** ERROR SGM0 *** Sigma0 : {} '.format(lon0,lat0,calcflagn)
                    if np.any(PRES[kA0-Nabv:kA0+Nabv+1,j2[len(pr0)],i2[len(pr0)]] > BOTD[j2[len(pr0)],i2[len(pr0)]]):
                        errorcode = errorcode + b0
                        if calcflag0 == 'no_data':
                            below0[i,j] = 1/6
                        elif calcflag0 == 'ok':
                            below0[i,j] = 3/6
                        else:
                            below0[i,j] = 5/6
                    errorcode = errorcode + '\n'
                    file.write(errorcode)
                if calcflag0 == 'no_data':
                    errors0[i,j] = 1/6
                elif calcflag0 == 'ok':
                    errors0[i,j] = 3/6
                else:
                    errors0[i,j] = 5/6
                Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = prn[len(prn)-1]-prn[0], abs(np.max(prn)-np.min(prn)), abs(np.max(san)-np.min(san)), abs(np.max(cTn)-np.min(cTn))
                Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = nan, nan,nan,nan
            else:
                print(i+1, 'out of ', imax,'and', j+1, 'out of ', jmax)
                Hghtsn[i,j], Δpn[i,j], Δsan[i,j], ΔcTn[i,j] = prn[len(prn)-1]-prn[0], abs(np.max(prn)-np.min(prn)), abs(np.max(san)-np.min(san)), abs(np.max(cTn)-np.min(cTn))
                Hghts0[i,j], Δp0[i,j], Δsa0[i,j], ΔcT0[i,j] = pr0[len(pr0)-1]-pr0[0], abs(np.max(pr0)-np.min(pr0)), abs(np.max(sa0)-np.min(sa0)), abs(np.max(cT0)-np.min(cT0))
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
    avgf, avgt, avgs = np.nanmean(Hghtsn), np.nanmean(Hghts0), np.nanmean(Hghts0-Hghtsn)

    Hghts = np.zeros((len(x),len(y)))
    while i < len(x):
        while j < len(y):
            Hghts[i,j] = gsw.alpha_on_beta (SALI[kAmax,j,i], TEMP[kAmax,j,i], PRES[kAmax,j,i])
            j += 1
        print(i + 1, 'out of ', len(x))
        i += 1
        j = 0

    SALI = SALI[0,:,:] > 0

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(1,1,1)
    ax1.set_facecolor('k')
    plt.contourf(x,y,SALI, colors ='w')
    plot = plt.pcolor(lon,lat,np.transpose(Hghtsn), cmap='bwr', vmin = -10, vmax = 10)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig1.suptitle('Pitch of ellipses with Rx = {}, Ry = {} \n  for initial presure = {}. Average = {} \n Neutral calculation'.format(Rx,Ry,PRES[kAmax,0,0],avgf))
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
    fig2.suptitle('Pitch of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Average = {} \n Sigma0 calculation'.format(Rx,Ry,PRES[kAmax,0,0],avgt))
    clb = plt.colorbar(plot)
    clb.set_label('Pitch of trajectory')
    ax2.set(xlim=(xmin,xmax))
    ax2.set(ylim=(ymin,ymax))

    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot(1,1,1)
    ax3.set_facecolor('k')
    clmaperror = cm.get_cmap('brg', 3)
    plt.contourf(x,y,SALI, colors ='w')
    plot = plt.pcolor(lon,lat,np.transpose(errorsn), cmap = clmaperror,vmin = 0, vmax = 1)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig3.suptitle('Error codes NEUTRAL')
    clb = plt.colorbar(plot)
    clb.set_ticks([1/6,3/6,5/6])
    clb.set_ticklabels(['No data', 'N2 neg', 'Crit neg'])
    ax3.set(xlim=(xmin,xmax))
    ax3.set(ylim=(ymin,ymax))

    fig4 = plt.figure(4)
    ax4 = fig4.add_subplot(1,1,1)
    ax4.set_facecolor('k')
    clmaperror = cm.get_cmap('brg', 3)
    plt.contourf(x,y,SALI, colors ='w')
    plot = plt.pcolor(lon,lat,np.transpose(errors0), cmap = clmaperror,vmin = 0, vmax = 1)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig4.suptitle('Error codes SIGMA0')
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
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig5.suptitle('Difference of pitch of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Average = {} \n Sigma0 calculation'.format(Rx,Ry,PRES[kAmax,0,0],avgs))
    clb = plt.colorbar(plot)
    clb.set_label('Pitch of trajectory')
    ax5.set(xlim=(xmin,xmax))
    ax5.set(ylim=(ymin,ymax))

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()

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
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        fig6.suptitle("Range of pressure of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δp = {} \n Neutral calculation".format(Rx,Ry,PRES[kAmax,0,0],np.nanmax(Δpn)))
        clb = plt.colorbar(plot)
        clb.set_label('Range of pressure ')

        fig7 = plt.figure(7)
        ax7 = fig7.add_subplot(1,1,1)
        ax7.set_facecolor('k')
        ax7.set(xlim=(xmin,xmax))
        ax7.set(ylim=(ymin,ymax))
        plt.contourf(x,y,SALI, colors ='w')
        plot = plt.pcolor(lon,lat,np.transpose(Δp0), cmap='Blues')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        fig7.suptitle("Range of pressure of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δp  = {} \n Sigma 0 calculation".format(Rx,Ry,PRES[kAmax,0,0],np.nanmax(Δp0)))
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
        fig8.suptitle("Range of salinity of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δsa  = {} \n Neutral calculation".format(Rx,Ry,PRES[kAmax,0,0],np.nanmax(Δsan)))
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
        fig9.suptitle("Range of salinity of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max Δsa  = {} \n Sigma 0 calculation".format(Rx,Ry,PRES[kAmax,0,0],np.nanmax(Δsa0)))
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
        fig10.suptitle("Range of conserved temperature of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max ΔCT  = {} \n Neutral calculation".format(Rx,Ry,PRES[kAmax,0,0],np.nanmax(ΔcTn)))
        clb = plt.colorbar(plot)
        clb.set_label('Range of conserved temperature ')

        fig11 = plt.figure(11)
        ax11 = fig11.add_subplot(1,1,1)
        ax11.set_facecolor('k')
        ax11.set(xlim=(xmin,xmax))
        ax11.set(ylim=(ymin,ymax))
        plt.contourf(x,y,SALI, colors ='w')
        plot = plt.pcolor(lon,lat,np.transpose(ΔcT0), cmap='Blues')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        fig11.suptitle("Range of conserved temperature of ellipses with Rx = {}, Ry = {} \n for initial presure = {}. Max ΔCT  = {} \n Sigma 0 calculation".format(Rx,Ry,PRES[kAmax,0,0],np.nanmax(ΔcT0)))
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


    return()
