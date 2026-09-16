"""
Subroutines useful for WHGC climato
(c) Arnaud, May 2021
"""

import numpy as np
import gsw
import matplotlib.pyplot as plt
import numpy.ma as ma


# Constants
km=1000
day = 24*3600
OM = 2*np.pi/day
RADIUS = 6371*km
gravi = 9.81
dbar = 1.e04

def isOCEAN(botd0):
    """
    =1 if ocean grid point
    =nan if not
    Based on the mask of the variable BOTD
    Input: value of BOTD at the gridpoint considered botd0 = BOTD[j,i]
    NB: includes marginal seas (e.g., Caspian Sea)
    """
    zbot0 = -botd0 #to convert into a depth
    if ma.getmask(zbot0)==True or zbot0>=0: #land grid points
        mask = np.nan
    else:
        mask = 1
    return mask


def isNA(lon,lat):
    # = 1 if in North Atlantic
    # = nan if not
    # lat in degree north (-90 to 90)
    # lon in degree east (0 to 360)
    # NATL = north of 20N and south of 70N
    # NB: does not separate between land and ocean so should run isOCEAN with it
    mask = np.nan
    if lat >= 20 and lon >= 260: #biggest chunk, including portion of Gulf of Mexico
        mask = 1
    if lat >= 50 and lon <=20: #high lat portion, excluding Barents Sea
        mask = 1
    if lat >=30 and lat <=45 and lon >= 351:
        mask = np.nan
    if lat > 70:
        mask = np.nan
    return mask

def isSO(lon,lat):
    # = 1 if in Southern Ocean
    # = nan if not
    # lat in degree north (-90 to 90)
    # lon in degree east (0 to 360)
    # SO = south of 30S
    # NB: does not separate between land and ocean so should run isOCEAN with it
    mask = np.nan
    if lat <= -30:
        mask = 1
    return mask

def isACC(lon,lat):
    # = 1 if in ACC band (cuurently 65S-40S)
    # = nan if not
    # lat in degree north (-90 to 90)
    # lon in degree east (0 to 360)
    # ACC = south of 40S
    # NB: does not separate between land and ocean so should run isOCEAN with it
    mask = np.nan
    if lat <= -40 and lat >=-65:
        mask = 1
    return mask

def isAR(lon,lat):
    # = 1 if in Arctic
    # = nan if not
    # lat in degree north (-90 to 90)
    # lon in degree east (0 to 360)
    # AR = north of 60N
    # NB: does not separate between land and ocean so should run isOCEAN with it
    mask = np.nan
    if lat >= 60:
        mask = 1
    return mask

def split_list(VARlist,lay_index):
    # VARlist contains elements [(0,VAR), (0,VAR), (1,VAR), (2, VAR), ...]
    # i.e. (index of layer, VAR at that location in layer)
    # Nlay is the overall number of layers
    # The subroutine returns the list corresponding to layer lay_index only
    # VARlist_lay = [VAR, ..., VAR]
    M = len(VARlist)
    VARlist_lay = []
    VARlist_sort = sorted(VARlist, key=lambda VARlist: VARlist[0])  #sort the list based on the 1st element of the "duet" (.,.) i.e. the layer index
    Nduet = len(VARlist_sort)
    m = 0
    A = VARlist_sort[m] #1st duet
    while A[0] != lay_index and m<= Nduet-2:
        m = m+1
        A = VARlist_sort[m] #next duet
    else:
        while A[0] == lay_index:
              VARlist_lay.append(A[1])
              if m==M-1:
                 break
              else:
                  m = m+1
                  A = VARlist_sort[m]
    return VARlist_lay


def volum_pdf(pv,vol,pvmin,pvmax,delta_pv):      
    """
    Computes the pdf of a variable (e.g. PV) along an isopycnal: fraction of total volume occupied by a given PV bin on an ispycnal
    pv = pv(lat,lon) or 1D array is the PV matrix in PVU
    vol = vol(lat,lon) or 1D array but same format as PV is the volume occupied by each cell
    pvmin, pvmax, delta_pv = define the resolution of the PDF
    nb: does not take into account data where pv<pvmin or pv>pvmax
    """
    volhist = []
    bin_centre = []
    i = pvmin
    bin_edges = [i,]
    while i < pvmax:
           j = i+delta_pv
           a = np.where(pv<j,pv,-100)
           b = np.where(a>=i,vol,np.nan)
           volhist.append(np.nansum(b)) #volume occupied by a given PV bin
           bin_centre.append((i+j)/2) #center of pv bin
           bin_edges.append(j)
           i = j
           
    volpdf = volhist/np.sum(volhist)
    return volhist, volpdf, bin_centre, bin_edges

def calc_OMgeo(phiA,zA,ds):
    # Computes the lat phiB=phiA+dphi (in degrees) and height zB=zA+dz (in m)
    # of a column AB parallel to the axis of rotation with length ds measured
    # along this axis
   
    dphi = (ds*180/np.pi)*np.cos(phiA*np.pi/180)/RADIUS
    dz = ds*np.sin(phiA*np.pi/180)
    return phiA+dphi, zA+dz

def calc_gamman_pressureWGHC(ct0i,sa0i,ox0i,p0i,ct1g,sa1g,ox1g,p1g,Bool_sig0,dummy):
    # Computes the depth of neutral surface, defined by the conservative temp and absolute salinity
    # at cast0 at p = pg0[ko], on the adjacent cast1 (p=p1).
    #
    # Input arrays are centred on ko, i.e., typically N = len(pg) =5,
    # corresponding to 2 points above ko and 2 points below ko. This means
    # that there is an assumption that the pressure grids p0g and p1g are not
    # too dissimilar and a test is performed to check that. Should not be a pbm for
    # gridded data like wghc.
    #
    # Interpolates cast1 profiles by a parabola and solve the resulting quadratic
    # equation for P1 (as suggested in MCDougall 1987). Computes an approximate value
    # of N2 for the two solutions to remove non physical one.
    #
    # if Bool_sig0 == True then computes sigma0 surface by computing alpha and beta with pref=0
    # if dummy = 1 then proceed to debugging the code
    #
    # Returns pressure of neutral surface at cast 1 as well as the ct, sa and ox increment  
    #
    # NB The variable ox can be anything. For example the depth z
   
    if dummy == 1:
        DEBUG = 1
    else:
        DEBUG = 0
    p1 = np.nan   #global p1    #seems required to avoid "local variable referenced before assignment" error message
    dct = np.nan
    dsa = np.nan
    dox = np.nan
    flag = 'ok'

    ################################# Values on neutral surface (cast0)
    #p01 = p1g[ko]
    #if np.abs(p0i-p01)>50: #threshold in dbar
    #    print('!! calc_gamman_pressureWGHC...warning: pressure chunk far from p0i !!')
    #    #quit()
    # flag
    if np.isnan(ct0i)== 1 or np.isnan(sa0i)==1:
        p1 = np.nan     #no data
        dct = np.nan
        dsa = np.nan
        dox = np.nan
        flag = 'no_data'
        return p1, dct, dsa, dox, flag
   
    # thermo coeffs
    if Bool_sig0 == True:
       alpha = gsw.alpha(sa0i,ct0i,0)
       beta = gsw.beta(sa0i,ct0i,0)
    else:
       alpha = gsw.alpha(sa0i,ct0i,p0i)
       beta = gsw.beta(sa0i,ct0i,p0i)    
    sg0i = gsw.sigma0(sa0i,ct0i)
    sg1g = gsw.sigma0(sa1g,ct1g)

    ################################# New profile (cast1)
    # fit profile 1
    polyct1 = np.polyfit(p1g-p0i,alpha*ct1g,2)
    polysa1 = np.polyfit(p1g-p0i,beta*sa1g,2)
    polyox1 = np.polyfit(p1g-p0i,ox1g,2)
    a = polyct1[0]
    b = polyct1[1]
    c = polyct1[2]
    d = polysa1[0]
    e = polysa1[1]
    f = polysa1[2]
    g = alpha*ct0i - beta*sa0i
    # parameters
    la = (b-e)/(a-d)
    mu = (c-f-g)/(a-d)
    crit = .25*la**2-mu
    # New pressure of neutral surface
    if crit < 0:
        p1 = np.nan     #no mathematical solution
        dct = np.nan
        dsa = np.nan
        dox = np.nan
        flag = 'crit_neg'
        if DEBUG==1:
           print('Neutral surf subroutine... CRIT is NEGATIVE!!')
           print('Neutral surf subroutine... Initial pressure: ', p0i)
           print('Neutral surf subroutine... Critical coeff: ', crit)
           print('Neutral surf subroutine... Coeff a and d: ', a, d)
           print('Neutral surf subroutine... Coeff lambda and mu: ', la, mu)
           print('')
           p1g_poly = np.arange(p1g[0],p1g[-1],5)
           ct1g_poly = np.polyval(polyct1,p1g_poly-p0i)/alpha
           sa1g_poly = np.polyval(polysa1,p1g_poly-p0i)/beta
           ox1g_poly = np.polyval(polyox1,p1g_poly-p0i)
           sg1g_poly = gsw.sigma0(sa1g_poly,ct1g_poly)
           plt.figure()
           plt.plot(ox1g,-p1g,'r-o',ox1g_poly,-p1g_poly,'r--',[ox0i, ox0i],[-p0i-50, -p0i+50],'k--')
           plt.title('Oxygen')
           plt.figure()
           plt.plot(sg1g,-p1g,'r-o',sg1g_poly,-p1g_poly,'r--',[sg0i, sg0i],[-p0i-50, -p0i+50],'k--')
           plt.title('Sigma0')
           plt.figure()
           plt.plot(ct1g,-p1g, 'b-o',ct1g_poly,-p1g_poly,'b--',ct0i,-p0i,'r*')
           plt.title('Conserved temperature')
           plt.figure()
           plt.title('Salinity')
           plt.plot(sa1g,-p1g, 'b-o',sa1g_poly,-p1g_poly,'b--',sa0i,-p0i,'r*')
           plt.figure()
           plt.title('Salinity ConservedT')
           plt.scatter(sa1g,ct1g)
           plt.plot(sa0i,ct0i,'r*')
           plt.show()
        return p1, dct, dsa, dox, flag
    else:
        root1 = p0i -.5*la + np.sqrt(crit)
        root2 = p0i -.5*la - np.sqrt(crit)
        N2prop1 = 2*(root1-p0i)*(d-a)+e-b  
        N2prop2 = 2*(root2-p0i)*(d-a)+e-b  
        if N2prop1>0 and N2prop2<0:    #use approximate N2 to eliminate non physical branch
            p1 = root1
        if N2prop1<0 and N2prop2>0:
            p1 = root2
        if N2prop1<0 and N2prop2<0:
            p1 = np.nan
            dct = np.nan
            dsa = np.nan
            dox = np.nan
            if DEBUG==1:
                print('N1 AND N2 NEGATIVE')
                print('Neutral surf subroutine... Initial pressure: ', p0i)
                print('Neutral surf subroutine... Critical coeff: ', crit)
                print('Neutral surf subroutine... Coeff a and d: ', a, d)
                print('Neutral surf subroutine... Coeff lambda and mu: ', la, mu)
                print('Neutral surf subroutine... root1 and root2: ',root1, root2)
                print('Neutral surf subroutine... N2(root1) and N2(root2): ',N2prop1, N2prop2)
                print('Neutral surf subroutine... New pressure: ', p1)
                print('Neutral surf subroutine... CT increment: ',dct)
                print('Neutral surf subroutine... SA increment: ',dsa)
                print('Neutral surf subroutine... OX increment: ',dox)
                print('Neutral surf subroutine... sigma0 increment: ',dsg)
                print('')
                p1g_poly = np.arange(p1g[0],p1g[-1],5)
                ct1g_poly = np.polyval(polyct1,p1g_poly-p0i)/alpha
                sa1g_poly = np.polyval(polysa1,p1g_poly-p0i)/beta
                ox1g_poly = np.polyval(polyox1,p1g_poly-p0i)
                sg1g_poly = gsw.sigma0(sa1g_poly,ct1g_poly)
                plt.figure(75)
                plt.plot(ox1g,-p1g,'r-o',ox1g_poly,-p1g_poly,'r--',[ox0i, ox0i],[-p0i-50, -p0i+50],'k--')
                plt.title('Oxygen')
                plt.figure(76)
                plt.plot(sg1g,-p1g,'r-o',sg1g_poly,-p1g_poly,'r--',[sg0i, sg0i],[-p0i-50, -p0i+50],'k--')
                plt.figure(77)
                plt.plot(ct1g,-p1g, 'b-o',ct1g_poly,-p1g_poly,'b--',ct0i,-p0i,'r*')
                plt.figure(78)
                plt.plot(sa1g,-p1g, 'b-o',sa1g_poly,-p1g_poly,'b--',sa0i,-p0i,'r*')
                plt.figure(79)
                plt.scatter(sa1g,ct1g)
                plt.plot(sa0i,ct0i,'r*')
                plt.show()
            return p1, dct, dsa, dox
        if N2prop1>0 and N2prop2>0:
            dist1 = np.abs(root1-p0i)
            dist2 = np.abs(root2-p0i)
            if dist1 < dist2:
                p1 = root1
            else:
                p1 = root2
    # sa, ox and ct at p1 on cast 1
    sa1f = np.polyval(polysa1,p1-p0i)/beta
    ct1f = np.polyval(polyct1,p1-p0i)/alpha
    ox1f = np.polyval(polyox1,p1-p0i)
    sg1f = gsw.sigma0(sa1f,ct1f)
    # change in ct and sa along neutral surface from cast 0 to cast 1
    dct = ct1f-ct0i
    dsa = sa1f-sa0i
    dox = ox1f-ox0i
    dsg = sg1f-sg0i
    # outputs & graphics (for testing/debugging)
    '''
    if DEBUG==1:
        print('')
        print('Neutral surf subroutine... Initial pressure: ', p0i)
        print('Neutral surf subroutine... Critical coeff: ', crit)
        print('Neutral surf subroutine... Coeff a and d: ', a, d)
        print('Neutral surf subroutine... Coeff lambda and mu: ', la, mu)
        print('Neutral surf subroutine... root1 and root2: ',root1, root2)
        print('Neutral surf subroutine... N2(root1) and N2(root2): ',N2prop1, N2prop2)
        print('Neutral surf subroutine... New pressure: ', p1)
        print('Neutral surf subroutine... CT increment: ',dct)
        print('Neutral surf subroutine... SA increment: ',dsa)
        print('Neutral surf subroutine... OX increment: ',dox)
        print('Neutral surf subroutine... sigma0 increment: ',dsg)
        print('')
        p1g_poly = np.arange(p1g[0],p1g[-1],5)
        ct1g_poly = np.polyval(polyct1,p1g_poly-p0i)/alpha
        sa1g_poly = np.polyval(polysa1,p1g_poly-p0i)/beta
        ox1g_poly = np.polyval(polyox1,p1g_poly-p0i)
        sg1g_poly = gsw.sigma0(sa1g_poly,ct1g_poly)
        plt.figure(75)
        plt.plot(ox1g,-p1g,'r-o',ox1g_poly,-p1g_poly,'r--',[ox0i, ox0i],[-p0i-50, -p0i+50],'k--')
        plt.title('Oxygen')
        plt.figure(76)
        plt.plot(sg1g,-p1g,'r-o',sg1g_poly,-p1g_poly,'r--',[sg0i, sg0i],[-p0i-50, -p0i+50],'k--')
        plt.figure(77)
        plt.plot(ct1g,-p1g, 'b-o',ct1g_poly,-p1g_poly,'b--',ct0i,-p0i,'r*')
        plt.figure(78)
        plt.plot(sa1g,-p1g, 'b-o',sa1g_poly,-p1g_poly,'b--',sa0i,-p0i,'r*')
        plt.figure(79)
        plt.scatter(sa1g,ct1g)
        plt.plot(sa0i,ct0i,'r*')
        plt.show()
        quit()
    '''
    return p1, dct, dsa, dox, flag   #depth of the neutral surface p1 on cst 1 and change in ct, sa and ox

def calc_f(lat):
    """ Computes Coriolis parameter at latitude lat (in degrees)
    """
    f = 2*OM*np.sin(lat*np.pi/180)
    return f



def calc_mvector(sa,ct,pr,zz,x,y,io,jo,ko):
    """Computes the vector normal to the local neutral plane both directly from the local vertical gradients
    (=mvec_dir) and by fitting a plane using 8 points surrounding a central point (=mvec_fit).
   
    ---> INPUTS: 9 "casts" in total forming a "cube":
         sa,ct,pr,zz all of dimensions (Nz,3,3)
         x and y are the lon and latitude 1D-arrays
         io,jo,ko are the indices of the central point
         The ordering is as follows:
              list (j,i): [central point, (0,0), (1,0), (2,0) ..., (0,1)]
              i.e. clockwise from sw corner of the reference point.
              Label points with ref as 0 and then clockwise, so that
              the labelling fits that of the lists x_n, y_n, etc.
 
                     3    4     5
                     2  ref=0   6
                     1    8     7

    --> OUTPUTS: mvec_fit, mvec_dir = the coordinates of the mvector in the local
    Cartesian coordinate system (mx, my, 1)

    """
    import sys
    FIG = 0
    dx = RADIUS*np.cos(y[jo]*np.pi/180)*0.5*np.pi/180 #distance between two gridpoints  
    dy = RADIUS*0.5*np.pi/180
    dx_perdeg = RADIUS*np.cos(y[jo]*np.pi/180)*np.pi/180 #distance per degree
    dy_perdeg = RADIUS*np.pi/180
   
    ####################################################################### 1. Direct calculation
    alpha = gsw.alpha(sa[ko,1,1],ct[ko,1,1],pr[ko,1,1])
    beta = gsw.beta(sa[ko,1,1],ct[ko,1,1],pr[ko,1,1])
    gradCTx = (ct[ko,1,2]-ct[ko,1,0])/(2*dx)
    gradCTy = (ct[ko,2,1]-ct[ko,0,1])/(2*dy)
    gradCTz = (ct[ko+1,1,1]-ct[ko-1,1,1])/(zz[ko+1,1,1]-zz[ko-1,1,1])
    gradCT = np.array([gradCTx,gradCTy,gradCTz]) #row vector
    gradSAx = (sa[ko,1,2]-sa[ko,1,0])/(2*dx)
    gradSAy = (sa[ko,2,1]-sa[ko,0,1])/(2*dy)
    gradSAz = (sa[ko+1,1,1]-sa[ko-1,1,1])/(zz[ko+1,1,1]-zz[ko-1,1,1])
    gradSA = np.array([gradSAx,gradSAy,gradSAz])
    N2ct = gravi*alpha*gradCTz
    N2sa = -gravi*beta*gradSAz
    N2 = N2ct+N2sa
    mvec_dir = (alpha*gradCT-beta*gradSA)*gravi/N2  #z-component = 1 as in McDougall88

    ######################################################################### 2. Fitting the neutral plane
    # reference point
    ct0i = ct[ko,1,1]
    sa0i = sa[ko,1,1]
    zz0i = zz[ko,1,1]
    pr0i = pr[ko,1,1]
    sg0i = gsw.sigma0(sa0i,ct0i)
    rho0i = gsw.rho(sa0i,ct0i,pr0i)
    x_n = [x[io],]
    y_n = [y[jo],]
    ct_n = [ct0i,]
    sa_n = [sa0i,]
    zz_n = [zz0i,]
    pr_n = [pr0i,]

    # param for neutral surf calc
    Nabbe = 2                                   #to compute new depth of neutral surf
    indchunk = np.arange(ko-Nabbe,ko+Nabbe+1,1) #index centred on ko

    ############## 1 = cast (0,0)
    # specifics
    ct1c = ct[indchunk,0,0]
    sa1c = sa[indchunk,0,0]
    zz1c = zz[indchunk,0,0]
    p1c = pr[indchunk,0,0]
    x_n.append(x[io-1])
    y_n.append(y[jo-1])
    # automatic  
    p1, dct, dsa, dzz, flag = calc_gamman_pressureWGHC(ct0i,sa0i,zz0i,pr0i,ct1c,sa1c,zz1c,p1c,False,0) #neutral surf. no debugging
    pr_n.append(p1)
    ct_n.append(ct0i+dct)
    sa_n.append(sa0i+dsa)
    zz_n.append(zz0i+dzz)

    ############## 2 = cast (1, 0)
    # specifics
    ct1c = ct[indchunk,1,0]
    sa1c = sa[indchunk,1,0]
    zz1c = zz[indchunk,1,0]
    p1c = pr[indchunk,1,0]
    x_n.append(x[io-1])
    y_n.append(y[jo])
    # automatic
    p1, dct, dsa, dzz, flag = calc_gamman_pressureWGHC(ct0i,sa0i,zz0i,pr0i,ct1c,sa1c,zz1c,p1c,False,0) #neutral surf. no debugging
    pr_n.append(p1)
    ct_n.append(ct0i+dct)
    sa_n.append(sa0i+dsa)
    zz_n.append(zz0i+dzz)

    ############## 3 = cast (2, 0)
    # specifics
    ct1c = ct[indchunk,2,0]
    sa1c = sa[indchunk,2,0]
    zz1c = zz[indchunk,2,0]
    p1c = pr[indchunk,2,0]
    x_n.append(x[io-1])
    y_n.append(y[jo+1])
    # automatic
    p1, dct, dsa, dzz, flag = calc_gamman_pressureWGHC(ct0i,sa0i,zz0i,pr0i,ct1c,sa1c,zz1c,p1c,False,0) #neutral surf. no debugging
    pr_n.append(p1)
    ct_n.append(ct0i+dct)
    sa_n.append(sa0i+dsa)
    zz_n.append(zz0i+dzz)

    ############## 4 = cast (2, 1)
    # specifics
    ct1c = ct[indchunk,2,1]
    sa1c = sa[indchunk,2,1]
    zz1c = zz[indchunk,2,1]
    p1c = pr[indchunk,2,1]
    x_n.append(x[io])
    y_n.append(y[jo+1])
    # automatic
    p1, dct, dsa, dzz, flag = calc_gamman_pressureWGHC(ct0i,sa0i,zz0i,pr0i,ct1c,sa1c,zz1c,p1c,False,0) #neutral surf. no debugging
    pr_n.append(p1)
    ct_n.append(ct0i+dct)
    sa_n.append(sa0i+dsa)
    zz_n.append(zz0i+dzz)

    ############## 5 = cast (2, 2)
    # specifics
    ct1c = ct[indchunk,2,2]
    sa1c = sa[indchunk,2,2]
    zz1c = zz[indchunk,2,2]
    p1c = pr[indchunk,2,2]
    x_n.append(x[io+1])
    y_n.append(y[jo+1])
    # automatic
    p1, dct, dsa, dzz, flag = calc_gamman_pressureWGHC(ct0i,sa0i,zz0i,pr0i,ct1c,sa1c,zz1c,p1c,False,0) #neutral surf. no debugging
    pr_n.append(p1)
    ct_n.append(ct0i+dct)
    sa_n.append(sa0i+dsa)
    zz_n.append(zz0i+dzz)

    ############## 6 = cast (1, 2)
    # specifics
    ct1c = ct[indchunk,1,2]
    sa1c = sa[indchunk,1,2]
    zz1c = zz[indchunk,1,2]
    p1c = pr[indchunk,1,2]
    x_n.append(x[io+1])
    y_n.append(y[jo])
    # automatic
    p1, dct, dsa, dzz, flag = calc_gamman_pressureWGHC(ct0i,sa0i,zz0i,pr0i,ct1c,sa1c,zz1c,p1c,False,0) #neutral surf. no debugging
    pr_n.append(p1)
    ct_n.append(ct0i+dct)
    sa_n.append(sa0i+dsa)
    zz_n.append(zz0i+dzz)

    ############## 7 = cast (0, 2)
    # specifics
    ct1c = ct[indchunk,0,2]
    sa1c = sa[indchunk,0,2]
    zz1c = zz[indchunk,0,2]
    p1c = pr[indchunk,0,2]
    x_n.append(x[io+1])
    y_n.append(y[jo-1])
    # automatic
    p1, dct, dsa, dzz, flag = calc_gamman_pressureWGHC(ct0i,sa0i,zz0i,pr0i,ct1c,sa1c,zz1c,p1c,False,0) #neutral surf. no debugging
    pr_n.append(p1)
    ct_n.append(ct0i+dct)
    sa_n.append(sa0i+dsa)
    zz_n.append(zz0i+dzz)

    ############## 8 = cast (0, 1)
    # specifics
    ct1c = ct[indchunk,0,1]
    sa1c = sa[indchunk,0,1]
    zz1c = zz[indchunk,0,1]
    p1c = pr[indchunk,0,1]
    x_n.append(x[io])
    y_n.append(y[jo-1])
    # automatic
    p1, dct, dsa, dzz, flag = calc_gamman_pressureWGHC(ct0i,sa0i,zz0i,pr0i,ct1c,sa1c,zz1c,p1c,False,0) #neutral surf. no debugging
    pr_n.append(p1)
    ct_n.append(ct0i+dct)
    sa_n.append(sa0i+dsa)
    zz_n.append(zz0i+dzz)

    #### Fit a plane to the neutral surface
    xc = np.mean(x_n)
    yc = np.mean(y_n)
    zc = np.mean(zz_n)
    x_n = (x_n - xc)*dx_perdeg
    y_n = (y_n - yc)*dy_perdeg
    z_n = zz_n - zc
    MAT = np.array([x_n,y_n,z_n])
    try: #to handle error message when SVD does not converge or others
        V, S, U = np.linalg.svd(MAT)   #S(1)>S(2)>S(3)
    except:
        print("### calc_mvec subroutine: By Jove!", sys.exc_info()[0], "occurred.")
        V = np.nan*np.ones((3,9),dtype=float)
    v1 = V[:,0]
    v2 = V[:,1]
    v3 = V[:,2]
    MAT = MAT.T #for display
    #NB: if the SVD is performed as below on the transposed MAT, can't get the right answer although it should be the same...?!
    #U, S, V = np.linalg.svd(MAT)   #S(1)>S(2)>S(3)
    #v1 = V[:,0]
    #v2 = V[:,1]
    #v3 = V[:,2]
    mvec_fit = v3 / v3[2]       #m-vector with unit z-component
    umvec_fit = v3  #3rd right singular vector (unit vector)
    umvec_dir = mvec_dir / np.sqrt(np.dot(mvec_dir,mvec_dir))      #unit vector for the direct method
    #print('#### calc_mvec subroutine ###',mvec_dir, mvec_fit)

    ############################################################ Figures
    if FIG==1:
        xp = np.array([-50,50,50,-50, -50])*km  #row vector of x-positions
        yp = np.array([-50, -50,50,50, -50])*km
        zp1 = -(v1[0]*xp+v1[1]*yp)/v1[2]  #corresponding z positions, to plot the plane
        zp2 = -(v2[0]*xp+v2[1]*yp)/v2[2]  #corresponding z positions, to plot the plane
        zp3 = -(v3[0]*xp+v3[1]*yp)/v3[2]  #corresponding z positions, to plot the plane
        zpd = -(umvec_dir[0]*xp+umvec_dir[1]*yp)/umvec_dir[2]     #corresponding z positions, to plot the plane
        refpoint = np.array([0, 0, 0])
        plt.figure()
        #axs = plt.axes(projection="3d")
        axs = plt.subplot(111,projection='3d')
        axs.scatter(MAT[:,0]/km, MAT[:,1]/km, MAT[:,2], color='k')
        axs.plot_wireframe(xp/km, yp/km, zp3, color='m')
        axs.plot_wireframe(xp/km, yp/km, zpd, color='g')
        #axs.quiver([0,],[0,],[0,],[umvec_fit[0]*km,],[umvec_fit[1]*km,],[umvec_fit[2],],color='m')
        #axs.quiver([0,],[0,],[0,],[umvec_dir[0]*km,],[umvec_dir[1]*km,],[umvec_dir[2],],color='g')
        axs.set_xlabel('east west direction (km)')
        axs.set_ylabel('north south direction (km)')
        axs.set_zlabel('vertical direction (m)')
        plt.show()
    return mvec_fit, mvec_dir

def calc_stab(sa,ct,pr,zz,ko):
    """ Computes N2 and contributions from salinity and temperature, as well as the thermal and
    haline expansion coeff, and in situ density
    INPUT --> same as calc_mvec
    OUTPUT --> N2, N2ct, N2sa, all in s-2
               alpha [1/T], beta [1/S], rh0i [kg/m3, e.g. 1025]
    """
    ct0i = ct[ko,1,1]
    sa0i = sa[ko,1,1]
    zz0i = zz[ko,1,1]
    pr0i = pr[ko,1,1]
    sg0i = gsw.sigma0(sa0i,ct0i)
    rho0i = gsw.rho(sa0i,ct0i,pr0i)
    alpha = gsw.alpha(sa[ko,1,1],ct[ko,1,1],pr[ko,1,1])
    beta = gsw.beta(sa[ko,1,1],ct[ko,1,1],pr[ko,1,1])
    gradCTz = (ct[ko+1,1,1]-ct[ko-1,1,1])/(zz[ko+1,1,1]-zz[ko-1,1,1])
    gradSAz = (sa[ko+1,1,1]-sa[ko-1,1,1])/(zz[ko+1,1,1]-zz[ko-1,1,1])
    N2ct = gravi*alpha*gradCTz
    N2sa = -gravi*beta*gradSAz
    N2 = N2ct+N2sa
    return N2, N2ct, N2sa, alpha, beta, rho0i


def calc_Bvector(sa,ct,pr,zz,x,y,io,jo,ko):
    """ Computes the baroclinicity vector B, a 1D array (Bx, By, Bz) in s-2, at the central location of the "cube"
    Assumes that the dominant term in the baroclinicity is the slope of teh neutral surfaces
    INPUT  --> same as calc_mvec
    OUTPUT --> Baroc = (Bx, By, Bz) in s-2
               umvec = unit vector normal to B (parallel to mvec)
               unvec = unit vector normal to B (parallel to gradn_P)
    """
    ### Neutral plane calc.
    mvec_fit, mvec_dir = calc_mvector(sa,ct,pr,zz,x,y,io,jo,ko)
 
    ### Stability
    N2, N2ct, N2sa, alpha, beta, rho = calc_stab(sa,ct,pr,zz,ko)

    ### Computes gradn_P
    gradn_Px = mvec_fit[0]*rho*gravi #in Pa/m
    gradn_Py = mvec_fit[1]*rho*gravi #in Pa/m
    gradn_P = np.array([gradn_Px, gradn_Py, 0])
   
    #### Computes Baroclinicity vector
    Baroc = np.cross(gradn_P,mvec_fit)*N2/(rho*gravi)

    ###Unit vectors orthogonal to Baroclinicity vector
    umvec = mvec_fit / np.sqrt(np.dot(mvec_fit,mvec_fit) )
    unvec = gradn_P / np.sqrt(np.dot(gradn_P,gradn_P) )
    return Baroc, umvec, unvec


def calc_Bvector_approx(sa,ct,pr,zz,x,y,io,jo,ko,LINeos):
    """Computes the baroclinicity vector at (io,jo,ko) directly from the horizontal gradients of CT and SA.
    The calculation neglects entirely the horizontal pressure gradient contribution and so returns a 3D vector
    with zero vertical component (Bx, By, Bz=0).
    LINeos is a vector:
        *If LINeos[0]=True then alpha and beta are not functions of the actual SA and CT but
        use a ref value specified by LINeos[1]=CTref and LINeos[2]=SAref.
        *If LINeos[0]=False then the full calculation is performed and LINeos[1], LINeos[2] are not used
   
    ---> INPUTS: 9 "casts" in total forming a "cube":
         sa,ct,pr,zz all of dimensions (Nz,3,3)
         x and y are the lon and latitude 1D-arrays
         io,jo,ko are the indices of the central point
         The ordering is as follows:
              list (j,i): [central point, (0,0), (1,0), (2,0) ..., (0,1)]
              i.e. clockwise from sw corner of the reference point.
              Label points with ref as 0 and then clockwise, so that
              the labelling fits that of the lists x_n, y_n, etc.
 
                     3    4     5
                     2  ref=0   6
                     1    8     7

    --> OUTPUTS: the Baroclinic vector (Bx, By, Bz=0), buoyancy gradient vector "grad_buo",
                 and a local buoyancy variable "buo"

    NB: does not use all the data in the "cube" but easier to use this input format

    """
    #import sys
    dx = RADIUS*np.cos(y[jo]*np.pi/180)*0.5*np.pi/180 #distance between two gridpoints  
    dy = RADIUS*0.5*np.pi/180
    CTref = 3.64  #glob=3.64
    SAref = 34.72   #glob=34.72

    #### Computes horizontal density gradients
    if LINeos[0] == False:
        alpha = gsw.alpha(sa[ko,1,1],ct[ko,1,1],pr[ko,1,1])
        beta = gsw.beta(sa[ko,1,1],ct[ko,1,1],pr[ko,1,1])
    else:
        CTref = LINeos[1] #overrides CTref
        SAref = LINeos[2] #overrides SAref
        alpha = gsw.alpha(SAref,CTref,pr[ko,1,1])
        beta = gsw.beta(SAref,CTref,pr[ko,1,1])
    gradCTx = (ct[ko,1,2]-ct[ko,1,0])/(2*dx)
    gradCTy = (ct[ko,2,1]-ct[ko,0,1])/(2*dy)
    gradCT = np.array([gradCTx,gradCTy,0*gradCTx]) #row vector
    gradSAx = (sa[ko,1,2]-sa[ko,1,0])/(2*dx)
    gradSAy = (sa[ko,2,1]-sa[ko,0,1])/(2*dy)
    gradSA = np.array([gradSAx,gradSAy,0*gradSAx])
    grad_buo = (alpha*gradCT-beta*gradSA)  #z-component =0

    ### Local buoyancy variable
    rhoo = gsw.rho(SAref,CTref,pr[ko,1,1])
    buo = rhoo*(alpha*(ct[ko,1,1]-CTref) - beta*(sa[ko,1,1]-SAref))/gravi

    ### Gravity
    gravi_vec = np.array([0,0,-gravi])

    #### Computes Baroclinicity vector
    Baroc = np.cross(gravi_vec,grad_buo)
   
    return Baroc, grad_buo, buo


def calc_dV(botd0,k0,zvec,lat0):
    """
    computes the volume occupied by gridbox at lon0,lat0, depth=zvec[ko] (<0) and with bottom depth BOTD[jo,io]=botd0
    """
    DEBUG = 0
    Nz = len(zvec)
    dx = RADIUS*np.cos(lat0*np.pi/180)*0.5*np.pi/180  #in m
    dy = RADIUS*0.5*np.pi/180  #in m
    dS = dx*dy
    zbot0 = -botd0 #to convert into a depth
    ktop = 0  #index of level just below the sea surface (1st level in dataset is z=0)
    crit = zvec-zbot0
    if ma.getmask(zbot0)==True or zbot0>=0: #land grid points
        dV = np.nan
    else:
        indok = np.where(crit>0)  #levels above seafloor
        indok = indok[0]
        kbot =  indok[-1] #index of level just above seafloor
        if DEBUG ==1:
            print('#####SUBROUTINE calc_Vol!!', k0,zvec[k0],kbot,zvec[kbot],zvec[kbot+1],zbot0)
        if k0 == ktop:
            dz = -(zvec[0]+zvec[1])/2
        elif k0 == kbot:
            dz = (zvec[kbot]+zvec[kbot-1])/2 - zbot0
        elif k0 > kbot:  #below seafloor
            dz = np.nan
        else:
            dz = (zvec[k0-1]-zvec[k0+1])/2
        dV = dx*dy*dz   #in m3
        if DEBUG==1:
            print('#####SUBROUTINE calc_Vol!!',dz)
    return dV, dS


def calc_PGzeta(Baroc,f):
    """
    computes the vorticity of a fluid element in the PG limit (-dv/dz, du/dz, f)

    INPUT--> Baroc: Baroclinicity vector, 1D-array (Bx, By, Bz) in s-2  
             f: the Coriolis parameter at the given location

    OUTPUT--> zeta = 1D array (zetax,zetay,zetaz) in s-1
    """
    zeta = np.array([Baroc[1]/f, -Baroc[0]/f, f])
    return zeta


def calc_HelMat(Baroc,f):
    """
    computes the "Helmholtz matrix" M

    INPUT--> Baroc: Baroclinicity vector, 1D-array (Bx, By, Bz) in s-2  
             f: the Coriolis parameter at the given location

    OUTPUT--> M = the "Helmholtz matrix" (3X3, each element in units of s-1)
                  currently only includes the dudz and dv/dz components
    """
    M = np.zeros((3,3))
    M[0,2] = -Baroc[0]/f
    M[1,2] = -Baroc[1]/f
    return M

def calc_stretchtilt(M,ds):
    """
    computes the rate of change of the shape of the filament ds between times t and t+dt from the "Helmholtz" matrix M

    INPUT--> M = 3X3 matrix
             ds = 1D array (dsx,dsy,dsz)

    OUTPUT--> the matrix product M ds
    """
    dsdt = np.matmul(M,ds)
    return dsdt