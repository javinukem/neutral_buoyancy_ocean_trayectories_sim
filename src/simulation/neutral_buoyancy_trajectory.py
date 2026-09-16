
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
