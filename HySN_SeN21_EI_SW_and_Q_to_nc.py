import sys, os, glob
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Helene B. Erlandsen
from datetime import date
today = str(date.today())
import numpy as np
import pandas as pd
import xarray as xr
import xarray.ufuncs as xu
from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster
client = Client(processes=False)
#lc=LocalCluster()
#w=lc.start_worker(ncores=3)

import gc #garbage collector
from pyproj import Proj
#The interpolation could be done with scipy interp,
#Basemap interp or xarray interp
#choosing to use xesmf cause saves weights, needs pip install 
import xesmf as xe #pip install / conda pip install
from pyproj import Proj

#import cProfile # if you want to print run time of script
#pr = cProfile.Profile()

import matplotlib.pyplot as plt


datadir="data/"
EIpath=datadir+'ERAI/'
#savedir='/lustre/storeB/users/helenebe/HySN/'
savedir=datadir+"HySN/"

startyr=1981
endyr=2016
#Make precip_24hours_means ie 06 UTC to 06 UTC set RRday to True
RRday=True

if RRday==True:
    dayshh=6
    lab='right'
else:
    dayshh=0
    lab='left'

#Constants 

T0NVE=273.1
#-----------------------------------------------------
#          EC constants from IFS manual and mars wiki
#          Td is calculated for saturation above water
ecg=9.80665 #gravitational constant
T0=273.16 #K
Rd=287.0597 #JK 
Rv=461.5250 #JK
eps=Rd/Rv
R=287.# J/kg/K

#Buck(1981) used to calculate Td in Era-interim
aw=611.21
bw=17.502
cw=240.97
#----------------------------------------------------
#AERKi to calculate Tf (freeze point temperature) 
ai=aw
bi=22.587
ci=273.86

# SW regridding
tau=0.720
alpha=-1.5E-5

#SeNorge initialization -----------------------------
SeNorgegeop='http://thredds.met.no/thredds/dodsC/senorge/geoinfo/seNorge2_dem_UTM33.nc'
dx=1000.

SeNorge21T2p='http://thredds.met.no/thredds/dodsC/senorge/seNorge2_1/TEMP1d/seNorge_v2_1_TEMP1d_grid_'
SeNyear=1979 # random

SeN=xr.open_dataset(SeNorge21T2p+str(SeNyear)+'.nc',chunks={'time': 30})#cache=False) #
projs=SeN['UTM_Zone_33'].proj4
myP=Proj(projs)

#maybe make def below ------------------------------------------
Xcorners=np.arange(SeN['X'].data[0]-dx/2., SeN['X'].data[-1]+3*dx/2., dx)
Ycorners=np.flipud(np.arange(SeN['Y'].data[-1]-dx/2., SeN['Y'].data[0]+3*dx/2., dx))
Lon2, Lat2 = myP(*np.meshgrid(SeN['X'].data,SeN['Y'].data),inverse=True)
Lon2b, Lat2b = myP(*np.meshgrid(Xcorners,Ycorners),inverse=True) #

lons=np.asarray(Lon2)
lats=np.asarray(Lat2)
SeN.coords['lat'] = (('Y','X'),Lat2)
SeN.coords['lon'] = (('Y','X'),Lon2)
#SeN.set_coords(['lat','lon'])

SeN.coords['Xb'] = (Xcorners)
SeN.coords['Yb'] = (Ycorners)
SeN.set_coords(['Xb','Yb'])

SeN.coords['lat_b'] = (('Yb','Xb'),Lat2b)
SeN.coords['lon_b'] = (('Yb','Xb'),Lon2b)
SeN.set_coords(['lat_b','lon_b'])

write_to_file=False
if write_to_file:
    mapSN=SeN.drop('mean_temperature')#][0,:,:]
    mapSN=mapSN.drop('time')
    mapSN.set_coords(['lat_b','lon_b'])
    mapSN.to_netcdf('SeN_UTM33_with_crns_4_HySN.nc', mode='w', format='NETCDF4',
                    group=None, engine='netcdf4',encoding={'lon': {'_FillValue':False},
                                                           'lat': {'_FillValue':False},
                                                           'lon_b': {'_FillValue':False},
                                                           'lat_b': {'_FillValue':False}})
#-------------

SeNoro=xr.open_dataset('http://thredds.met.no/thredds/dodsC/senorge/geoinfo/seNorge2_dem_UTM33.nc')
SeNorge_oro=SeNoro['elevation']
SeNorgemask=SeNorge_oro.isnull() #True/1 for water, False/0 for land
SeNorgelsm=SeNorge_oro.notnull()*1.0 #0. for water, 1. for land
SeN.coords['mask'] = (('Y', 'X'), SeNorgelsm)
SeN['orography'] = (('Y', 'X'), SeNorge_oro)

#SeNorge is indexed as precip days with the stamp reflecting
#the accumulated 24 hour precip until the date at 06 hours  

#Era init-------------------------------------------------------------------------
#Era-Interim data
#Era-Interim 6-hourly SW and LW files, from forecasts at 00 and 12 UTC
Era_rad_f00 = xr.open_dataset(EIpath+'LW_SW_fc_EI_f00_1979_2017.nc',cache=False)
Era_rad_f12 = xr.open_dataset(EIpath+'LW_SW_fc_EI_f12_1979_2017.nc',cache=False)
#Era-Interim 2-meter temperature (T2) and 2-meter dew point temperature, from analysis fields
Era_termo_a =  xr.open_dataset(EIpath+'t2m_d2m_ps_anaEI_1979_2017.nc',cache=False)
#Era-Interim orography
Era_oro = xr.open_dataset(EIpath+'oro_veg_lmaskEClim.nc')

# Era-Interim grid

eclat=Era_rad_f00['latitude'][:]
eclon=Era_rad_f00['longitude'][:]
dl=abs(eclat.diff('latitude').data[0])
lonc=np.arange(eclon.data[0]-dl/2., eclon.data[-1]+3*dl/2., dl)
latc=np.flipud(np.arange(eclat.data[-1]-dl/2., eclat.data[0]+3*dl/2., dl))

def fixcoords(var):
    var=var.rename({'longitude': 'lon', 'latitude': 'lat'})
    var.coords['lat_b'] = (latc)
    var.coords['lon_b'] = (lonc)
    var.set_coords(['lat_b','lon_b'])
    return var

Era_rad_f00=fixcoords(Era_rad_f00)
Era_rad_f12=fixcoords(Era_rad_f12)
Era_oro=fixcoords(Era_oro)
Era_termo_a=fixcoords(Era_termo_a)

ecsgeo=Era_oro['z'][0]  
ecmaskb=Era_oro['lsm'][0] #True/1 for water
ecnamask=ecmaskb.where(ecmaskb==1.,np.nan)
ecmask=ecnamask.notnull() 
ECoro=ecsgeo/ecg
ECoro.attrs['units'] = 'meters'
ECoro.name='orography'

Era_oro.coords['mask'] = (('lat', 'lon'), ecmaskb)
Era_oro['orography']=(['lat', 'lon'], ECoro.data, ECoro.attrs)
savEra=False
if savEra:
    Era_oro.to_netcdf('./data/ERAI/Era_with_crns.nc',encoding={'lon': {'_FillValue':False},'lat': {'_FillValue':False},'lon_b': {'_FillValue':False},'lat_b': {'_FillValue':False}})

#Making daily means with prcipitation days 06-06UTC:
T2E=Era_termo_a['t2m'].resample(time='24h', base=dayshh,label=lab).mean(dim='time')
T2dE=Era_termo_a['d2m'].resample(time='24h', base=dayshh,label=lab).mean(dim='time')
psE=Era_termo_a['sp'].resample(time='24h', base=dayshh,label=lab).mean(dim='time')

T2E[0, :,:]=np.nan # first day is not 24 hours for radiation
T2E=T2E.dropna('time',how='all')
T2dE[0, :,:]=np.nan # first day is not 24 hours for radiation
T2dE=T2dE.dropna('time',how='all')
psE[0, :,:]=np.nan # first day is not 24 hours for radiation
psE=psE.dropna('time',how='all')

#Era-Interim Vapour pressure
e2w=aw*np.exp(bw*(T2dE-T0)/(T2dE-T0 +cw))#/1000. #for water, ref ifs copernicus knowledge base

#Era-Interim saturation vapour pressure
e2si=ai*np.exp(bi*(T2E-T0)/(T2E-T0 +ci))
e2sw=aw*np.exp(bw*(T2E-T0)/(T2E-T0 +cw))

#do not need to use enhancement factors because they are cancelled for RH 
RHEw=e2w/e2sw*100.
RHEi=e2w/e2si*100.
#check RH not more than 100%
RHEi=RHEi.where(RHEi<100.,100.)
RHEw=RHEw.where(RHEw<100.,100.)



#Incident SW at two forecast times ---------------------------
# 00 has 12-24, 12 has 24(00)-12 to have spin up in diag fields
SW00acc=Era_rad_f00['ssrd']
SW12acc=Era_rad_f12['ssrd']
LW00acc=Era_rad_f00['strd']
LW12acc=Era_rad_f12['strd'] 

#The fluxes are in two files from two forecast times. Want half a day from each file
#previous 6 hours average inst SW flux starts at ind 1 not 0. units are [J/m2]:
#to get [J/s/m2=W/m2] divide with 3600*6

#will drop first date and the second date now holds the past 6hrs mean flux
SW00=SW00acc.diff('time')/(6.*60.*60.) 
SW12=SW12acc.diff('time')/(6.*60.*60.)
LW00=LW00acc.diff('time')/(6.*60.*60.) 
LW12=LW12acc.diff('time')/(6.*60.*60.)

#Every third acc. flux holds the diff between last time
# in last forecast and newest in next: delete 2::3
SW00[2::3, :,:]=np.nan 
SW00=SW00.dropna('time',how='all')
SW12[2::3, :,:]=np.nan
SW12=SW12.dropna('time',how='all')
SW=SW00.combine_first(SW12) #combine the two

#Every third acc. flux holds the diff between last time
# in last forecast and newest in next: delete 2::3
LW00[2::3, :,:]=np.nan 
LW00=LW00.dropna('time',how='all')
LW12[2::3, :,:]=np.nan
LW12=LW12.dropna('time',how='all')
LW=LW00.combine_first(LW12) #combine the two

SW=SW.resample(time='24h', base=dayshh,label=lab).mean(dim='time') #precip day mean
SW[0, :,:]=np.nan # first day is not 24 hours
SW=SW.dropna('time',how='all')

LW=LW.resample(time='24h', base=dayshh,label=lab).mean(dim='time') #precip day mean
LW[0, :,:]=np.nan # first day is not 24 hours
LW=LW.dropna('time',how='all')

# alligns the SW data with the rest
# specific to my two EC-files that this is enough
SW,dump=xr.align(SW,RHEw, join='inner')
LW,dump=xr.align(LW,RHEw, join='inner')


#Regridding this is based on
#pip install --upgrade git+https://github.com/JiaweiZhuang/xESMF.git@masking
#(https://github.com/JiaweiZhuang/xESMF/issues/22)
#regridder_consnorm=xe.Regridder(Era_oro, SeN, 'conservative_normed',reuse_weights=True)

regridder=xe.Regridder(Era_oro, SeN, 'bilinear',reuse_weights=True)

regridder_s2d = xe.Regridder(Era_oro, SeN, 'nearest_s2d',reuse_weights=True)

dr_mask_bil=regridder(Era_oro['mask'])
extrapolate=dr_mask_bil.where(dr_mask_bil!=0,-99999)
dr_mask=regridder_s2d(Era_oro['mask'])

#after much thought using bilinear on all except near coast for consistency
#1:since Era-I is already bil-interp from MARS
#2:since going from 0.75deg to 1km --> diff only a few points
#3:on all vars and params for consistency:
#i.e. if Eras T2 is bilinearly interp and then oro is not,
# the ds will be somewhat inconsistent as T2=b+lapse*z
dr_oro_bil=regridder(Era_oro['orography'])
dr_oro_s2d=regridder_s2d(Era_oro['orography'])
dr_oro=dr_oro_bil.where(extrapolate!=-99999, dr_oro_s2d)
gc.collect()

#--DZ--------
dz=SeN['orography']-dr_oro #dz.plot();plt.show()
#dzmask=ma.masked_array(dz,Commask) #masked to common land points
#print('mean elevation diff (SeNorge-Erai):', ma.mean(dzmask))


def regrid2step(data):
    dr_out_bil=regridder(data)
    #dr_out_cons_norm=regridder_consnorm(data)
    dr_outs2d=regridder_s2d(data)
    #outside mapped places nearest s2d:  
    dr_out=dr_out_bil.where(extrapolate != -99999, dr_outs2d)
    #add coords
    dr_out.coords['X'] = (SeN['X'].data)
    dr_out.coords['Y'] = (SeN['Y'].data)
    dr_out.coords['mask']=SeN['mask'].astype('i2')
    dr_out.name=data.name
    dr_out=dr_out.astype('f4')
    return dr_out

def writetotmp_nc(data,pn):
    filname='tmp/'+data.name+pn+'.nc'
    data.to_netcdf(filname,encoding={data.name:{'dtype': 'f4', '_FillValue': -9999}})
    #faster IO w/o zlib, but takes much more space 
    
# To work with pc memory, and due to issues with dask+opendap
# the data is sliced in time
# Three bolks a year ['p1', 'p2','p3']
#p1ul='<120'; p2ll='>=120'; p2ul='<140'; p3ll='>=140'
loweri = {}; upperi={}
loweri['p1']=0
loweri['p2']=120
loweri['p3']=240
upperi['p1']=120
upperi['p2']=240
upperi['p3']=None

for filename in glob.glob('tmp/*') :
    os.remove( filename )
    
for yr in np.arange(startyr,endyr): #make loop
    SeNT=xr.open_dataset(SeNorge21T2p+str(yr)+'.nc',cache=False) #temperature in Celsius
    if RRday: 
        SeNT.coords['time']=SeNT.indexes['time']+pd.Timedelta(hours=6) #adding the six hours

    SeNT.encoding={'time': dict(unlimited=False)}

    for pn in ['p1', 'p2','p3']:
        SeNTu=SeNT['mean_temperature'][loweri[pn]:upperi[pn],:,:]
        SeNTu,SWp=xr.align(SeNTu,SW, join='inner')
        SeNTu,LWp=xr.align(SeNTu,LW, join='inner')
        SeNTu,psp=xr.align(SeNTu,psE, join='inner')
        SeNTu,RHwp=xr.align(SeNTu,RHEw, join='inner')
        RHwp.name='RHw'
        SeNTu,RHip=xr.align(SeNTu,RHEi, join='inner')
        RHip.name='RHi'
        SeNTu,T2p=xr.align(SeNTu,T2E, join='inner')
        SeNTu,e2p=xr.align(SeNTu,e2w, join='inner')
        e2p.name='e2'

        writetotmp_nc(SeNTu,pn=pn)

        #--Regridding some variables----------
        print('regridding '+pn)
        Mps=regrid2step(psp)
        writetotmp_nc(Mps,pn=pn)
        Epsn=Mps.name
        del Mps
        MT2=regrid2step(T2p)
        writetotmp_nc(MT2,pn=pn)
        ET2n=MT2.name
        del MT2
        MRHw=regrid2step(RHwp)
        writetotmp_nc(MRHw,pn=pn)
        ERHwn=MRHw.name
        del MRHw
        gc.collect()
        MRHi=regrid2step(RHip)
        writetotmp_nc(MRHi,pn=pn)
        ERHin=MRHi.name
        del MRHi
        MSW=regrid2step(SWp)
        writetotmp_nc(MSW,pn=pn)
        ESWn=MSW.name
        del MSW
        MLW=regrid2step(LWp)
        writetotmp_nc(MLW,pn=pn)
        ELWn=MLW.name
        del MLW
        Me2=regrid2step(e2p)
        writetotmp_nc(Me2,pn=pn)
        Ee2n=Me2.name
        del Me2
        print('regridded '+pn+' '+str(yr))
        gc.collect()
    del SeNTu

    gc.collect()

    # START ------------------------------------------------------------
    print('Starting downscaling')

    #The first part is based on Cosgrove and similar to the PGF/WATCH downscaling

    SeNT=xr.open_mfdataset('tmp/mean_temperature'+'*.nc',chunks={'time': 10})
    SeNTr=SeNT['mean_temperature']
    SeNTr.attrs={}

    EraT=xr.open_mfdataset('tmp/'+T2E.name+'*.nc',chunks={'time': 10})
    EraTr=EraT[T2E.name]

    #calculating hypso adjusted pressure
    Eps=xr.open_mfdataset('tmp/'+psE.name+'*.nc',chunks={'time': 10})
    Epsr=Eps[psE.name]

    pSN=Epsr/xu.exp(ecg*dz/(Rd*((SeNTr+T0NVE)+EraTr)/2.)) #Kelvin

    #Make DataArray to nice dataset
    pSN.name='sp'
    pSN.attrs['longname']='surface pressure'
    pSN.attrs['unit']='Pa'
    pSN.attrs['notes']='Instantaneous, average of the last 24 hours, sampled every six hours.'
    pSN.attrs['date_created'] = today
    pSN.attrs['license'] = 'Norwegian Licence for Open Government Data (NLOD), https://data.norge.no/nlod/en/1.0'
    pSN.attrs['creator_name'] = 'Helene B. Erlandsen'
    pSN.attrs['proj4']=projs

    delayed_obj = pSN.to_netcdf(savedir+'HySN_Surface_Pressure_'+str(yr)+'.nc',
                                compute=False,unlimited_dims=None,
                                encoding={pSN.name:{'dtype': 'i4','_FillValue': -9999,'zlib':True},
                                          'lon': {'dtype': 'f4','_FillValue':False},
                                          'lat': {'dtype': 'f4','_FillValue':False},
                                          'X': {'dtype': 'f4','_FillValue':False},
                                          'Y': {'dtype': 'f4','_FillValue':False},
                                          'mask':{'dtype':'i2','_FillValue': -9999},
                                          'time':{'dtype':'f8'}}) 

    print('Writing downscaled surface pressure for year '+str(yr) +' to file:'+
          savedir+'HySN_Surface_Pressure_'+str(yr)+'.nc')
    with ProgressBar():
          results = delayed_obj.compute()

    gc.collect()

    ##-------------------------------------------------------------------
    #Calculating adjusted vapour pressure assuming constant RH with height
    ERHw=xr.open_mfdataset('tmp/RHw'+'*.nc',chunks={'time': 10})
    ERHwr=ERHw['RHw']
    ERHi=xr.open_mfdataset('tmp/RHi'+'*.nc',chunks={'time': 10})
    ERHir=ERHi['RHi']

    #in Celsius:
    TdH0w=cw*(xu.log(ERHwr/100.) + bw*SeNTr/(cw+SeNTr))/(bw-xu.log(ERHwr/100.)-bw*SeNTr/(cw+SeNTr))
    TdH0i=ci*(xu.log(ERHir/100.) + bi*SeNTr/(ci+SeNTr))/(bi-xu.log(ERHir/100.)-bi*SeNTr/(ci+SeNTr)) 

    fw=(1.00071*xu.exp(0.0000045*pSN/100.)) #Ald Esk 1996 eq 17
    fi=(0.99882*xu.exp(0.000008*pSN/100.)) #Ald Esk 1996 eq 18

    eHw=fw*aw*xu.exp(bw*TdH0w/(cw+TdH0w))# Buck
    eHi=fi*ai*xu.exp(bi*TdH0i/(ci+TdH0i))# AERK ice eq 24 want temp in C, p in hPa     
    eH=eHi.where(SeNTr<0, eHw) #test if tmean is lower than zero

    qH=eps*eH/(pSN-(1-eps)*eH) #are all in Pa
    qH.name='huss'
    qH.attrs['longname']='Near-Surface Specific Humidity'
    qH.attrs[ 'unit']='kg kg-1'
    qH.attrs['conversion_formula']='eps*e/(ps(1-eps)*e)'
    qH.attrs['notes']='Instantaneous, average of the last 24 hours, sampled every six hours.'
    qH.attrs['date_created'] = today
    qH.attrs['license'] = 'Norwegian Licence for Open Government Data (NLOD), https://data.norge.no/nlod/en/1.0'
    qH.attrs['creator_name'] = 'Helene B. Erlandsen'
    qH.attrs['proj4']=projs

    delayed_obj = qH.to_netcdf(savedir+'HySN_Near_Surface_Specific_Humidity_'+str(yr)+'.nc',
                               compute=False,unlimited_dims=None,
                               encoding={qH.name:{'dtype': 'f8','_FillValue': -9999.,
                                                  'least_significant_digit':6,'zlib':True},
                                         'lon': {'dtype': 'f4','_FillValue':False},
                                         'lat': {'dtype': 'f4','_FillValue':False},
                                         'X': {'dtype': 'f4','_FillValue':False},
                                         'Y': {'dtype': 'f4','_FillValue':False},
                                         'mask':{'dtype':'i2','_FillValue': -9999},
                                         'time':{'dtype':'f8'}})

    print('Writing downscaled specific humidity for year '+str(yr) +' to file:'
          +savedir+'HySN_Near_Surface_Specific_Humidity_'+str(yr)+'.nc')
    with ProgressBar():
          results = delayed_obj.compute()

    #Longwave clear sky scaling 
    Ee=xr.open_mfdataset('tmp/e2'+'*.nc',chunks={'time': 10})
    Eer=Ee['e2']
    ET2=xr.open_mfdataset('tmp/t2m'+'*.nc',chunks={'time': 10})
    ET2r=ET2['t2m']
    ELW=xr.open_mfdataset('tmp/strd'+'*.nc',chunks={'time': 10})
    ELWr=ELW['strd']

    #check power in dask xr
    epsE=1.08*(1-xu.exp(-((Eer/100.)**(ET2r/2016.))))#hpa correct
    epsH=1.08*(1-xu.exp(-((eH/100.)**((SeNTr+T0NVE)/2016.)))) #satterlund 1969 
    sca=(epsH/epsE)*((SeNTr+T0NVE)/ET2r)**4
    LWSN=sca*ELW

    LWSN=LWSN.rename({'strd':'rlds'})
    #LWSN.name='rlds' 
    LWSN['rlds'].attrs['longname']='Surface Downwelling Longwave Radiation'
    LWSN['rlds'].attrs['unit']='W m-2'
    LWSN['rlds'].attrs['notes']='Average over the last 24 hours. Positive downwards'
    LWSN.attrs['date_created'] = today
    LWSN.attrs['license'] = 'Norwegian Licence for Open Government Data (NLOD), https://data.norge.no/nlod/en/1.0'
    LWSN.attrs['creator_name'] = 'Helene B. Erlandsen'
    LWSN.attrs['proj4']=projs

    delayed_obj = LWSN.to_netcdf(savedir+'HySN_Surface_Downwelling_Longwave_Radiation_'+str(yr)+'.nc',
                                 compute=False,unlimited_dims=None,
                                 encoding={'rlds':{'dtype': 'f8','_FillValue': -9999.,
                                                   'least_significant_digit':2,'zlib':True},
                                           'lon': {'dtype': 'f4','_FillValue':False},
                                           'lat': {'dtype': 'f4','_FillValue':False},
                                           'X': {'dtype': 'f4','_FillValue':False},
                                           'Y': {'dtype': 'f4','_FillValue':False},
                                           'mask':{'dtype':'i2','_FillValue': -9999},
                                           'time':{'dtype':'f8'}})

    print('Writing downscaled longwave radiation for year '+str(yr) +
          ' to file:'+savedir+'HySN_Surface_Downwelling_Longwave_Radiation_'+str(yr)+'.nc')
    with ProgressBar():
          results = delayed_obj.compute()

    ###---------------------------------------------------------------
    ### Cosgrove end ------------------------------------------------

    # SW clear sky scaling T&R --------------------------------------

    #tau_H=np.power(tau,pSN/101300.)+alpha*eH
    tau_H=tau**(pSN/101300.)+alpha*eH
    tau_E=tau**(Epsr/101300.)+alpha*Eer
    Taut_ratio=tau_H/tau_E

    ESW=xr.open_mfdataset('tmp/ssrd'+'*.nc',chunks={'time': 10})
    ESWr=ESW['ssrd']

    SWSN=(Taut_ratio**2)*ESWr
    SWSN.name='rsds'
    SWSN.attrs['longname']='Surface Downwelling Shortwave Radiation'
    SWSN.attrs['unit']='W m-2'
    SWSN.attrs['notes']='Average over the last 24 hours. Positive downwards'
    SWSN.attrs['date_created'] = today
    SWSN.attrs['license'] = 'Norwegian Licence for Open Government Data (NLOD),https://data.norge.no/nlod/en/1.0'
    SWSN.attrs['creator_name'] = 'Helene B. Erlandsen'
    SWSN.attrs['proj4']=projs
    #SWSN.attrs['creator_email'] = 'hebe@nve.no'
    #ncfile.Proj4_string=


    delayed_obj=SWSN.to_netcdf(savedir+'HySN_Surface_Downwelling_Shortwave_Radiation_'+str(yr)+'.nc',
                               compute=False,unlimited_dims=None,
                               encoding={'rsds':{'dtype': 'f8','_FillValue': -9999.,
                                                 'least_significant_digit':1,'zlib':True},
                                         'lon': {'dtype': 'f4','_FillValue':False},
                                         'lat': {'dtype': 'f4','_FillValue':False},
                                         'X': {'dtype': 'f4','_FillValue':False},
                                         'Y': {'dtype': 'f4','_FillValue':False},
                                         'mask':{'dtype':'i2','_FillValue': -9999},
                                         'time':{'dtype':'f8'}})

    print('Writing downscaled shortwave radiation for year '+str(yr) +
          ' to file:'+savedir+'HySN_Surface_Downwelling_Shortwave_Radiation_'+str(yr)+'.nc')
    with ProgressBar():
          results = delayed_obj.compute()

    for filename in glob.glob('tmp/*') :
        os.remove( filename )      

print('done all')


 
