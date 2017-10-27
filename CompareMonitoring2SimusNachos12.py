#!/usr/bin/env python
#
"""
This script is making some plots from annual means of NACHSO12 outputs using cartopy
"""

## imports

import sys
import numpy as np
import numpy.ma as ma
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def all_plots(case1,case2,year,**kwargs):

  dirmean1='/scratch/cnt0024/hmg2840/albert7a/NACHOS12.L75/NACHOS12.L75-'+case1+'-MEAN/1d/'+year+'/'
  dirmean2='/scratch/cnt0024/hmg2840/albert7a/NACHOS12.L75/NACHOS12.L75-'+case2+'-MEAN/1d/'+year+'/'
  dirplot='/scratch/cnt0024/hmg2840/albert7a/NACHOS12.L75/PLOTS/NACHOS12.L75-'+case1+'-'+case2+'/python/'

  title="NACHOS12.L75 "+case1+"-"+case2+" "+year

  file1flxT=dirmean1+'NACHOS12.L75-'+case1+'_y'+year+'.1d_flxT.nc'
  file1T=dirmean1+'NACHOS12.L75-'+case1+'_y'+year+'.1d_gridT.nc'
  file1EKE=dirmean1+'NACHOS12.L75-'+case1+'_y'+year+'.1d_EKE.nc'
  file1MXL03=dirmean1+'NACHOS12.L75-'+case1+'_y'+year+'m03.1d_MXL.nc'
  file1MXL09=dirmean1+'NACHOS12.L75-'+case1+'_y'+year+'m09.1d_MXL.nc'
  file1ICE03=dirmean1+'NACHOS12.L75-'+case1+'_y'+year+'m03.1d_icemod3.nc'
  file1ICE09=dirmean1+'NACHOS12.L75-'+case1+'_y'+year+'m09.1d_icemod3.nc'
  file1PSI=dirmean1+'NACHOS12.L75-'+case1+'_y'+year+'.1d_PSI.nc'

  file2flxT=dirmean2+'NACHOS12.L75-'+case2+'_y'+year+'.1d_flxT.nc'
  file2T=dirmean2+'NACHOS12.L75-'+case2+'_y'+year+'.1d_gridT.nc'
  file2EKE=dirmean2+'NACHOS12.L75-'+case2+'_y'+year+'.1d_EKE.nc'
  file2MXL03=dirmean2+'NACHOS12.L75-'+case2+'_y'+year+'m03.1d_MXL.nc'
  file2MXL09=dirmean2+'NACHOS12.L75-'+case2+'_y'+year+'m09.1d_MXL.nc'
  file2ICE03=dirmean2+'NACHOS12.L75-'+case2+'_y'+year+'m03.1d_icemod3.nc'
  file2ICE09=dirmean2+'NACHOS12.L75-'+case2+'_y'+year+'m09.1d_icemod3.nc'
  file2PSI=dirmean2+'NACHOS12.L75-'+case2+'_y'+year+'.1d_PSI.nc'

  ds1T=xr.open_dataset(file1T)
  tem1=ds1T.votemper[0]
  sal1=ds1T.vosaline[0]
  ssh1=ds1T.sossheig[0]
  lat=ds1T.nav_lat
  lon=ds1T.nav_lon

  ds2T=xr.open_dataset(file2T)
  tem2=ds2T.votemper[0]
  sal2=ds2T.vosaline[0]
  ssh2=ds2T.sossheig[0]
  lat2=ds2T.nav_lat
  lon2=ds2T.nav_lon
  
  ds1MXL03=xr.open_dataset(file1MXL03)
  mxl103_rho010=ds1MXL03.somxl010[0]
  mxl103_rho030=ds1MXL03.somxl030[0]
  mxl103_t02=ds1MXL03.somxlt02[0]
  ds1MXL09=xr.open_dataset(file1MXL09)
  mxl109_rho010=ds1MXL09.somxl010[0]
  mxl109_rho030=ds1MXL09.somxl030[0]
  mxl109_t02=ds1MXL09.somxlt02[0]
  
  ds2MXL03=xr.open_dataset(file2MXL03)
  mxl203_rho010=ds2MXL03.somxl010[0]
  mxl203_rho030=ds2MXL03.somxl030[0]
  mxl203_t02=ds2MXL03.somxlt02[0]
  ds2MXL09=xr.open_dataset(file2MXL09)
  mxl209_rho010=ds2MXL09.somxl010[0]
  mxl209_rho030=ds2MXL09.somxl030[0]
  mxl209_t02=ds2MXL09.somxlt02[0]
  
  ds1EKE=xr.open_dataset(file1EKE)
  eke1=ds1EKE.voeke[0,0]
  
  ds2EKE=xr.open_dataset(file2EKE)
  eke2=ds2EKE.voeke[0,0]
  
  ds1PSI=xr.open_dataset(file1PSI)
  psi1=ds1PSI.sobarstf[0]
  
  ds2PSI=xr.open_dataset(file2PSI)
  psi2=ds2PSI.sobarstf[0]
  
  ds1flxT=xr.open_dataset(file1flxT)
  Heat1=ds1flxT.sohefldo[0]
  WaterFlx1=ds1flxT.sowaflup[0]
  WaterDmp1=ds1flxT.sowafld[0]
  
  ds2flxT=xr.open_dataset(file2flxT)
  Heat2=ds2flxT.sohefldo[0]
  WaterFlx2=ds2flxT.sowaflup[0]
  WaterDmp2=ds2flxT.sowafld[0]

  ds1ICE03=xr.open_dataset(file1ICE03)
  iconc103=ds1ICE03.siconc[0]
  ivolu103=ds1ICE03.sivolu[0]
  ds1ICE09=xr.open_dataset(file1ICE09)
  iconc109=ds1ICE09.siconc[0]
  ivolu109=ds1ICE09.sivolu[0]

  ds2ICE03=xr.open_dataset(file2ICE03)
  iconc203=ds2ICE03.siconc[0]
  ivolu203=ds2ICE03.sivolu[0]
  ds2ICE09=xr.open_dataset(file2ICE09)
  iconc209=ds2ICE09.siconc[0]
  ivolu209=ds2ICE09.sivolu[0]
  
  
  def plot_glob(fig,sub,var,vmin,vmax,unit,name,pal):
    ax = fig.add_subplot(sub,projection=ccrs.Orthographic(central_longitude=-30,
                                                    central_latitude=35))
    cmap = plt.get_cmap(pal)
    cmap.set_under(color='grey')
    pcolor=ax.pcolormesh(lon,lat,ma.masked_invalid(var),transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    ax.set_global()
    ax.add_feature(cfeature.LAND,facecolor='grey')
    ax.coastlines()
    cbar=plt.colorbar(pcolor,orientation='vertical',fraction=0.026,pad=0.1)
    cbar.ax.tick_params(labelsize=20)
    ax.set_title(name+' '+unit,size=17)
  
  def plot_glob_diff(fig,sub,var1,var2,vmin,vmax,unit,name,pal):
    ax = fig.add_subplot(sub,projection=ccrs.Orthographic(central_longitude=-30,
                                                    central_latitude=35))
    cmap = plt.get_cmap(pal)
    cmap.set_under(color='grey')
    pcolor=ax.pcolormesh(lon,lat,ma.masked_invalid(var1-var2),transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    ax.set_global()
    ax.add_feature(cfeature.LAND,facecolor='grey')
    ax.coastlines()
    cbar=plt.colorbar(pcolor,orientation='vertical',fraction=0.026,pad=0.1)
    cbar.ax.tick_params(labelsize=20)
    ax.set_title(name+' '+unit,size=17)

  def plot_atl(fig,sub,var,vmin,vmax,unit,name,pal):
    ax = fig.add_subplot(sub,projection=ccrs.PlateCarree(central_longitude=-30))
    cmap = plt.get_cmap(pal)
    ax.set_extent([-100, 50, 0, 70])
    cmap.set_under(color='grey')
    pcolor=ax.pcolormesh(lon,lat,ma.masked_invalid(var),transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    ax.add_feature(cfeature.LAND,facecolor='grey')
    ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

    fig.subplots_adjust(right=0.8)
    ax.text(-0.07, 0.55, 'Latitude (in degree)', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)
    ax.text(0.5, -0.2, 'Longitude (in degree)', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)
    cbar = plt.colorbar(pcolor,orientation='horizontal',shrink=0.75)
    ax.set_title(name+' '+unit,size=17,y=1.08)
    
  def plot_atl_cont(fig,sub,var,unit,name,vmin,vmax,pal):
    ax = fig.add_subplot(sub,projection=ccrs.PlateCarree(central_longitude=-30))
    ax.set_extent([-100, 50, 0, 70])
    cmap = plt.get_cmap(pal)
    cmap.set_under(color='grey')
    pcolor=ax.pcolormesh(lon,lat,ma.masked_invalid(var),transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    pcont=ax.contour(lon,lat,ma.masked_invalid(var),10,colors='k',transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND,facecolor='black')
    ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='grey', alpha=0.5, linestyle='--')

    fig.subplots_adjust(right=0.8)
    ax.text(-0.07, 0.55, 'Latitude (in degree)', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)
    ax.text(0.5, -0.2, 'Longitude (in degree)', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)
    cbar = plt.colorbar(pcolor,orientation='horizontal',shrink=0.75)
    
    ax.set_title(name+' '+unit,size=17,y=1.08)
    
  def plot_natl(fig,sub,var,vmin,vmax,unit,name,pal):
    ax = fig.add_subplot(sub,projection=ccrs.PlateCarree(central_longitude=-30))
    ax.set_extent([-100, 50, 50, 70])
    cmap = plt.get_cmap(pal)
    cmap.set_under(color='grey')
    pcolor=ax.pcolormesh(lon,lat,ma.masked_invalid(var),transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    ax.add_feature(cfeature.LAND,facecolor='grey')
    ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')

    fig.subplots_adjust(right=0.8)
    ax.text(-0.07, 0.55, 'Latitude (in degree)', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)
    ax.text(0.5, -0.4, 'Longitude (in degree)', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)
    cbar = plt.colorbar(pcolor,orientation='horizontal',shrink=0.75)
    ax.set_title(name+' '+unit,size=17,y=1.19)
    

  print dirplot

# Tous les plots

# Tous les plots glob

# Eke, SSH,T et S

  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,10000*eke1,0,2500,'',case1,'jet')
  plot_glob(fig,132,10000*eke2,0,2500,'',case2,'jet')
  plot_glob(fig,133,10000*eke1-10000*eke2,-500,500,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Surf EKE 1e4m2s '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_eke0_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,ssh1,-2.5,-0.7,'',case1,'tab20b')
  plot_glob(fig,132,ssh2,-2.5,-0.7,'',case2,'tab20b')
  plot_glob(fig,133,ssh1-ssh2,-0.2,0.2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 SSH m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_ssh_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,tem1[0],-2,30,'',case1,'jet')
  plot_glob(fig,132,tem2[0],-2,30,'',case2,'jet')
  plot_glob(fig,133,tem1[0]-tem2[0],-2,2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Surf Temperature deg C '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_t0_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,sal1[0],30,40,'',case1,'jet')
  plot_glob(fig,132,sal2[0],30,40,'',case2,'jet')
  plot_glob(fig,133,sal1[0]-sal2[0],-1,1,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Surf Salinity PSU '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_s0_'+year+'.png')
  plt.close()

#MXL 
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,mxl103_rho010,0,1500,'',case1,'jet')
  plot_glob(fig,132,mxl203_rho010,0,1500,'',case2,'jet')
  plot_glob(fig,133,mxl103_rho010-mxl203_rho010,-100,100,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March MXL rho010 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_mxl03_rho010_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,mxl109_rho010,0,200,'',case1,'jet')
  plot_glob(fig,132,mxl209_rho010,0,200,'',case2,'jet')
  plot_glob(fig,133,mxl109_rho010-mxl209_rho010,-50,50,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept MXL rho010 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_mxl09_rho010_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,mxl103_rho030,0,1500,'',case1,'jet')
  plot_glob(fig,132,mxl203_rho030,0,1500,'',case2,'jet')
  plot_glob(fig,133,mxl103_rho030-mxl203_rho030,-100,100,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March MXL rho030 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_mxl03_rho030_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,mxl109_rho030,0,200,'',case1,'jet')
  plot_glob(fig,132,mxl209_rho030,0,200,'',case2,'jet')
  plot_glob(fig,133,mxl109_rho030-mxl209_rho030,-50,50,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept MXL rho030 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_mxl09_rho030_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,mxl103_t02,0,1500,'',case1,'jet')
  plot_glob(fig,132,mxl203_t02,0,1500,'',case2,'jet')
  plot_glob(fig,133,mxl103_t02-mxl203_t02,-100,100,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March MXL t02 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_mxl03_t02_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,mxl109_t02,0,200,'',case1,'jet')
  plot_glob(fig,132,mxl209_t02,0,200,'',case2,'jet')
  plot_glob(fig,133,mxl109_t02-mxl209_t02,-50,50,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept MXL t02 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_mxl09_t02_'+year+'.png')
  plt.close()
 
#flx 
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,Heat1,-400,400,'',case1,'jet')
  plot_glob(fig,132,Heat2,-400,400,'',case2,'jet')
  plot_glob(fig,133,Heat1-Heat2,-100,100,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Net Heat Flux '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_heat_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,86400*WaterFlx1,-9,7,'',case1,'jet')
  plot_glob(fig,132,86400*WaterFlx2,-9,7,'',case2,'jet')
  plot_glob(fig,133,86400*WaterFlx1-86400*WaterFlx2,-2,2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Water Flux '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_wflx_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,86400*WaterDmp1,-7,7,'',case1,'jet')
  plot_glob(fig,132,86400*WaterDmp2,-7,7,'',case2,'jet')
  plot_glob(fig,133,86400*WaterDmp1-86400*WaterDmp2,-2,2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Water Damping '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_wdmp_'+year+'.png')
  plt.close()
  
#ice
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,iconc103,0,1,'',case1,'jet')
  plot_glob(fig,132,iconc203,0,1,'',case2,'jet')
  plot_glob(fig,133,iconc103-iconc203,-0.5,0.5,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March Ice concentration % '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_iceconc03_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,ivolu103,0,4,'',case1,'jet')
  plot_glob(fig,132,ivolu203,0,4,'',case2,'jet')
  plot_glob(fig,133,ivolu103-ivolu203,-0.5,0.5,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March Ice Volume m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_icevolu03_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,iconc109,0,1,'',case1,'jet')
  plot_glob(fig,132,iconc209,0,1,'',case2,'jet')
  plot_glob(fig,133,iconc109-iconc209,-0.2,0.2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept Ice concentration % '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_iceconc09_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_glob(fig,131,ivolu109,0,4,'',case1,'jet')
  plot_glob(fig,132,ivolu209,0,4,'',case2,'jet')
  plot_glob(fig,133,ivolu109-ivolu209,-0.5,0.5,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept Ice Volume m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_glob_icevolu09_'+year+'.png')
  plt.close()
  
# Tous les plots Atlantique

# T & S

  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,tem1[0],-2,30,'',case1,'jet')
  plot_atl(fig,132,tem2[0],-2,30,'',case2,'jet')
  plot_atl(fig,133,tem1[0]-tem2[0],-2,2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Surf Temperature deg C '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_t0_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,sal1[0],30,40,'',case1,'jet')
  plot_atl(fig,132,sal2[0],30,40,'',case2,'jet')
  plot_atl(fig,133,sal1[0]-sal2[0],-1,1,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Surf Salinity PSU '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_s0_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,tem1[30],-2,30,'',case1,'jet')
  plot_atl(fig,132,tem2[30],-2,30,'',case2,'jet')
  plot_atl(fig,133,tem1[30]-tem2[30],-2,2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 200m Temperature deg C '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_t200_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,sal1[30],30,40,'',case1,'jet')
  plot_atl(fig,132,sal2[30],30,40,'',case2,'jet')
  plot_atl(fig,133,sal1[30]-sal2[30],-2,2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 200m Salinity PSU '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_s200_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,tem1[46],-2,30,'',case1,'jet')
  plot_atl(fig,132,tem2[46],-2,30,'',case2,'jet')
  plot_atl(fig,133,tem1[46]-tem2[46],-2,2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 1000m Temperature deg C '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_t1000_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,sal1[46],30,40,'',case1,'jet')
  plot_atl(fig,132,sal2[46],30,40,'',case2,'jet')
  plot_atl(fig,133,sal1[46]-sal2[46],-1,1,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 1000m Salinity PSU '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_s1000_'+year+'.png')
  plt.close()

#PSI  
  fig = plt.figure(figsize=(22,7))
  plot_atl_cont(fig,131,1e-7*psi1,'',case1,-4,4,'Blues')
  plot_atl_cont(fig,132,1e-7*psi2,'',case2,-4,4,'Blues')
  plot_atl_cont(fig,133,1e-7*psi1-1e-7*psi2,'',case1+'-'+case2,-1,1,'bwr')
  fig.suptitle('NACHOS12.L75 Stream function '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_psi_'+year+'.png')
  plt.close()

#flx  
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,Heat1,-400,400,'',case1,'jet')
  plot_atl(fig,132,Heat2,-400,400,'',case2,'jet')
  plot_atl(fig,133,Heat1-Heat2,-100,100,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Net Heat Flux '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_heat_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,86400*WaterFlx1,-9,7,'',case1,'jet')
  plot_atl(fig,132,86400*WaterFlx2,-9,7,'',case2,'jet')
  plot_atl(fig,133,86400*WaterFlx1-86400*WaterFlx2,-2,2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Water Flux '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_wflx_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,86400*WaterDmp1,-7,7,'',case1,'jet')
  plot_atl(fig,132,86400*WaterDmp2,-7,7,'',case2,'jet')
  plot_atl(fig,133,86400*WaterDmp1-86400*WaterDmp2,-2,2,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Water Damping '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_wdmp_'+year+'.png')
  plt.close()
  
#MXL
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,mxl103_rho010,0,1500,'',case1,'jet')
  plot_atl(fig,132,mxl203_rho010,0,1500,'',case2,'jet')
  plot_atl(fig,133,mxl103_rho010-mxl203_rho010,-100,100,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March MXL rho010 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_mxl03_rho010_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_atl(fig,131,mxl109_rho010,0,200,'',case1,'jet')
  plot_atl(fig,132,mxl209_rho010,0,200,'',case2,'jet')
  plot_atl(fig,133,mxl109_rho010-mxl209_rho010,-50,50,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept MXL rho010 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_atl_mxl09_rho010_'+year+'.png')
  plt.close()
 
#Tous les plots N Atlantique
 
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,iconc103,0,1,'',case1,'jet')
  plot_natl(fig,132,iconc203,0,1,'',case2,'jet')
  plot_natl(fig,133,iconc103-iconc203,-0.5,0.5,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March Ice concentration % '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_iceconc03_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,ivolu103,0,4,'',case1,'jet')
  plot_natl(fig,132,ivolu203,0,4,'',case2,'jet')
  plot_natl(fig,133,ivolu103-ivolu203,-0.5,0.5,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March Ice Volume m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_icevolu03_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,iconc109,0,1,'',case1,'jet')
  plot_natl(fig,132,iconc209,0,1,'',case2,'jet')
  plot_natl(fig,133,iconc109-iconc209,-0.5,0.5,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept Ice concentration % '+case1+'-'+case2+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_iceconc09_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,ivolu109,0,4,'',case1,'jet')
  plot_natl(fig,132,ivolu209,0,4,'',case2,'jet')
  plot_natl(fig,133,ivolu109-ivolu209,-0.5,0.5,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept Ice Volume m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_icevolu09_'+year+'.png')
  plt.close()
  
#MXL
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,mxl103_rho010,0,1500,'',case1,'jet')
  plot_natl(fig,132,mxl203_rho010,0,1500,'',case2,'jet')
  plot_natl(fig,133,mxl103_rho010-mxl203_rho010,-100,100,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March MXL rho010 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_mxl03_rho010_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,mxl109_rho010,0,200,'',case1,'jet')
  plot_natl(fig,132,mxl209_rho010,0,200,'',case2,'jet')
  plot_natl(fig,133,mxl109_rho010-mxl209_rho010,-50,50,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept MXL rho010 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_mxl09_rho010_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,mxl103_rho030,0,1500,'',case1,'jet')
  plot_natl(fig,132,mxl203_rho030,0,1500,'',case2,'jet')
  plot_natl(fig,133,mxl103_rho030-mxl203_rho030,-100,100,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March MXL rho030 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_mxl03_rho030_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,mxl109_rho030,0,200,'',case1,'jet')
  plot_natl(fig,132,mxl209_rho030,0,200,'',case2,'jet')
  plot_natl(fig,133,mxl109_rho030-mxl209_rho030,-50,50,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept MXL rho030 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_mxl09_rho030_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,mxl103_t02,0,1500,'',case1,'jet')
  plot_natl(fig,132,mxl203_t02,0,1500,'',case2,'jet')
  plot_natl(fig,133,mxl103_t02-mxl203_t02,-100,100,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 March MXL t02 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_mxl03_t02_'+year+'.png')
  plt.close()
  fig = plt.figure(figsize=(22,7))
  plot_natl(fig,131,mxl109_t02,0,200,'',case1,'jet')
  plot_natl(fig,132,mxl209_t02,0,200,'',case2,'jet')
  plot_natl(fig,133,mxl109_t02-mxl209_t02,-50,50,'',case1+'-'+case2,'bwr')
  fig.suptitle('NACHOS12.L75 Sept MXL t02 m '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case1+'-'+case2+'_natl_mxl09_t02_'+year+'.png')
  plt.close()

def script_parser():
    """Customized parser.
    """
    from optparse import OptionParser
    usage = "usage: %prog CASE1 CASE2 year"
    parser = OptionParser(usage=usage)
    return parser


def main():
    parser = script_parser()
    (options, args) = parser.parse_args()
    if len(args) < 3: # print the help message if number of args is not 2.
        parser.print_help()
        sys.exit()
    optdic = vars(options)
    if len(args) == 3:
      case1 = args[0]
      case2 = args[1]
      year = args[2]
      all_plots(case1,case2,year,**optdic)

if __name__ == '__main__':
    sys.exit(main() or 0)




