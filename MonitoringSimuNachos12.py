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


def all_plots(case,year,**kwargs):

  dirmean='/scratch/cnt0024/hmg2840/albert7a/NACHOS12.L75/NACHOS12.L75-'+case+'-MEAN/1d/'+year+'/'
  dirplot='/scratch/cnt0024/hmg2840/albert7a/NACHOS12.L75/PLOTS/NACHOS12.L75-'+case+'/python/'

  title="NACHOS12.L75-"+case+" "+year

  fileflxT=dirmean+'NACHOS12.L75-'+case+'_y'+year+'.1d_flxT.nc'
  fileT=dirmean+'NACHOS12.L75-'+case+'_y'+year+'.1d_gridT.nc'
  fileEKE=dirmean+'NACHOS12.L75-'+case+'_y'+year+'.1d_EKE.nc'
  fileMXL03=dirmean+'NACHOS12.L75-'+case+'_y'+year+'m03.1d_MXL.nc'
  fileMXL09=dirmean+'NACHOS12.L75-'+case+'_y'+year+'m09.1d_MXL.nc'
  fileICE03=dirmean+'NACHOS12.L75-'+case+'_y'+year+'m03.1d_icemod3.nc'
  fileICE09=dirmean+'NACHOS12.L75-'+case+'_y'+year+'m09.1d_icemod3.nc'
  filePSI=dirmean+'NACHOS12.L75-'+case+'_y'+year+'.1d_PSI.nc'

  dsT=xr.open_dataset(fileT)
  tem=dsT.votemper[0]
  sal=dsT.vosaline[0]
  ssh=dsT.sossheig[0]
  lat=dsT.nav_lat
  lon=dsT.nav_lon
  
  dsMXL03=xr.open_dataset(fileMXL03)
  mxl03_rho010=dsMXL03.somxl010[0]
  mxl03_rho030=dsMXL03.somxl030[0]
  mxl03_t02=dsMXL03.somxlt02[0]
  dsMXL09=xr.open_dataset(fileMXL09)
  mxl09_rho010=dsMXL09.somxl010[0]
  mxl09_rho030=dsMXL09.somxl030[0]
  mxl09_t02=dsMXL09.somxlt02[0]
  
  dsEKE=xr.open_dataset(fileEKE)
  eke=dsEKE.voeke[0,0]
  
  dsPSI=xr.open_dataset(filePSI)
  psi=dsPSI.sobarstf[0]
  
  dsflxT=xr.open_dataset(fileflxT)
  Heat=dsflxT.sohefldo[0]
  WaterFlx=dsflxT.sowaflup[0]
  WaterDmp=dsflxT.sowafld[0]

  dsICE03=xr.open_dataset(fileICE03)
  iconc03=dsICE03.siconc[0]
  ivolu03=dsICE03.sivolu[0]
  dsICE09=xr.open_dataset(fileICE09)
  iconc09=dsICE09.siconc[0]
  ivolu09=dsICE09.sivolu[0]
  
  
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
    ax.text(0.5, -0.1, 'Longitude (in degree)', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)
    cbar = plt.colorbar(pcolor,orientation='horizontal',shrink=0.75)
    ax.set_title(name+' '+unit,size=17,y=1.08)
    
  def plot_atl_cont(fig,sub,var,unit,name,vmin,vmax,pal):
    ax = fig.add_subplot(sub,projection=ccrs.PlateCarree(central_longitude=-30))
    ax.set_extent([-100, 50, 0, 70])
    cmap = plt.get_cmap(pal)
    cmap.set_under(color='grey')
    pcolor=ax.pcolormesh(lon,lat,ma.masked_invalid(var),transform=ccrs.PlateCarree(),cmap=plt.cm.Blues,vmin=vmin,vmax=vmax)
    pcont=ax.contour(lon,lat,ma.masked_invalid(var),10,colors='k',transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND,facecolor='black')
    ax.coastlines()
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='grey', alpha=0.5, linestyle='--')

    fig.subplots_adjust(right=0.8)
    ax.text(-0.07, 0.55, 'Latitude (in degree)', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)
    ax.text(0.5, -0.1, 'Longitude (in degree)', va='bottom', ha='center',
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
    ax.text(0.5, -0.1, 'Longitude (in degree)', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)
    cbar = plt.colorbar(pcolor,orientation='horizontal',shrink=0.75)
    ax.set_title(name+' '+unit,size=17,y=1.09)
    

  print dirplot

# Tous les plots

# Tous les plots glob

# Eke, SSH,T et S

  fig = plt.figure(figsize=(15,15))
  plot_glob(fig,221,10000*eke,0,2500,'1e4m2s','Surf EKE','jet')
  plot_glob(fig,222,ssh,-2.5,-0.7,'m','SSH','tab20b')
  plot_glob(fig,223,tem[0],-2,30,'deg C','Surf Temperature','jet')
  plot_glob(fig,224,sal[0],30,40,'PSU','Surf Salinity','jet')
  fig.suptitle('NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_glob_eke0-ssh-t0-s0_'+year+'.png')

#MXL 
  fig = plt.figure(figsize=(30,20))
  plot_glob(fig,231,mxl03_rho010,0,1500,'m','March MXL rho010','jet')
  plot_glob(fig,234,mxl09_rho010,0,200,'m','Sept MXL rho010','jet')
  plot_glob(fig,232,mxl03_rho030,0,1500,'m','March MXL rho030','jet')
  plot_glob(fig,235,mxl09_rho030,0,200,'m','Sept MXL rho030','jet')
  plot_glob(fig,233,mxl03_t02,0,1500,'m','March MXL t02','jet')
  plot_glob(fig,236,mxl09_t02,0,200,'m','Sept MXL t02','jet')
  fig.suptitle('NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_glob_mxl_'+year+'.png')
 
#flx 
  fig = plt.figure(figsize=(30,7))
  plot_glob(fig,131,Heat,-400,400,'','Net Heat Flux','jet')
  plot_glob(fig,132,86400*WaterFlx,-9,7,'','Water Flux','jet')
  plot_glob(fig,133,86400*WaterDmp,-7,7,'','Water Damping','jet')
  fig.suptitle('NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_glob_flxt_'+year+'.png')
  
#ice
  fig = plt.figure(figsize=(15,15))
  plot_glob(fig,221,iconc03,0,1,'%','March Ice concentration','jet')
  plot_glob(fig,222,ivolu03,0,4,'m','March Ice Volume','jet')
  plot_glob(fig,223,iconc09,0,1,'%','Sept Ice concentration','jet')
  plot_glob(fig,224,ivolu09,0,4,'m','Sept Ice Volume','jet')
  fig.suptitle('NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_glob_ice_'+year+'.png')
  
# Tous les plots Atlantique

# T & S

  fig = plt.figure(figsize=(30,15))
  plot_atl(fig,231,tem[0],-2,30,'deg C','Surf Temperature','jet')
  plot_atl(fig,234,sal[0],30,40,'PSU','Surf Salinity','jet')
  plot_atl(fig,232,tem[30],-2,30,'deg C','200m Temperature','jet')
  plot_atl(fig,235,sal[30],30,40,'PSU','200m Salinity','jet')
  plot_atl(fig,233,tem[46],-2,30,'deg C','1000m Temperature','jet')
  plot_atl(fig,236,sal[46],30,40,'PSU','1000m Salinity','jet')
  fig.suptitle('NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_atl_t-s-0-200-1000_'+year+'.png')

#PSI  
  fig = plt.figure(figsize=(10,7))
  plot_atl_cont(fig,111,1e-7*psi,'','',-4,4,'jet')
  fig.suptitle('Stream function NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_atl_psi_'+year+'.png')

#flx  
  fig = plt.figure(figsize=(30,7))
  plot_atl(fig,131,Heat,-400,400,'','Net Heat Flux','jet')
  plot_atl(fig,132,86400*WaterFlx,-9,7,'','Water Flux','jet')
  plot_atl(fig,133,86400*WaterDmp,-7,7,'','Water Damping','jet')
  fig.suptitle('NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_atl_flxt_'+year+'.png')
  
#MXL
  fig = plt.figure(figsize=(30,10))
  plot_atl(fig,121,mxl03_rho010,0,1500,'m','March MXL rho010','jet')
  plot_atl(fig,122,mxl09_rho010,0,200,'m','Sept MXL rho010','jet')
  fig.suptitle('NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_atl_mxl_'+year+'.png')
 
#Tous les plots N Atlantique
 
  fig = plt.figure(figsize=(30,15))
  plot_natl(fig,221,iconc03,0,1,'%','March Ice concentration','jet')
  plot_natl(fig,222,ivolu03,0,4,'m','March Ice Volume','jet')
  plot_natl(fig,223,iconc09,0,1,'%','Sept Ice concentration','jet')
  plot_natl(fig,224,ivolu09,0,4,'m','Sept Ice Volume','jet')
  fig.suptitle('NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_natl_ice_'+year+'.png')
  
#MXL
  fig = plt.figure(figsize=(30,25))
  plot_natl(fig,321,mxl03_rho010,0,1500,'m','March MXL rho010','jet')
  plot_natl(fig,322,mxl09_rho010,0,200,'m','Sept MXL rho010','jet')
  plot_natl(fig,323,mxl03_rho030,0,1500,'m','March MXL rho030','jet')
  plot_natl(fig,324,mxl09_rho030,0,200,'m','Sept MXL rho030','jet')
  plot_natl(fig,325,mxl03_t02,0,1500,'m','March MXL t02','jet')
  plot_natl(fig,326,mxl09_t02,0,200,'m','Sept MXL t02','jet')
  fig.suptitle('NACHOS12.L75-'+case+' '+year, fontsize=25)
  plt.savefig(dirplot+'NACHOS12.L75-'+case+'_natl_mxl_'+year+'.png')

def script_parser():
    """Customized parser.
    """
    from optparse import OptionParser
    usage = "usage: %prog CASE year"
    parser = OptionParser(usage=usage)
    return parser


def main():
    parser = script_parser()
    (options, args) = parser.parse_args()
    if len(args) < 2: # print the help message if number of args is not 2.
        parser.print_help()
        sys.exit()
    optdic = vars(options)
    if len(args) == 2:
      case = args[0]
      year = args[1]
      all_plots(case,year,**optdic)

if __name__ == '__main__':
    sys.exit(main() or 0)




