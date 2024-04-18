###################################################################### CHJ #####
## Name		: plot_fv3lam_mrms.py
## Language	: Python 3.7
## Usage	: Plot MRMS reflectivity data for comparison
## Input files  : MRMS_MergedReflectivityQCComposite_00.50.XXX.grib2
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/07/24: Chan-Hoo Jeon : Preliminary version
## V001: 2020/07/28: Chan-Hoo Jeon : Add 3D plot option
## V002: 2020/07/29: Chan-Hoo Jeon : Add cross-sectional plot option
## V003: 2021/03/05: Chan-Hoo Jeon : Simplify the script
###################################################################### CHJ #####

import os,sys,time
import pygrib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import itertools
from cartopy.mpl.patch import geos_to_path
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection,LineCollection
import matplotlib.gridspec as gridspec

# HPC machine ('hera','orion')
machine='hera'

print(' You are on', machine)

# Path to Natural Earth Data-set for background plot
if machine=='hera':
    path_NE='/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/NaturalEarth'
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig/"
elif machine=='orion':
    path_NE='/home/chjeon/tools/NaturalEarth'
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
else:
    sys.exit('ERROR: path to Natural Earth Data or output dir. is not set !!!')

cartopy.config['data_dir']=path_NE
os.environ["CARTOPY_USER_BACKGROUNDS"]=path_NE+'/raster_files'

plt.switch_backend('agg')

# Global variables ======================================== CHJ =====
# ..... Case-dependent input :: should be changed case-by-case .....
# ******
# INPUT
# ******

s_date='2019070100'
s_time='03'

# Path to the directory where the input file is located.
dnm_data='/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/00_DATA/MRMS/mrms_'+s_date+'/new_grib2/'

# Input file name
fnm_in_com='QCComposite_00.50_'+s_date[0:8]+s_time+'00_new.grib2'
fnm_in_3d='ref3D_'+s_date[0:8]+s_time+'00_new.grib2'

# Flag for composite reflectivity
plt_com='yes'
# Flag for cross-sectional reflectivity
plt_crx='yes'
# target lon,lat for cross-section/3d (in degrees)
crx_lon=-99
crx_lat=43.7
# Distance from the target point (in degrees)
crx_dist=4
# Flag for 3d reflectivity
plt_3d='yes'
# plotting layers' index numbers (only for 3d):
plt_hght_i=[0,2,4,6]

# Title and output file name:
out_title_com='MRMS::Composite Reflectivity::'+s_date+'/'+s_time
out_fname_com='fv3lam_mrms_comRefl_'+s_date+'_'+s_time
out_fname_crx='fv3lam_mrms_xzRefl_'+s_date+'_'+s_time
out_fname_3d='fv3lam_mrms_3dRefl_'+s_date+'_'+s_time

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'

# high-resolution background image ('on', 'off')
back_img='off'



# Main part (will be called at the end) ==================== CHJ =====
def main():
# ========================================================== CHJ =====

   if plt_com=='yes':
      plot_composite()
   if plt_crx=='yes':
      plot_xsect()
   if plt_3d=='yes':
      plot_3d()



# Plot: composite reflectivity ============================= CHJ =====
def plot_composite():
# ========================================================== CHJ =====

    print(' ===== MRMS:: Composite Reflectivity ===============================')

    # open the data file
    fname=os.path.join(dnm_data,fnm_in_com)
    try: grbs=pygrib.open(fname)
    except: raise Exception('Could NOT find the file',fname)

    for grb in grbs:
#        print(grb.name)
        print(grb.typeOfLevel)
        print(grb.level)
#        print(grb.validDate)
#        print(grb.analDate)
#        print(grb.Nx)
#        print(grb.Ny)
#        print(grb.shortName)
#        print(grb.missingValue)

    grbv=grbs.select(typeOfLevel="heightAboveSea",level=500)[0]
    sval=grbv.values

    #print(type(sval))
    # replace values for missing data
    sval[sval<0]=0.0

    Nxd=grbv.Nx
    Nyd=grbv.Ny
    lat,lon=grbv.latlons()

    nm_svar='Composite reflectivity (dBZ)'
    cs_cmap,cs_norm=new_cmap()
    lb_ext='max'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=2
    cmap_range='designed'

    # Highest and lowest longitudes and latitudes for plot extent
#        lon_min=-128.39
#        lon_max=-66.62
#        lat_min=25.12
#        lat_max=49.23
    lon_min=np.min(lon)
    lon_max=np.max(lon)
    lat_min=np.min(lat)
    lat_max=np.max(lat)

    print(' lon_min=',lon_min,', lon_max=',lon_max)
    print(' lat_min=',lat_min,', lat_max=',lat_max)

    esp=1    
    extent=[lon_min-esp,lon_max+esp,lat_min-esp,lat_max+esp]
#    extent=[lon_min,lon_max,lat_min,lat_max]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])
 
    # Max and Min of the field
    fmax=np.nanmax(sval)
    fmin=np.nanmin(sval)

    print(' fld_min=',fmin)
    print(' fld_max=',fmax)

    # Make the colormap range symmetry
    print(' cmap range=',cmap_range)
    if cmap_range=='symmetry': 
        tmp_cmp=max(abs(fmax),abs(fmin))
        cs_min=round(-tmp_cmp,n_rnd)
        cs_max=round(tmp_cmp,n_rnd)
    elif cmap_range=='round':
        cs_min=round(fmin,n_rnd)
        cs_max=round(fmax,n_rnd)
    elif cmap_range=='real':
        cs_min=fmin
        cs_max=fmax
    elif cmap_range=='fixed':
        cs_min=-10.0
        cs_max=10.0
    elif cmap_range=='designed':
        cs_min=5
        cs_max=None
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_min=',cs_min)
    print(' cs_max=',cs_max)
 
    # Plot field
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.PlateCarree(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title_com,fontsize=9)
    cs=ax.pcolormesh(lon,lat,sval,cmap=cs_cmap,rasterized=True,
            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree(),norm=cs_norm)

    gl=ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                    linewidth=0.2,color='chocolate',alpha=0.3,linestyle='-')
    gl.xlabels_top=False
    gl.ylabels_right=False
    gl.xlocator=mticker.FixedLocator([-130,-120,-110,-100,-90,-80,-70,-60])
    gl.ylocator=mticker.FixedLocator([20,30,40,50,60])
    gl.xformatter=LONGITUDE_FORMATTER
    gl.yformatter=LATITUDE_FORMATTER
    gl.xlabel_style={'size':5,'color':'black'}
    gl.ylabel_style={'size':5,'color':'black'}

    if cmap_range=='designed':
        cs.cmap.set_under('white',alpha=0.)
        cs.cmap.set_over('black')
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(nm_svar,fontsize=8)

    # Output figure
    ndpi=300
    out_file(out_fname_com,ndpi) 




# Plot: Cross-sectional reflectivity ======================= CHJ =====
def plot_xsect():
# ========================================================== CHJ =====
    global extent

    print(' ===== MRMS:: Composite Reflectivity ===============================')

    # open the data file
    fname=os.path.join(dnm_data,fnm_in_com)
    try: grbs=pygrib.open(fname)
    except: raise Exception('Could NOT find the file',fname)

    for grb in grbs:
#        print(grb.name)
        print(grb.typeOfLevel)
        print(grb.level)
#        print(grb.validDate)
#        print(grb.analDate)
#        print(grb.Nx)
#        print(grb.Ny)
#        print(grb.shortName)
#        print(grb.missingValue)

    grbv=grbs.select(typeOfLevel="heightAboveSea",level=500)[0]
    sval_com=grbv.values

    # replace values for missing data
    sval_com[sval_com<0]=0.0

    Njx=grbv.Nx
    Niy=grbv.Ny
    lat,lon=grbv.latlons()


    print(' ===== MRMS:: Cross-sectional Reflectivity ===============================')

    # open the data file
    fname=os.path.join(dnm_data,fnm_in_3d)
    try: grbs=pygrib.open(fname)
    except: raise Exception('Could NOT find the file',fname)

    for grb in grbs:
#        print(grb.name)
        print(grb.typeOfLevel)
        print(grb.level)
#        print(grb.validDate)
#        print(grb.analDate)
#        print(grb.Nx)
#        print(grb.Ny)
#        print(grb.shortName)
#        print(grb.missingValue)


    hghts=np.concatenate([np.arange(500,3000,250),np.arange(3000,9000,500),np.arange(9000,13000,1000)])
    hght_km=hghts/1000
    print(hghts)
    print(hght_km)
    Nkz=len(hghts)
    print(Niy,Njx,Nkz)

    sval_3d=np.zeros((Niy,Njx,Nkz))

#    i_tmp=0
    for i in range(0,Nkz):
        print(i,hghts[i])
        grbv=grbs.select(typeOfLevel="heightAboveSea",level=hghts[i])[0]
        sval=grbv.values
        # replace values for missing data
        sval[sval<0]=0
        sval_3d[:,:,i]=sval
            
#        if i_tmp==0:
#            Nxd=grbv.Nx
#            Nyd=grbv.Ny
#            lat,lon=grbv.latlons()
#            i_tmp=1

    print(sval_3d.shape)

    # Find the nearest longitude/latitude from the target point
    # lon: 0~360 -> -180~180
    lon_pm=lon-360.0
    tmp=np.abs(lon_pm[0,:]-crx_lon)
    idx_lon=np.argwhere(tmp==np.min(tmp))
    print(' Nearest longitude at',idx_lon)
    tmp=np.abs(lat[:,0]-crx_lat)
    idx_lat=np.argwhere(tmp==np.min(tmp))
    print(' Nearest latitude at',idx_lat)

    nst_tg_lon=np.squeeze(lon_pm[0,idx_lon])
    nst_tg_lat=np.squeeze(lat[idx_lat,0]) 
    print(' Target (lon,lat)=',crx_lon,crx_lat)
    print(' Nearest (lon,lat)=',nst_tg_lon,nst_tg_lat)

    crx_xz_s_lon=np.squeeze(lon_pm[0,idx_lon])-crx_dist
    crx_xz_s_lat=np.squeeze(lat[0,idx_lon])
    crx_xz_e_lon=np.squeeze(lon_pm[0,idx_lon])+crx_dist
    crx_xz_e_lat=np.squeeze(lat[0,idx_lon])
    crx_yz_s_lon=np.squeeze(lon_pm[idx_lat,0])
    crx_yz_s_lat=np.squeeze(lat[idx_lat,0])-crx_dist
    crx_yz_e_lon=np.squeeze(lon_pm[idx_lat,0])
    crx_yz_e_lat=np.squeeze(lat[idx_lat,0])+crx_dist

    print(' xz line from',crx_xz_s_lon,crx_xz_s_lat,' to',crx_xz_e_lon,crx_xz_e_lat)
    print(' yz line from',crx_yz_s_lon,crx_yz_s_lat,' to',crx_yz_e_lon,crx_yz_e_lat)


    crx2d_lon=np.squeeze(sval_3d[idx_lat,:,:])
    crx2d_lat=np.squeeze(sval_3d[:,idx_lon,:])

    print(crx2d_lon.shape)
    print(crx2d_lat.shape)

    lon2d=np.squeeze(lon_pm[idx_lat,:])
    lat2d=np.squeeze(lat[:,idx_lon])

    print(lon2d.shape)
    print(lat2d.shape)


    crx_xz_lon,crx_xz_hgt=np.meshgrid(lon2d,hght_km)
    crx_yz_lat,crx_yz_hgt=np.meshgrid(lat2d,hght_km)

    print(crx_xz_lon.shape)
    print(crx_xz_hgt.shape)
    print(crx_yz_lat.shape)
    print(crx_yz_hgt.shape)

    crx2d_lon=np.transpose(crx2d_lon)
    crx2d_lat=np.transpose(crx2d_lat)

    print(crx2d_lon.shape)
    print(crx2d_lat.shape)
 
    nm_svar='Reflectivity (dBZ)'
    cs_cmap,cs_norm=new_cmap()
    lb_ext='max'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=2
    cmap_range='designed'

    # Highest and lowest longitudes and latitudes for plot extent
    #lon_min=-128.39
    #lon_max=-66.62
    #lat_min=25.12
    #lat_max=49.23
    lon_min=np.min(lon)
    lon_max=np.max(lon)
    lat_min=np.min(lat)
    lat_max=np.max(lat)

    print(' lon_min=',lon_min,', lon_max=',lon_max)
    print(' lat_min=',lat_min,', lat_max=',lat_max)

    esp=1    
    extent=[lon_min-esp,lon_max+esp,lat_min-esp,lat_max+esp]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])
 
    # Max and Min of the field
    fmax=np.nanmax(sval_3d)
    fmin=np.nanmin(sval_3d)

    print(' fld_min=',fmin)
    print(' fld_max=',fmax)

    # Make the colormap range symmetry
    print(' cmap range=',cmap_range)
    if cmap_range=='symmetry': 
        tmp_cmp=max(abs(fmax),abs(fmin))
        cs_min=round(-tmp_cmp,n_rnd)
        cs_max=round(tmp_cmp,n_rnd)
    elif cmap_range=='round':
        cs_min=round(fmin,n_rnd)
        cs_max=round(fmax,n_rnd)
    elif cmap_range=='real':
        cs_min=fmin
        cs_max=fmax
    elif cmap_range=='fixed':
        cs_min=-10.0
        cs_max=10.0
    elif cmap_range=='designed':
        cs_min=5
        cs_max=None
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_min=',cs_min)
    print(' cs_max=',cs_max)
 

    fig=plt.figure(figsize=(7,12))   #(width,height)
    spec_fig=fig.add_gridspec(3,1,hspace=0.25)

# subplot 1: composite reflectivity + target point lines
    ax1=fig.add_subplot(spec_fig[0,0],projection=ccrs.PlateCarree(c_lon))

    ax1.set_extent(extent,ccrs.PlateCarree())
    back_plot(ax1)

    ax1.set_title(out_title_com,fontsize=9)
    cs=ax1.pcolormesh(lon,lat,sval_com,cmap=cs_cmap,rasterized=True,
            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree(),norm=cs_norm)

    gl=ax1.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                    linewidth=0.2,color='chocolate',alpha=0.3,linestyle='-')
    gl.xlabels_top=False
    gl.ylabels_right=False
    gl.xlocator=mticker.FixedLocator([-130,-120,-110,-100,-90,-80,-70,-60])
    gl.ylocator=mticker.FixedLocator([20,30,40,50,60])
    gl.xformatter=LONGITUDE_FORMATTER
    gl.yformatter=LATITUDE_FORMATTER
    gl.xlabel_style={'size':7,'color':'black'}
    gl.ylabel_style={'size':7,'color':'black'}

    ax1.plot(nst_tg_lon,nst_tg_lat,'mx',markersize=5,transform=ccrs.Geodetic())
    ax1.text(nst_tg_lon+0.5,nst_tg_lat-1.5,'Target Point',fontsize=10,color='fuchsia',transform=ccrs.Geodetic())


    if cmap_range=='designed':
        cs.cmap.set_under('white',alpha=0.)
        cs.cmap.set_over('black')
    divider=make_axes_locatable(ax1)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label(nm_svar,fontsize=7)


# subplot 2: cross-section in the longitudinal direction
    ax2=fig.add_subplot(spec_fig[1,0])

    ax2.contourf(crx_xz_lon,crx_xz_hgt,crx2d_lon,np.linspace(5,70,14),cmap=cs_cmap,norm=cs_norm,
                 vmin=cs_min,vmax=cs_max,alpha=1.0)
    ax2.set_xlim(crx_xz_s_lon,crx_xz_e_lon)
    ax2.set_ylim(0,12)
    ax2.set_title('Longitudinal cross-section (X-Z) centered at the target point',fontsize=9)
    divider=make_axes_locatable(ax2)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label(nm_svar,fontsize=7)

    ax2.grid(color='chocolate',linestyle='-',linewidth=0.2,alpha=0.3)
    lb_font_sz=7
    ax2.set_xlabel('Longitude ($^\circ$)',fontsize=lb_font_sz)
    ax2.set_ylabel('Height (km)',fontsize=lb_font_sz)
    ax2.xaxis.set_tick_params(labelsize=lb_font_sz)
    ax2.yaxis.set_tick_params(labelsize=lb_font_sz)


# subplot 3: cross-section in the latitudinal direction
    ax3=fig.add_subplot(spec_fig[2,0],sharey=ax2)

    ax3.contourf(crx_yz_lat,crx_yz_hgt,crx2d_lat,np.linspace(5,70,14),cmap=cs_cmap,norm=cs_norm,
                 vmin=cs_min,vmax=cs_max,alpha=1.0)
    ax3.set_xlim(crx_yz_s_lat,crx_yz_e_lat)
    ax3.set_ylim(0,12)
    ax3.set_title('Latitudinal cross-section (Y-Z) centered at the target point',fontsize=9)
    divider=make_axes_locatable(ax3)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label(nm_svar,fontsize=7)

    ax3.grid(color='chocolate',linestyle='-',linewidth=0.2,alpha=0.3)
    ax3.set_xlabel('Latitude ($^\circ$)',fontsize=lb_font_sz)
    ax3.set_ylabel('Height (km)',fontsize=lb_font_sz)
    ax3.xaxis.set_tick_params(labelsize=lb_font_sz)
    ax3.yaxis.set_tick_params(labelsize=lb_font_sz)


    # Output figure
    ndpi=300
    out_file(out_fname_crx,ndpi)




# Plot: 3-D reflectivity =================================== CHJ =====
def plot_3d():
# ========================================================== CHJ =====
    global extent,max_zaxis

    print(' ===== MRMS:: 3-D Reflectivity ===============================')

    # open the data file
    fname=os.path.join(dnm_data,fnm_in_3d)
    try: grbs=pygrib.open(fname)
    except: raise Exception('Could NOT find the file',fname)

    for grb in grbs:
#        print(grb.name)
        print(grb.typeOfLevel)
        print(grb.level)
#        print(grb.validDate)
#        print(grb.analDate)
#        print(grb.Nx)
#        print(grb.Ny)
#        print(grb.shortName)
#        print(grb.missingValue)

    Njx=grb.Nx
    Niy=grb.Ny

    hghts=np.concatenate([np.arange(500,3000,250),np.arange(3000,9000,500),np.arange(9000,13000,1000)])
    hght_km=hghts/1000
    print(hghts)
    print(hght_km)
    Nkz=len(hghts)
    print(Niy,Njx,Nkz)

    sval_3d=np.zeros((Niy,Njx,Nkz))

#    plt_hght_i=[0,1,2,3]
    max_zaxis=hght_km[plt_hght_i[-1]]

    print('max_hght=',max_zaxis)

    i_tmp=0
    for i in plt_hght_i:
        print(i,hghts[i])
        grbv=grbs.select(typeOfLevel="heightAboveSea",level=hghts[i])[0]
        sval=grbv.values
        # replace values for missing data
        sval[sval<0]=0
        sval_3d[:,:,i]=sval
            
        if i_tmp==0:
            Nxd=grbv.Nx
            Nyd=grbv.Ny
            lat,lon=grbv.latlons()
            i_tmp=1

    print(sval_3d.shape)

    # Find the nearest longitude/latitude from the target point
    # lon: 0~360 -> -180~180
    lon_pm=lon-360.0
    tmp=np.abs(lon_pm[0,:]-crx_lon)
    idx_lon=np.argwhere(tmp==np.min(tmp))
    print(' Nearest longitude at',idx_lon)
    tmp=np.abs(lat[:,0]-crx_lat)
    idx_lat=np.argwhere(tmp==np.min(tmp))
    print(' Nearest latitude at',idx_lat)

    tmp=np.abs(lon_pm[0,:]-(crx_lon-2*crx_dist))
    idx_lon_s=np.squeeze(np.argwhere(tmp==np.min(tmp)))
    print(' Longitude limit from',idx_lon_s)
    tmp=np.abs(lon_pm[0,:]-(crx_lon+2*crx_dist))
    idx_lon_e=np.squeeze(np.argwhere(tmp==np.min(tmp)))
    print(' Longitude limit to',idx_lon_e)

    # latitude has descending order: high -> low
    tmp=np.abs(lat[:,0]-(crx_lat+crx_dist))
    idx_lat_s=np.squeeze(np.argwhere(tmp==np.min(tmp)))
    print(' Latitude limit from',idx_lat_s)
    tmp=np.abs(lat[:,0]-(crx_lat-crx_dist))
    idx_lat_e=np.squeeze(np.argwhere(tmp==np.min(tmp)))
    print(' Latitude limit to',idx_lat_e)


    lon_box=lon[idx_lat_s:idx_lat_e,idx_lon_s:idx_lon_e]
    lat_box=lat[idx_lat_s:idx_lat_e,idx_lon_s:idx_lon_e]
    s3d_box=sval_3d[idx_lat_s:idx_lat_e,idx_lon_s:idx_lon_e,:]

    print(lon_box.shape)
    print(lat_box.shape)
    print(s3d_box.shape)
    
    nm_svar='Reflectivity (dBZ)'
    cs_cmap,cs_norm=new_cmap_3d()
    lb_ext='max'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=2
    cmap_range='designed'

    # Highest and lowest longitudes and latitudes for plot extent
#    if domain=='HRRR':
#        lon_min=-128.39
#        lon_max=-66.62
#        lat_min=25.12
#        lat_max=49.23
#    else:
#        lon_min=np.min(lon)
#        lon_max=np.max(lon)
#        lat_min=np.min(lat)
#        lat_max=np.max(lat)
    lon_min=crx_lon-2*crx_dist
    lon_max=crx_lon+2*crx_dist
    lat_min=crx_lat-crx_dist
    lat_max=crx_lat+crx_dist



    print(' lon_min=',lon_min,', lon_max=',lon_max)
    print(' lat_min=',lat_min,', lat_max=',lat_max)

    esp=1    
    extent=[lon_min-esp,lon_max+esp,lat_min-esp,lat_max+esp]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])
 
    # Max and Min of the field
    fmax=np.nanmax(sval)
    fmin=np.nanmin(sval)

    print(' fld_min=',fmin)
    print(' fld_max=',fmax)

    # Make the colormap range symmetry
    print(' cmap range=',cmap_range)
    if cmap_range=='symmetry': 
        tmp_cmp=max(abs(fmax),abs(fmin))
        cs_min=round(-tmp_cmp,n_rnd)
        cs_max=round(tmp_cmp,n_rnd)
    elif cmap_range=='round':
        cs_min=round(fmin,n_rnd)
        cs_max=round(fmax,n_rnd)
    elif cmap_range=='real':
        cs_min=fmin
        cs_max=fmax
    elif cmap_range=='fixed':
        cs_min=-10.0
        cs_max=10.0
    elif cmap_range=='designed':
        cs_min=0
        cs_max=70
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_min=',cs_min)
    print(' cs_max=',cs_max)
 
    
    fig=plt.figure()
    ax3d=fig.add_axes([0,0,1,1],projection='3d')
    # adjust axis aspect ratio
    ax3d.get_proj=lambda: np.dot(Axes3D.get_proj(ax3d), np.diag([1,1,1.5,1]))

    proj_ax=plt.figure().add_axes([0,0,1,1],projection=ccrs.PlateCarree())

    i_tmp=1
    plt_hmax=len(plt_hght_i)
    for il in plt_hght_i: 
        print(i_tmp,'/',plt_hmax)
        cs=proj_ax.contourf(lon_box,lat_box,s3d_box[:,:,il],np.linspace(0,70,15),cmap=cs_cmap,norm=cs_norm,
                    transform=ccrs.PlateCarree(),vmin=cs_min,vmax=cs_max,alpha=0.9)
        ax3d.projection=proj_ax.projection
        add_contourf3d(ax3d,cs,hght_km[il])
        i_tmp=i_tmp+1

    back_plot_3d(proj_ax,ax3d)
                       
    ax3d.view_init(20,-90)    #(elev,angle)

    ax3d.xaxis.pane.fill=False
    ax3d.yaxis.pane.fill=False
    ax3d.zaxis.pane.fill=False
    ax3d.grid(False)

    lb_font_sz=7
    ax3d.set_xlabel('Longitude',fontsize=lb_font_sz)
    ax3d.set_ylabel('Latitude',fontsize=lb_font_sz)
    ax3d.set_zlabel('Height (km)',fontsize=lb_font_sz) #,rotation=180)
    ax3d.xaxis.set_tick_params(labelsize=lb_font_sz)
    ax3d.yaxis.set_tick_params(labelsize=lb_font_sz)
    ax3d.zaxis.set_tick_params(labelsize=lb_font_sz)

    plt.close(proj_ax.figure)


    # Output figure
    ndpi=300
    out_file(out_fname_3d,ndpi)



# new colormap option ====================================== CHJ =====
def new_cmap():
# ========================================================== CHJ =====
    c_lvls=np.linspace(5,70,14)
    c_list=['turquoise','dodgerblue','mediumblue','lime','limegreen','green','#EEEE00','#EEC900','darkorange','red','firebrick','darkred','fuchsia']
    new_cmap=colors.ListedColormap(c_list)
    new_norm=colors.BoundaryNorm(c_lvls,new_cmap.N)

    return new_cmap,new_norm
       
  

# new colormap option for 3d =============================== CHJ =====
def new_cmap_3d():
# ========================================================== CHJ =====
    c_lvls=np.linspace(0,70,15)
    c_list=['ghostwhite','turquoise','dodgerblue','mediumblue','lime','limegreen','green','#EEEE00','#EEC900','darkorange','red','firebrick','darkred','fuchsia']
    new_cmap=colors.ListedColormap(c_list)
    new_norm=colors.BoundaryNorm(c_lvls,new_cmap.N)

    return new_cmap,new_norm

 

# Background plot ========================================== CHJ =====
def back_plot(ax):
# ========================================================== CHJ =====
    fline_wd=0.5  # line width
    falpha=0.3 # transparency

    # natural_earth
    land=cfeature.NaturalEarthFeature('physical','land',back_res,
                      edgecolor='face',facecolor=cfeature.COLORS['land'],
                      alpha=falpha)
    lakes=cfeature.NaturalEarthFeature('physical','lakes',back_res,
                      edgecolor='blue',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)
    coastline=cfeature.NaturalEarthFeature('physical','coastline',
                      back_res,edgecolor='blue',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)
    states=cfeature.NaturalEarthFeature('cultural','admin_1_states_provinces',
                      back_res,edgecolor='black',facecolor='none',
                      linewidth=fline_wd,linestyle=':',alpha=falpha)
    borders=cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                      back_res,edgecolor='red',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)

    # high-resoultion background image
    if back_img=='on':
        ax.background_img(name='NE', resolution='high')
 
    # On/off features
#    ax.add_feature(land)
    ax.add_feature(lakes)
    ax.add_feature(states)
    ax.add_feature(borders)
    ax.add_feature(coastline)



# Background plot for 3D =================================== CHJ =====
def back_plot_3d(proj_ax,ax3d):
# ========================================================== CHJ =====
    fline_wd=0.5  # line width
    falpha=0.3 # transparency

    # natural_earth
#    land=cfeature.NaturalEarthFeature('physical','land',back_res)
    lakes=cfeature.NaturalEarthFeature('physical','lakes',back_res)
    coastline=cfeature.NaturalEarthFeature('physical','coastline',back_res)
    states=cfeature.NaturalEarthFeature('cultural','admin_1_states_provinces',back_res)
    borders=cfeature.NaturalEarthFeature('cultural','admin_0_countries',back_res)


    clip_geom=proj_ax._get_extent_geom().buffer(0)
    
    add_feature3d(ax3d,lakes,clip_geom,0,'blue','-',fline_wd,falpha)
    add_feature3d(ax3d,states,clip_geom,0,'black',':',fline_wd,falpha)
    add_feature3d(ax3d,borders,clip_geom,0,'red','-',fline_wd,falpha)
    add_feature3d(ax3d,coastline,clip_geom,0,'blue','-',fline_wd,falpha)



# add_contourf3d ========================================== CHJ ====
def add_contourf3d(ax3d,contour_set,zlev):
# ========================================================= CHJ ===
    proj_ax = contour_set.collections[0].axes
    for collection in contour_set.collections:
        paths = collection.get_paths()
        trans_to_proj = collection.get_transform() - proj_ax.transData
        paths = [trans_to_proj.transform_path(path) for path in paths]
        verts = [path.vertices for path in paths]
        codes = [path.codes for path in paths]
        pc = PolyCollection([])
        pc.set_verts_and_codes(verts, codes)

        pc.set_facecolor(collection.get_facecolor())
        pc.set_edgecolor(collection.get_edgecolor())
        pc.set_alpha(collection.get_alpha())

        ax3d.add_collection3d(pc,zs=zlev)

    proj_ax.autoscale_view()

    ax3d.set_xlim(*proj_ax.get_xlim())
    ax3d.set_ylim(*proj_ax.get_ylim())
    ax3d.set_zlim((0,max_zaxis))



# add_feature3d =========================================== CHJ =====
def add_feature3d(ax3d,feature,clip_geom,zbase,lcol,lsty,lwid,la):
# ========================================================= CHJ =====
    concat = lambda iterable: list(itertools.chain.from_iterable(iterable))

    target_projection = ax3d.projection
    geoms = list(feature.geometries())

    if target_projection != feature.crs:
        geoms = [target_projection.project_geometry(geom, feature.crs)
                 for geom in geoms]

    if clip_geom:
        geoms = [geom.intersection(clip_geom) for geom in geoms]

    paths = concat(geos_to_path(geom) for geom in geoms)

    kwargs = feature.kwargs
    if kwargs.get('edgecolor') == 'face':
        kwargs['edgecolor'] = kwargs['facecolor']

    polys = concat(path.to_polygons(closed_only=False) for path in paths)

    lc=LineCollection(polys,colors=lcol,linestyle=lsty,linewidths=lwid,alpha=la)
#    if kwargs.get('facecolor', 'none') == 'none':
#        lc = LineCollection(polys, **kwargs)
#    else:
#        lc = PolyCollection(polys, closed=False, **kwargs)

    ax3d.add_collection3d(lc, zs=zbase)




# Output file ============================================= CHJ =====
def out_file(out_file,ndpi):
# ========================================================= CHJ =====
    # Output figure
    plt.savefig(out_fig_dir+out_file+'.png',dpi=ndpi,bbox_inches='tight')
    plt.close('all')



# Main call ================================================ CHJ =====
if __name__=='__main__':
    main()

