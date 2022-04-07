###################################################################### CHJ #####
## Name		: plot_fv3lam_multigrid.py
## Language	: Python 3.7
## Usage	: Plot boundaries of multiple grids on the map
## Input files  : grid.tile7.haloX.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2022/04/06: Chan-Hoo Jeon : Preliminary version
###################################################################### CHJ #####

import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

# HPC machine ('hera','orion')
machine='hera'

print(' You are on', machine)

#### Machine-specific input data ==================================== CHJ =====
# cartopy.config: Natural Earth data for background
# out_fig_dir: directory where the output files are created
# mfdt_kwargs: mfdataset argument

if machine=='hera':
    cartopy.config['data_dir']='/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/NaturalEarth'
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig/"
    mfdt_kwargs={'parallel':False}
elif machine=='orion':
    cartopy.config['data_dir']='/home/chjeon/tools/NaturalEarth'
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
    mfdt_kwargs={'parallel':False,'combine':'by_coords'}
else:
    sys.exit('ERROR: Required input data are NOT set !!!')

plt.switch_backend('agg')

# Case-dependent input =============================================== CHJ =====
# Path to the directory where the grid file is located.
dir_grid1="/scratch2/NCEPDEV/stmp3/Chan-hoo.Jeon/expt_dirs/uwm_aqm_hrrr25/grid/"
dir_grid2="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/ufs_srw_app/srw_dev_test/expt_dirs/grid_RRFS_CONUS_25km_ics_FV3GFS_lbcs_FV3GFS_suite_GFS_v16/grid/"

# Grid file
fnm_grid1="C403_grid.tile7.halo4.nc"
fnm_grid2="C403_grid.tile7.halo4.nc"

# Domain name
domain_nm1='GSD_HRRR'
domain_nm2='RRFS_CONUS'

# output title and file names
out_grd_title='FV3LAM::multi-grids'
out_grd_fname='fv3lam_multigrid_'+domain_nm1+'_vs_'+domain_nm2

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'


# Main part (will be called at the end) ============================== CHJ =====
def main():
# ==================================================================== CHJ =====

    grid_info(dir_grid1,dir_grid2,fnm_grid1,fnm_grid2)
    grid_bndr_plot()


# Grid boundary plot ================================================= CHJ =====
def grid_info(dir_grid1,dir_grid2,fnm_grid1,fnm_grid2):
# ==================================================================== CHJ =====

    global grd_x1,grd_y1,grd_x2,grd_y2
    global lon_min,lon_max,lat_min,lat_max
    global npx1,npy1,npx2,npy2

    print(' ===== GRID 1 ====================================')
    # open grid file
    fname=os.path.join(dir_grid1,fnm_grid1)
    try: grd1=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(grd1)

    # Extract longitudes, latitudes, and others
    grd_x1=np.ma.masked_invalid(grd1['x'].data)
    grd_y1=np.ma.masked_invalid(grd1['y'].data)

    # array size
    (nyp1,nxp1)=grd_x1.shape
    print('super-grid array size (nyp1,nxp1)=',grd_x1.shape)
    npx1=int(nxp1/2)
    npy1=int(nyp1/2)
    print('orography array size (npy1,npx1)=',npy1,npx1)

    print(' ===== GRID 2 ====================================')
    # open grid file
    fname=os.path.join(dir_grid2,fnm_grid2)
    try: grd2=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(grd2)

    # Extract longitudes, latitudes, and others
    grd_x2=np.ma.masked_invalid(grd2['x'].data)
    grd_y2=np.ma.masked_invalid(grd2['y'].data)

    # array size
    (nyp2,nxp2)=grd_x2.shape
    print('super-grid array size (nyp2,nxp2)=',grd_x2.shape)
    npx2=int(nxp2/2)
    npy2=int(nyp2/2)
    print('orography array size (npy2,npx2)=',npy2,npx2)


    # Hightest/Lowest longitudes and latitudes for text
    lon1_min=np.min(grd_x1)
    lon1_max=np.max(grd_x1)
    lat1_min=np.min(grd_y1)
    lat1_max=np.min(grd_y1)
   
    lon2_min=np.min(grd_x2)
    lon2_max=np.max(grd_x2)
    lat2_min=np.min(grd_y2)
    lat2_max=np.min(grd_y2)

    lon_min=min(lon1_min,lon2_min)
    lon_max=max(lon1_max,lon2_max)
    lat_min=min(lat1_min,lat2_min)
    lat_max=min(lat1_max,lat2_max)




# Grid boundary plot ================================================= CHJ =====
def grid_bndr_plot():
# ==================================================================== CHJ =====


    print(' ===== boundary ===== GRID 1 ====================================')
    # Boundary: 1C (1st column of the array) 
    grd1_B1C_lon=grd_x1[:,0]
    grd1_B1C_lat=grd_y1[:,0]
    # Boundary: 1R (1st row of the array)
    grd1_B1R_lon=grd_x1[0,:]
    grd1_B1R_lat=grd_y1[0,:]
    # Boundary: xC (last column of the array)
    grd1_BxC_lon=grd_x1[:,-1]
    grd1_BxC_lat=grd_y1[:,-1]
    # Boundary: xR (last row of the array)
    grd1_BxR_lon=grd_x1[-1,:]
    grd1_BxR_lat=grd_y1[-1,:]

    print(' ===== boundary ===== GRID 1 ====================================')
    # Boundary: 1C (1st column of the array) 
    grd2_B1C_lon=grd_x2[:,0]
    grd2_B1C_lat=grd_y2[:,0]
    # Boundary: 1R (1st row of the array)
    grd2_B1R_lon=grd_x2[0,:]
    grd2_B1R_lat=grd_y2[0,:]
    # Boundary: xC (last column of the array)
    grd2_BxC_lon=grd_x2[:,-1]
    grd2_BxC_lat=grd_y2[:,-1]
    # Boundary: xR (last row of the array)
    grd2_BxR_lon=grd_x2[-1,:]
    grd2_BxR_lat=grd_y2[-1,:]
  
    print(' ===== plot boundaries ====================================')

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(
                        central_longitude=-107,central_latitude=53)))

    ax.set_title(out_grd_title, fontsize=9)

    # Call background plot
    back_plot(ax)

    # Scatter plot (zorder: lowest-plot on bottom, highest-plot on top)
    sp_scale=0.1
    ax.scatter(grd1_B1C_lon,grd1_B1C_lat,transform=ccrs.PlateCarree(),c='red',s=sp_scale,zorder=3)
    ax.scatter(grd1_B1R_lon,grd1_B1R_lat,transform=ccrs.PlateCarree(),c='red',s=sp_scale,zorder=3)
    ax.scatter(grd1_BxC_lon,grd1_BxC_lat,transform=ccrs.PlateCarree(),c='red',s=sp_scale,zorder=3)
    ax.scatter(grd1_BxR_lon,grd1_BxR_lat,transform=ccrs.PlateCarree(),c='red',s=sp_scale,zorder=3)

    ax.scatter(grd2_B1C_lon,grd2_B1C_lat,transform=ccrs.PlateCarree(),c='blue',s=sp_scale,zorder=4)
    ax.scatter(grd2_B1R_lon,grd2_B1R_lat,transform=ccrs.PlateCarree(),c='blue',s=sp_scale,zorder=4)
    ax.scatter(grd2_BxC_lon,grd2_BxC_lat,transform=ccrs.PlateCarree(),c='blue',s=sp_scale,zorder=4)
    ax.scatter(grd2_BxR_lon,grd2_BxR_lat,transform=ccrs.PlateCarree(),c='blue',s=sp_scale,zorder=4)

    # Add text to each boundary
    tsize=7
    txt_x=grd1_B1R_lon[npx1]-0.5
    txt_y=grd1_B1R_lat[npy1]-0.5
    ax.text(txt_x,txt_y,domain_nm1,color='red',fontsize=tsize,transform=ccrs.PlateCarree())
    txt_x=grd2_BxR_lon[npx2]+2
    txt_y=grd2_BxR_lat[npy2]+2
    ax.text(txt_x,txt_y,domain_nm2,color='blue',fontsize=tsize,transform=ccrs.PlateCarree())

    # Output figure
    ndpi=300
    out_file(out_grd_fname,ndpi)




# Background plot ========================================== CHJ =====
def back_plot(ax):
# ========================================================== CHJ =====
    fline_wd=0.5  # line width
    falpha=0.3 # transparency

    # natural_earth
#    land=cfeature.NaturalEarthFeature('physical','land',back_res,
#                      edgecolor='face',facecolor=cfeature.COLORS['land'],
#          	      alpha=falpha)
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

#    ax.add_feature(land)
    ax.add_feature(lakes)
    ax.add_feature(states)
    ax.add_feature(borders)
    ax.add_feature(coastline)
   


# Output file ============================================= CHJ =====
def out_file(out_file,ndpi):
# ========================================================= CHJ =====
    # Output figure

    plt.savefig(out_fig_dir+out_file+'.png',dpi=ndpi,bbox_inches='tight')
    plt.close('all')



# Main call ================================================ CHJ =====
if __name__=='__main__':
    main()

