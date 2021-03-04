###################################################################### CHJ #####
## Name		: plot_fv3lam_gridonly.py
## Language	: Python 3.7
## Usage	: Plot regional FV3 super-grid on the map
## Input files  : grid.tile7.haloX.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/07/20: Chan-Hoo Jeon : Preliminary version
## V001: 2021/03/04: Chan-Hoo Jeon : Simplify the script
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

# Path to Natural Earth Data-set for background plot
if machine=='hera':
    cartopy.config['data_dir']='/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/NaturalEarth'
elif machine=='orion':
    cartopy.config['data_dir']='/home/chjeon/tools/NaturalEarth'
else:
    sys.exit('ERROR: path to Natural Earth Data is not set !!!')

plt.switch_backend('agg')

# Global variables ======================================== CHJ =====
# ..... Case-dependent input :: should be changed case-by-case .....
# ******
# INPUT
# ******
# Path to the directory where the grid file is located.
dnm_data="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/expt_dirs/test_community/2020122700/INPUT/"

# Grid file
fnm_in_grid='grid.tile7.halo4.nc'

# *******
# OUTPUT
# *******
# Path to directory
if machine=='hera':
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig/"
elif machine=='orion':
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
else:
    sys.exit('ERROR: path to output directory is not set !!!')

# output title and file names
out_grd_title='FV3LAM::grid'
out_grd_fname='fv3lam_grid_only'

# Colormap range option flag ('symmetry','roudn','real','fixed')
cmap_range_grd='round'

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'

# Machine-specific mfdataset arguments
if machine=='hera':
    mfdt_kwargs={'parallel':False}
elif machine=='orion':
    mfdt_kwargs={'parallel':False,'combine':'by_coords'}
else:
    mfdt_kwargs={'parallel':False}


# Main part (will be called at the end) ============================== CHJ =====
def main():
# ==================================================================== CHJ =====
    global grd_x,grd_y,npx,npy
     # open the grid file
    fname=os.path.join(dnm_data,fnm_in_grid)
    try: grd=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)

    print(' ===== GRID ====================================')
    print(grd)

    # Extract longitudes, latitudes, and others
    grd_x=np.ma.masked_invalid(grd['x'].data)
    grd_y=np.ma.masked_invalid(grd['y'].data)
    grd_area=np.ma.masked_invalid(grd['area'].data)

    # array size
    (nyp,nxp)=grd_x.shape
    print('grid array size=',grd_x.shape)
    (ny,nx)=grd_area.shape
    print('area array size=',grd_area.shape)

    npx=int(nxp/2)
    npy=int(nyp/2)
    print('orography array size (npy,npx)=',npy,npx)

    # Hightest/Lowest longitudes and latitudes for text
    lon_min=np.min(grd_x)
    lon_max=np.max(grd_x)
    lat_min=np.min(grd_y)
    lat_max=np.min(grd_y)

    # Plot grid sizes
    grid_dxy_plot(grd_area,lon_min,lat_min) 

    # Plot boundary
    grid_bndr_plot(nxp,nyp)
   


# Grid plot: dx/dy =================================================== CHJ =====
def grid_dxy_plot(grd_area,lon_min,lat_min):
# ==================================================================== CHJ =====
    global c_lon,c_lat,extent
    print(' ===== cell size ===== GRID =========================================')

    oro_x=np.zeros((npy,npx))
    oro_y=np.zeros((npy,npx))
    cell_area=np.zeros((npy,npx))
    cell_dxy=np.zeros((npy,npx))

    for iy in range(npy):
        for jx in range(npx):
            iy2=2*iy
            iy2p1=2*iy+1
            jx2=2*jx
            jx2p1=2*jx+1
            
            oro_x[iy,jx]=grd_x[iy2p1,jx2p1]
            oro_y[iy,jx]=grd_y[iy2p1,jx2p1]
            cell_area[iy,jx]=grd_area[iy2,jx2]+grd_area[iy2,jx2p1]+grd_area[iy2p1,jx2]+grd_area[iy2p1,jx2p1]
            cell_dxy[iy,jx]=math.sqrt(cell_area[iy,jx])/1000

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(oro_x)
    lon_max=np.max(oro_x)
    lat_min=np.min(oro_y)
    lat_max=np.max(oro_y)

    print(' ***** Ref. for lon1/lat1/lon2/lat2 in model_configure *****')
    print(' oro:lon-min(lon1)=',lon_min-360)
    print(' oro:lon-max(lon2)=',lon_max-360)
    print(' oro:lat-min(lat1)=',lat_min)
    print(' oro:lat-max(lat2)=',lat_max)

    print(' ***** npx/npy in input.nml/fv_core_nml *****')
    hcond=fnm_in_grid[-8:-3]
    if hcond=='halo0':
        print(' npx=',npx+1)
        print(' npy=',npy+1)
    elif hcond=='halo3':
        print(' npx=',npx-5)
        print(' npy=',npy-5)     
    elif hcond=='halo4':
        print(' npx=',npx-7)
        print(' npy=',npy-7)
    else:
        sys.exit('ERROR: wrong fnm_in_base !!!!!')

    # Plot extent
    extent=[lon_min-5,lon_max+5,lat_min-5,lat_max+3]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])


    # Max and Min of the field
    fmax=np.max(cell_dxy)
    fmin=np.min(cell_dxy)
    favg=np.average(cell_dxy)
    print(' fld_max=',fmax)
    print(' flx_min=',fmin)
    print(' flx_avg=',favg)

    # Set the colormap range
#    cmap_range_grd='round'
    n_rnd=2
    print(' cmap range=',cmap_range_grd)
    if cmap_range_grd=='symmetry':
        tmp_cmp=max(abs(fmax),abs(fmin))
        cs_min=round(-tmp_cmp,n_rnd)
        cs_max=round(tmp_cmp,n_rnd)
        cs_avg=round(favg,n_rnd)
    elif cmap_range_grd=='round':
        cs_min=round(fmin,n_rnd)
        cs_max=round(fmax,n_rnd)
        cs_avg=round(favg,n_rnd)
    elif cmap_range_grd=='real':
        cs_min=fmin
        cs_max=fmax
        cs_avg=favg
    elif cmap_range_grd=='fixed':
        cs_min=2.7
        cs_max=3.2
        cs_avg=round(favg,n_rnd)
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_min=',cs_min)
    print(' cs_max=',cs_max)
    print(' cs_avg=',cs_avg)


    nm_svar='Cell size (km)'
    cs_cmap='nipy_spectral_r'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())

    out_grd_dx_title=out_grd_title+'::Cell Size'
    ax.set_title(out_grd_dx_title, fontsize=9)

    # Call background plot
    back_plot(ax)

    cs=ax.pcolormesh(oro_x,oro_y,cell_dxy,cmap=cs_cmap,rasterized=True,
              vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())

    ref_txt='Max='+str(round(fmax,2))+', Min='+str(round(fmin,2))+', Avg='+str(round(favg,2))
    plt.text(lon_min-2,lat_min-2,ref_txt,horizontalalignment='left',
             transform=ccrs.Geodetic(),fontsize=9)

    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(nm_svar,fontsize=8)


    # Output figure
    out_grd_dx_fname=out_grd_fname+'_dxy'
    ndpi=300
    out_file(out_grd_dx_fname,ndpi)



# Grid boundary plot ================================================= CHJ =====
def grid_bndr_plot(nxp,nyp):
# ==================================================================== CHJ =====
    print(' ===== boundary ===== GRID ====================================')

    # Boundary: 1C (1st column of the array) 
    grd_B1C_lon=grd_x[:,0]
    grd_B1C_lat=grd_y[:,0]
    # Boundary: 1R (1st row of the array)
    grd_B1R_lon=grd_x[0,:]
    grd_B1R_lat=grd_y[0,:]
    # Boundary: xC (last column of the array)
    grd_BxC_lon=grd_x[:,-1]
    grd_BxC_lat=grd_y[:,-1]
    # Boundary: xR (last row of the array)
    grd_BxR_lon=grd_x[-1,:]
    grd_BxR_lat=grd_y[-1,:]

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
 
    out_grd_bndr_title=out_grd_title+'::Boundary'
    ax.set_title(out_grd_bndr_title, fontsize=9)

    # Call background plot
    back_plot(ax)

    # Scatter plot (zorder: lowest-plot on bottom, highest-plot on top)
    sp_scale=0.1
    ax.scatter(grd_B1C_lon,grd_B1C_lat,transform=ccrs.PlateCarree(),c='red',s=sp_scale,label='B1C',zorder=3)
    ax.scatter(grd_B1R_lon,grd_B1R_lat,transform=ccrs.PlateCarree(),c='blue',s=sp_scale,label='B1R',zorder=3)
    ax.scatter(grd_BxC_lon,grd_BxC_lat,transform=ccrs.PlateCarree(),c='purple',s=sp_scale,label='BxC',zorder=3)
    ax.scatter(grd_BxR_lon,grd_BxR_lat,transform=ccrs.PlateCarree(),c='green',s=sp_scale,label='BxR',zorder=3)
    # Add text to each boundary
    tsize=9
    ntxt=int(nyp/2)
    txt_x=grd_B1C_lon[ntxt]+1
    txt_y=grd_B1C_lat[ntxt]-1
    ax.text(txt_x,txt_y,'B1C',color='red',fontsize=tsize,transform=ccrs.PlateCarree())
    txt_x=grd_BxC_lon[ntxt]+1
    txt_y=grd_BxC_lat[ntxt]+1
    ax.text(txt_x,txt_y,'BxC',color='purple',fontsize=tsize,transform=ccrs.PlateCarree())
    ntxt=int(nxp/2)
    txt_x=grd_B1R_lon[ntxt]
    txt_y=grd_B1R_lat[ntxt]+1
    ax.text(txt_x,txt_y,'B1R',color='blue',fontsize=tsize,transform=ccrs.PlateCarree())
    txt_x=grd_BxR_lon[ntxt]
    txt_y=grd_BxR_lat[ntxt]+1
    ax.text(txt_x,txt_y,'BxR',color='green',fontsize=tsize,transform=ccrs.PlateCarree())
 
    # File name
    out_grd_bndr_fname=out_grd_fname+'_bndr'
    # Output figure
    ndpi=300
    out_file(out_grd_bndr_fname,ndpi)




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

