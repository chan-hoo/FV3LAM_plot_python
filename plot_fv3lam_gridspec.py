###################################################################### CHJ #####
## Name		: plot_fv3lam_gridspec.py
## Language	: Python 3.7
## Usage	: Plot an output, grid_spec, for fv3 regional modeling
## Input files  : grid_spec.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/04/24: Chan-Hoo Jeon : Preliminary version
## V001: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
## V002: 2021/03/05: Chan-Hoo Jeon : Simplify the script
## V003: 2021/07/06: Chan-Hoo Jeon : Add a plot for the entire domain
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
# Path to the directory where the input NetCDF file is located.
dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/ufs_srw_app/srw_dev_test/expt_dirs/inline_post/2019070100/"

# grid file name
fnm_input='grid_spec.nc'

# Domain name
domain_nm='RRFS_CONUS_25km'

# Grid point plot (every 'n_skip' rows/columns)
n_skip=5

# Flag for the reference grid (on/off)
i_ref='on'

# Reference grid file name (grid/orography w/ halo4)
dnm_ref=dnm_data+"INPUT/"
fnm_ref_grd='grid.tile7.halo4.nc'
fnm_ref_oro='oro_data.tile7.halo4.nc'

# basic forms of title and file name
out_grd_title='FV3LAM::GRID_SPEC::'+domain_nm
out_grd_fname='fv3lam_out_grdspec_'+domain_nm

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'


# Main part (will be called at the end) =================== CHJ =====
def main():
# ========================================================= CHJ =====
    
    global golon,golat,golont,golatt
    print(' ===== grid_spec ========================================')
    # open the data file
    fname=os.path.join(dnm_data,fnm_input)
    try: grdo=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(grdo)
    # Extract longitudes, and latitudes
    golon=np.ma.masked_invalid(grdo['grid_lon'].data)
    golat=np.ma.masked_invalid(grdo['grid_lat'].data)
    golont=np.ma.masked_invalid(grdo['grid_lont'].data)
    golatt=np.ma.masked_invalid(grdo['grid_latt'].data)

    grd_plot()

    if i_ref=='on':
        # Number of grid points in plotting
        n_gpt=20

        grd_ref(n_gpt)

# Grid plot =============================================== CHJ =====
def grd_plot():
# ========================================================= CHJ =====

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(golon)
    lon_max=np.max(golon)
    lat_min=np.min(golat)
    lat_max=np.max(golat)

    # Plot extent
    extent=[lon_min-0.1,lon_max+1,lat_min-0.1,lat_max+1]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])

    # corner grid points (every 'n_skip' rows/columns)
    grdx_slc=golon[::n_skip,::n_skip]
    grdy_slc=golat[::n_skip,::n_skip]

    # center grid points (every 'n_skip' rows/columns)
    grdtx_slc=golont[::n_skip,::n_skip]
    grdty_slc=golatt[::n_skip,::n_skip]

    if domain_nm[:7]=='RRFS_NA':
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(
                            central_longitude=-107,central_latitude=53)))
        ref_lon=-133.5
        ref_lat=lat_min-5.5
        lgd_loc='lower left'
    else:
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
        ax.set_extent(extent, ccrs.PlateCarree())
        ref_lon=lon_min-3
        ref_lat=lat_max
        lgd_loc='lower right'

    ax.set_title(out_grd_title, fontsize=8)

    # Call background plot
    back_plot(ax)

    # Scatter plot (zorder: lowest-plot on bottom, highest-plot on top)
    sp_scale=2
    s1=ax.scatter(grdx_slc,grdy_slc,transform=ccrs.PlateCarree(),marker='o',
         facecolors="none",edgecolors='red',linewidth=0.3,s=sp_scale,zorder=3)

    s2=ax.scatter(grdtx_slc,grdty_slc,transform=ccrs.PlateCarree(),marker='*',
         facecolors='blue',edgecolors='blue',linewidth=0.3,s=sp_scale,zorder=4)

    ref_txt='Every '+str(n_skip)+' (i,j)s'
    plt.text(ref_lon,ref_lat,ref_txt,transform=ccrs.Geodetic(),fontsize=6)
    plt.legend((s1,s2),('g-corner','g-center'),scatterpoints=1,loc=lgd_loc,ncol=1,fontsize=6)

    # Output figure
    ndpi=300
    out_file(out_grd_fname,ndpi)



# Grid comparison with reference grids ==================== CHJ =====
def grd_ref(n_gpt):
# ========================================================= CHJ =====

    print(' ===== REFERENCE: Super-grid (halo4) =======================')
    # open the data file
    fname=os.path.join(dnm_ref,fnm_ref_grd)
    try: refg=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(refg)
    # Extract grid info.
    ref_glon=np.ma.masked_invalid(refg['x'].data)
    ref_glat=np.ma.masked_invalid(refg['y'].data)
  
    print(' ===== REFERENCE: Orography-grid (halo4) ==================')
    # open the data file
    fname=os.path.join(dnm_ref,fnm_ref_oro)
    try: refo=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(refo)
    # Extract grid info.
    ref_olon=np.ma.masked_invalid(refo['geolon'].data)
    ref_olat=np.ma.masked_invalid(refo['geolat'].data)

    # Grid points (first 'n_gpt' rows/columns)
    golon_c=golon[:n_gpt,:n_gpt]
    golat_c=golat[:n_gpt,:n_gpt]
    golont_c=golont[:n_gpt,:n_gpt]
    golatt_c=golatt[:n_gpt,:n_gpt]

    n_gpt2=n_gpt*2
    nhalo=4
    nhalo2=nhalo*2
    ref_glon_c=ref_glon[:n_gpt2+nhalo2,:n_gpt2+nhalo2]
    ref_glat_c=ref_glat[:n_gpt2+nhalo2,:n_gpt2+nhalo2]
    
    ref_olon_c=ref_olon[:n_gpt+nhalo,:n_gpt+nhalo]
    ref_olat_c=ref_olat[:n_gpt+nhalo,:n_gpt+nhalo]

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(ref_glon_c)
    lon_max=np.max(ref_glon_c)
    lat_min=np.min(ref_glat_c)
    lat_max=np.max(ref_glat_c)

    # Plot extent
    extent=[lon_min-0.1,lon_max+1,lat_min-0.1,lat_max+1]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.set_title(out_grd_title, fontsize=9)

    # Call background plot
    back_plot(ax)

    # Scatter plot (zorder: lowest-plot on bottom, highest-plot on top)
    sp_scale=3
    # super-grid
    s1=ax.scatter(ref_glon_c,ref_glat_c,transform=ccrs.PlateCarree(),marker='o',
       facecolors="None",edgecolors='black',linewidth=0.3,s=sp_scale,zorder=3)
    # orography grid
    s2=ax.scatter(ref_olon_c,ref_olat_c,transform=ccrs.PlateCarree(),marker='s',
       facecolors="None",edgecolors='magenta',linewidth=0.3,s=sp_scale+1,zorder=4)
     # output grid
    s3=ax.scatter(golon_c,golat_c,transform=ccrs.PlateCarree(),marker='o',
       facecolors="red",edgecolors='red',linewidth=0.3,s=sp_scale,zorder=5)
    s4=ax.scatter(golont_c,golatt_c,transform=ccrs.PlateCarree(),marker='*',
       facecolors='blue',edgecolors='blue',linewidth=0.3,s=sp_scale,zorder=6)

    ref_txt='First '+str(n_gpt)+' grid points from Bottom-Left'
    plt.text(lon_min,lat_min-0.05,ref_txt,transform=ccrs.Geodetic(),fontsize=8)
    plt.legend((s1,s2,s3,s4),('super-grid','oro-grid','lon/lat','lont/latt'),scatterpoints=1,loc='upper right',ncol=2,fontsize=8)

    # Output figure
    out_grd_ref_fname=out_grd_fname+'_comref'
    ndpi=300
    out_file(out_grd_ref_fname,ndpi)

   

# Background plot ========================================== CHJ =====
def back_plot(ax):
# ========================================================== CHJ =====
    fline_wd=0.5  # line width
    falpha=0.3 # transparency

    # natural_earth
#    land=cfeature.NaturalEarthFeature('physical','land',back_res,
#                      edgecolor='face',facecolor=cfeature.COLORS['land'],
#                      alpha=falpha)
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

