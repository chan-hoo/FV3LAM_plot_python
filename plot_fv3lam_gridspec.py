###################################################################### CHJ #####
## Name		: plot_fv3lam_gridspec.py
## Language	: Python 3.7
## Usage	: Plot an output, grid_spec, for fv3 regional modeling
## Input files  : grid_spec.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/04/24: Chan-Hoo Jeon : Preliminary version
## V001: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
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
# Path to the directory where the input NetCDF file is located.
dnm_out="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/run_C96/"

# grid file name
fnm_input='grid_spec.nc'

# Domain name:
domain='CONUS'

# Grid resolution ('C96'/'C768'):
res='C96'

# Grid type ('ESG'/'GFDL')
gtype='GFDL'

# GFDL grid-refinement ratio (for ESG grid, refine=0)
if gtype=='ESG':
    refine=0
elif gtype=='GFDL':
    refine=3

# Number of grid points in plotting
n_gpt=20

# Input grid file name (grid/orography)
dnm_ref="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/C96/"
fnm_ref_grd=res+'_grid.tile7.halo0.nc'
fnm_ref_oro=res+'_oro_data.tile7.halo0.nc'


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

# basic forms of title and file name
if gtype=='ESG':
    out_grd_title='GRID_SPEC::'+domain+'(ESG)::'+res
    out_grd_fname='fv3_out_grdspec_'+domain+'_esg_'+res
elif gtype=='GFDL':
    out_grd_title='GRID_SPEC::'+domain+'(GFDL)::'+res+'(x'+str(refine)+')'
    out_grd_fname='fv3_out_grdspec_'+domain+'_gfdl_'+res

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'

# Machine-specific mfdataset arguments
if machine=='hera':
    mfdt_kwargs={'parallel':False}
elif machine=='orion':
    mfdt_kwargs={'parallel':False,'combine':'by_coords'}
else:
    mfdt_kwargs={'parallel':False}


# Main part (will be called at the end) =================== CHJ =====
def main():
# ========================================================= CHJ =====

    print(' ===== OUTPUT: grid_spec ================================')
    # open the data file
    fname=os.path.join(dnm_out,fnm_input)
    try: grdo=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(grdo)
    # Extract longitudes, and latitudes
    golon=np.ma.masked_invalid(grdo['grid_lon'].data)
    golat=np.ma.masked_invalid(grdo['grid_lat'].data)
    golont=np.ma.masked_invalid(grdo['grid_lont'].data)
    golatt=np.ma.masked_invalid(grdo['grid_latt'].data)

    print(' ===== REFERENCE: Super-grid (halo0) =======================')
    # open the data file
    fname=os.path.join(dnm_ref,fnm_ref_grd)
    try: refg=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(refg)
    # Extract grid info.
    ref_glon=np.ma.masked_invalid(refg['x'].data)
    ref_glat=np.ma.masked_invalid(refg['y'].data)
  
    print(' ===== REFERENCE: Orography-grid (halo0) ==================')
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
    ref_glon_c=ref_glon[:n_gpt2,:n_gpt2]
    ref_glat_c=ref_glat[:n_gpt2,:n_gpt2]
    
    ref_olon_c=ref_olon[:n_gpt,:n_gpt]
    ref_olat_c=ref_olat[:n_gpt,:n_gpt]

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(ref_glon_c)
    lon_max=np.max(ref_glon_c)
    lat_min=np.min(ref_glat_c)
    lat_max=np.max(ref_glat_c)

    # Plot extent
    extent=[lon_min-0.1,lon_max+0.1,lat_min-0.1,lat_max+0.1]
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

    ref_txt='First '+str(n_gpt)+' grid points from Top-Right'
    plt.text(lon_min,lat_min-0.05,ref_txt,transform=ccrs.Geodetic(),fontsize=8)
    plt.legend((s1,s2,s3,s4),('super-grid','oro-grid','lon/lat','lont/latt'),scatterpoints=1,loc='upper right',ncol=2,fontsize=8)

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

