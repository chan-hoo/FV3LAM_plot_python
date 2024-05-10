###################################################################### CHJ #####
## Name		: plot_landda_sfc_tile.py
## Language	: Python 3.7
## Usage	: Plot output of land-DA workflow
## Input files  : grid_spec.nc
## NOAA/EPIC
## History ===============================
## V000: 2024/02/01: Chan-Hoo Jeon : Preliminary version
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
path_orog="/scratch2/NAGAPE/epic/UFS_Land-DA/inputs/forcing/era5/orog_files/"
path_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/land_DA_test/workdir/mem000/"

# input file name
fnm_base_orog='oro_C96.mx100.tile'
fnm_base_data=''
fnm_date='20191221'

num_tiles=6

# Domain name
domain_nm='C96_25km'

# Grid point plot (every 'n_skip' rows/columns)
n_skip=5

# Flag for the reference grid (on/off)
i_ref='on'

# basic forms of title and file name
out_title='Land-DA::SFC::'+domain_nm
out_fname='landda_out_sfc_'+domain_nm

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'


# Main part (will be called at the end) =================== CHJ =====
def main():
# ========================================================= CHJ =====

    print(' ===== orography ========================================')
    # open the data file

    glon_all=[]
    glat_all=[]
    for it in range(num_tiles):
        itp=it+1
        fnm_orog=fnm_base_orog+str(itp)+'.nc'
        fname=os.path.join(path_orog,fnm_orog)

        try: orog=xr.open_mfdataset(fname,**mfdt_kwargs)
        except: raise Exception('Could NOT find the file',fname)
        print(orog)
        # Extract longitudes, and latitudes
        geolon=np.ma.masked_invalid(orog['geolon'].data)
        geolat=np.ma.masked_invalid(orog['geolat'].data)
        glon_all.append(geolon[None,:])
        glat_all.append(geolat[None,:])

    glon=np.vstack(glon_all)
    glat=np.vstack(glat_all)

    print(glon.shape)
    print(glat.shape)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(glon)
    lon_max=np.max(glon)
    lat_min=np.min(glat)
    lat_max=np.max(glat)

    # center of map
    c_lon=-77.0369

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
#    ax.set_extent(extent, ccrs.PlateCarree())
    ax.set_title(out_title, fontsize=6)

    # Call background plot
    back_plot(ax)

    # Scatter plot (zorder: lowest-plot on bottom, highest-plot on top)
    s_scale=0.05
    s_color=["r", "b", "c", "g", "y", "m"]

    # orography grid
    for it in range(num_tiles):  
        ax.scatter(glon[it,:,:],glat[it,:,:],transform=ccrs.PlateCarree(),marker='o',facecolors=s_color[it],s=s_scale,zorder=1)

    # Output figure
    ndpi=300
    out_file(out_fname,ndpi)

   

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

