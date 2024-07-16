#!/usr/bin/env python3

###################################################################### CHJ #####
## Name		: plot_landda_sfcanal.py
## Language	: Python 3.7
## Usage	: Plot sfcanal of land-DA analysis task
## Input files  : sfcanal.nc
## NOAA/EPIC
## History ===============================
## V000: 2024/07/16: Chan-Hoo Jeon : Preliminary version
###################################################################### CHJ #####

import os, sys
import numpy as np
import netCDF4 as nc
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

# HPC machine ('hera','orion')
machine='hera'

print(' You are on', machine)

#### Machine-specific input data ===================================== CHJ =====
# cartopy.config: Natural Earth data for background
# out_fig_dir: directory where the output files are created

if machine=='hera':
    cartopy.config['data_dir']='/scratch2/NAGAPE/epic/UFS_Land-DA_Dev/inputs/NaturalEarth'
    path_orog="/scratch2/NAGAPE/epic/UFS_Land-DA_Dev/inputs/FV3_fix_tiled/C96"
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig"
elif machine=='orion' or machine=='hercules':
    cartopy.config['data_dir']='/work/noaa/epic/UFS_Land-DA_Dev/inputs/NaturalEarth'
    path_orog="/work/noaa/epic/UFS_Land-DA_Dev/inputs/FV3_fix_tiled/C96"
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig"
else:
    sys.exit('ERROR: Required input data are NOT set in this machine !!!')

# Case-dependent input ============================================== CHJ =====

# input data:
path_dir='/scratch2/NAGAPE/epic/Chan-hoo.Jeon/global-workflow_test/comroot/snow_gsi/gdas.20211220/18/analysis/atmos'
fn_input='gdas.t18z.sfcanl.nc'

svars_list=['snod']

# Vertical layer number
zlyr_num=2

# basic forms of title and file name
out_title_base='Land-DA::SFCANAL::'
out_fn_base='landda_out_sfcanal_'

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'
# center of map
c_lon=-77.0369


# Main part (will be called at the end) ============================= CHJ =====
def main():
# =================================================================== CHJ =====

    global data_in,glon,glat
     # open the orography file
    fname=os.path.join(path_dir,fn_input)
    try: data_in=xr.open_mfdataset(fname)
    except: raise Exception('Could NOT find the file',fname)

    print(' ===== INPUT DATA ================================')
    print(data_in)

    # Extract longitudes, and latitudes
    glon=np.ma.masked_invalid(data_in['lon'].data)
    glat=np.ma.masked_invalid(data_in['lat'].data)

    for svar in svars_list:
        plot_svar(svar)


# data plot ========================================================== CHJ =====
def plot_svar(svar):
# ==================================================================== CHJ =====

    # Extract variable
    svar_data=np.ma.masked_invalid(data_in[svar].data)
    print(svar_data.shape)
    data2d=np.squeeze(svar_data,axis=0)

    var_max=np.max(data2d)
    var_min=np.min(data2d)
    print('var_max=',var_max)
    print('var_min=',var_min)
    var_max08=var_max*0.8
    var_min08=var_min*0.8
    print('var_max08=',var_max08)
    print('var_min08=',var_min08)

    cmap_range_opt='round'
    cs_cmap='gist_ncar_r'
    if cmap_range_opt=='symmetry':
        n_rnd=0
        tmp_cmp=max(abs(var_max08),abs(var_min08))
        cs_min=round(-tmp_cmp,n_rnd)
        cs_max=round(tmp_cmp,n_rnd)
        cbar_extend='both'
    elif cmap_range_opt=='round':
        n_rnd=1
        cs_min=round(var_min08,n_rnd)
        cs_max=round(var_max08,n_rnd)
        cbar_extend='both'
    elif cmap_range_opt=='real':
        cs_min=var_min
        cs_max=var_max
        cbar_extend='neither'
    elif cmap_range_opt=='fixed':
        cs_min=0.0
        cs_max=150.0
        cbar_extend='both'
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print('cs_max=',cs_max)
    print('cs_min=',cs_min)

    out_title=out_title_base+svar
    out_fn=out_fn_base+svar

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_title(out_title, fontsize=6)
    # Call background plot
    back_plot(ax)

    cs=ax.pcolormesh(glon,glat,data2d,cmap=cs_cmap,rasterized=True,
        vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())

    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=cbar_extend)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(svar,fontsize=6)
    # Output figure
    ndpi=300
    out_file(out_fn,ndpi)


# Background plot ==================================================== CHJ =====
def back_plot(ax):
# ==================================================================== CHJ =====
    fline_wd=0.5  # line width
    falpha=0.7 # transparency

    # natural_earth
    land=cfeature.NaturalEarthFeature('physical','land',back_res,
                      edgecolor='face',facecolor=cfeature.COLORS['land'],
                      alpha=falpha)
    lakes=cfeature.NaturalEarthFeature('physical','lakes',back_res,
                      edgecolor='blue',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)
    coastline=cfeature.NaturalEarthFeature('physical','coastline',
                      back_res,edgecolor='black',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)
    states=cfeature.NaturalEarthFeature('cultural','admin_1_states_provinces',
                      back_res,edgecolor='green',facecolor='none',
                      linewidth=fline_wd,linestyle=':',alpha=falpha)
    borders=cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                      back_res,edgecolor='red',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)

#    ax.add_feature(land)
#    ax.add_feature(lakes)
#    ax.add_feature(states)
#    ax.add_feature(borders)
    ax.add_feature(coastline)


# Output file ======================================================= CHJ =====
def out_file(out_file,ndpi):
# =================================================================== CHJ =====
    # Output figure
    fp_out=os.path.join(out_fig_dir,out_file)
    plt.savefig(fp_out+'.png',dpi=ndpi,bbox_inches='tight')
    plt.close('all')


# Main call ========================================================= CHJ =====
if __name__=='__main__':
    main()

