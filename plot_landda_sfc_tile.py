#!/usr/bin/env python3

###################################################################### CHJ #####
## Name		: plot_landda_sfc_tile.py
## Language	: Python 3.7
## Usage	: Plot sfc output of land-DA workflow analysis task
## Input files  : sfc_data.tile#.nc
## NOAA/EPIC
## History ===============================
## V000: 2024/06/14: Chan-Hoo Jeon : Preliminary version
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
machine='orion'

print(' You are on', machine)

#### Machine-specific input data ===================================== CHJ =====
# cartopy.config: Natural Earth data for background
# out_fig_dir: directory where the output files are created

if machine=='hera':
    cartopy.config['data_dir']='/scratch2/NAGAPE/epic/UFS_Land-DA_Dev/inputs/NaturalEarth'
    path_orog="/scratch2/NAGAPE/epic/UFS_Land-DA_Dev/inputs/orog_files"
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig"
elif machine=='orion' or machine=='hercules':
    cartopy.config['data_dir']='/work/noaa/epic/UFS_Land-DA_Dev/inputs/NaturalEarth'
    path_orog="/work/noaa/epic/UFS_Land-DA_Dev/inputs/orog_files"
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig"
else:
    sys.exit('ERROR: Required input data are NOT set in this machine !!!')

# Case-dependent input ============================================== CHJ =====

# SFC_1 data: pre_anal
path_sfc1="/work/noaa/epic/chjeon/landda_test/ptmp/test/tmp/DATA_SHARE/20000104"
fn_sfc1_base='20000104.000000.sfc_data.tile'
sfc1_out_txt='pre_anal'

# SFC_2 data: analysis
path_sfc2="/work/noaa/epic/chjeon/landda_test/ptmp/test/com/landda/v1.2.1/landda.20000104"
fn_sfc2_base='20000104.000000.sfc_data.tile'
fn_xainc_base='20000104.000000.xainc.sfc_data.tile'
sfc2_out_txt='analysis'

# basic forms of title and file name
out_title_base='Land-DA::SFC::'
out_fn_base='landda_out_sfc_'

# Geo files
fn_orog_base='oro_C96.mx100.tile'
# Number of tiles
num_tiles=6
# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'
# center of map
c_lon=-77.0369


# Main part (will be called at the end) ============================= CHJ =====
def main():
# =================================================================== CHJ =====
    # get lon, lat from orography
    get_geo()
    # get sfc data from dir1
    sfc1_data=get_sfc(path_sfc1,fn_sfc1_base,sfc1_out_txt)
    print(sfc1_data.shape)
    # get sfc data from dir2
    sfc2_data=get_sfc(path_sfc2,fn_sfc2_base,sfc2_out_txt)
    print(sfc2_data.shape)
    # get sfc increment data from dir2
    sfc_xainc_data=get_sfc(path_sfc2,fn_xainc_base,'xainc')
    print(sfc_xainc_data.shape)
    # compare and plot sfc1 and sfc2
    comp_sfc(sfc1_data,sfc2_data)


# Orography plot ==================================================== CHJ =====
def get_geo():
# =================================================================== CHJ =====

    global glon,glat

    print(' ===== geo data files ====================================')
    # open the data file
    glon_all=[]
    glat_all=[]
    for it in range(num_tiles):
        itp=it+1
        fn_orog=fn_orog_base+str(itp)+'.nc'
        fp_orog=os.path.join(path_orog,fn_orog)

        try: orog=xr.open_dataset(fp_orog)
        except: raise Exception('Could NOT find the file',fp_orog)
#        print(orog)
        # Extract longitudes, and latitudes
        geolon=np.ma.masked_invalid(orog['geolon'].data)
        geolat=np.ma.masked_invalid(orog['geolat'].data)
        glon_all.append(geolon[None,:])
        glat_all.append(geolat[None,:])

    glon=np.vstack(glon_all)
    glat=np.vstack(glat_all)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(glon)
    lon_max=np.max(glon)
    lat_min=np.min(glat)
    lat_max=np.max(glat)

    out_title=out_title_base+'GEO'
    out_fn=out_fn_base+'geo'

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
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
    out_file(out_fn,ndpi)


# Get sfc_data from files and plot ================================== CHJ =====
def get_sfc(path_sfc,fn_sfc_base,sfc_out_txt):
# =================================================================== CHJ =====

    print(' ===== sfc files :'+sfc_out_txt+' ========================')
    # open the data file
    sfc_data_all=[]
    for it in range(num_tiles):
        itp=it+1
        fn_sfc=fn_sfc_base+str(itp)+'.nc'
        fp_sfc=os.path.join(path_sfc,fn_sfc)

        try: sfc=xr.open_dataset(fp_sfc)
        except: raise Exception('Could NOT find the file',fp_sfc)
#        print(sfc)
        # Extract snow depth
        sfc_var_nm='snwdph'
        sfc_data=np.ma.masked_invalid(sfc[sfc_var_nm].data)
        if sfc_out_txt == 'xainc':
            sfc_data2d=np.squeeze(sfc_data,axis=(0,1))
            cmap_range_opt='symmetry'
            cs_cmap='seismic'
            nm_svar='\u0394'+sfc_var_nm
        else:
            sfc_data2d=np.squeeze(sfc_data,axis=0)
            cmap_range_opt='round'
            cs_cmap='gist_ncar_r'
            nm_svar=sfc_var_nm
        sfc_data_all.append(sfc_data2d[None,:])

    sfc_var=np.vstack(sfc_data_all)

    sfc_var_max=np.max(sfc_var)
    sfc_var_min=np.min(sfc_var)
    print('sfc_var_max=',sfc_var_max)
    print('sfc_var_min=',sfc_var_min)
    sfc_var_min05=sfc_var_min*0.5
    sfc_var_max05=sfc_var_max*0.5
    print('sfc_var_max05=',sfc_var_max05)
    print('sfc_var_min05=',sfc_var_min05)

    if cmap_range_opt=='symmetry':
        n_rnd=0
        tmp_cmp=max(abs(sfc_var_max05),abs(sfc_var_min05))
        cs_min=round(-tmp_cmp,n_rnd)
        cs_max=round(tmp_cmp,n_rnd)
        cbar_extend='both'
    elif cmap_range_opt=='round':
        n_rnd=0
        cs_min=round(sfc_var_min05,n_rnd)
        cs_max=round(sfc_var_max05,n_rnd)
        cbar_extend='both'
    elif cmap_range_opt=='real':
        cs_min=sfc_var_min
        cs_max=sfc_var_max
        cbar_extend='neither'
    elif cmap_range_opt=='fixed':
        cs_min=-5.0
        cs_max=5.0
        cbar_extend='both'
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print('cs_max=',cs_max)
    print('cs_min=',cs_min)

    out_title=out_title_base+sfc_out_txt
    out_fn=out_fn_base+sfc_out_txt

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_title(out_title, fontsize=6)
    # Call background plot
    back_plot(ax)

    for it in range(num_tiles):
        cs=ax.pcolormesh(glon[it,:,:],glat[it,:,:],sfc_var[it,:,:],cmap=cs_cmap,rasterized=True,
            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())

    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=cbar_extend)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(nm_svar,fontsize=6)

    # Output figure
    ndpi=300
    out_file(out_fn,ndpi)
    
    return sfc_var


# Compare and plot sfc_data ========================================= CHJ =====
def comp_sfc(sfc1_data,sfc2_data):
# =================================================================== CHJ =====

    print(' ===== compare sfc_data files ========================')




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

