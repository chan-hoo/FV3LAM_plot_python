#!/usr/bin/env python3

###################################################################### CHJ #####
## Name		: plot_landda_restart.py
## Language	: Python 3.7
## Usage	: Plot restart output file of land-DA workflow
## Input files  : ufs_land_restart.tile#.nc
## NOAA/EPIC
## History ===============================
## V000: 2024/09/26: Chan-Hoo Jeon : Preliminary version
###################################################################### CHJ #####

import os, sys
import numpy as np
import netCDF4 as nc
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig"
elif machine=='orion' or machine=='hercules':
    cartopy.config['data_dir']='/work/noaa/epic/UFS_Land-DA_Dev/inputs/NaturalEarth'
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig"
else:
    sys.exit('ERROR: Required input data are NOT set in this machine !!!')

# Case-dependent input ============================================== CHJ =====

# SFC_1 data:
path_data='/scratch2/NCEPDEV/naqfc/Chan-hoo.Jeon/landda_test/ptmp/test/com/landda/v2.0.0/landda.20000104'
fn_data_date='2000-01-05_00-00-00'
fn_data_base='ufs_land_restart.'+fn_data_date+'.tile'
fn_data_ext='.nc'
# variable
var_list=['snwdph','smc']
# soil level number
soil_lvl_num=1
# basic forms of title and file name
out_title_base='Land-DA::restart::'
out_fn_base='landda_out_restart_'
# Number of tiles
num_tiles=6
# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'
# center of map
c_lon=-77.0369


# Main part (will be called at the end) ============================= CHJ =====
def main():
# =================================================================== CHJ =====
    # get lon, lat
    get_geo()
    # plot restart file
    for var_nm in var_list:
        plot_data(path_data,fn_data_base,fn_data_ext,var_nm)
       

# geo lon/lat from orography ======================================== CHJ =====
def get_geo():
# =================================================================== CHJ =====

    global glon,glat

    print(' ===== geo data files ====================================')

    # open the data file
    for it in range(num_tiles):
        itp=it+1
        fn_data=fn_data_base+str(itp)+fn_data_ext
        fp_data=os.path.join(path_data,fn_data)
        try: data_raw=nc.Dataset(fp_data)
        except: raise Exception('Could NOT find the file',fp_data)
        if itp == 1:
            print(data_raw)
        # Extract geo data
        glon_data=np.ma.masked_invalid(data_raw.variables['grid_xt'])
        print('Dimension of glon(grid_xt)=',glon_data.shape)
        print('Tile',itp,',max:',np.max(glon_data))
        print('Tile',itp,',min:',np.min(glon_data))

        glat_data=np.ma.masked_invalid(data_raw.variables['grid_yt'])
        print('Dimension of glat(grid_yt)=',glat_data.shape)
        print('Tile',itp,',max:',np.max(glat_data))
        print('Tile',itp,',min:',np.min(glat_data))

        if itp == 1:
            ny,nx=glon_data.shape
            glon=np.zeros((num_tiles,ny,nx))
            glat=np.zeros((num_tiles,ny,nx))

        data_raw.close()        
        glon[it,:,:]=glon_data[:,:]
        glat[it,:,:]=glat_data[:,:]

    print('Dimension of glon=',glon.shape)
    print('Dimension of glon=',glat.shape)


# Get sfc_data from files and plot ================================== CHJ =====
def plot_data(path_data,fn_data_base,fn_data_ext,var_nm):
# =================================================================== CHJ =====

    print(' ===== data file: '+var_nm+' ========================')
    # open the data file
    for it in range(num_tiles):
        itp=it+1
        fn_data=fn_data_base+str(itp)+fn_data_ext
        fp_data=os.path.join(path_data,fn_data)
        try: data_raw=nc.Dataset(fp_data)
        except: raise Exception('Could NOT find the file',fp_data)
        # Extract valid variable
        var_data=np.ma.masked_invalid(data_raw.variables[var_nm])
        if var_nm == 'stc' or var_nm == 'smc' or var_nm == 'slc':
            print('Dimension of original data=',var_data.shape)
            var_data_2d=var_data[:,soil_lvl_num-1,:,:]
        else:
            var_data_2d=var_data                
 
        print('Dimension of data=',var_data_2d.shape)
        print('Tile',itp,',max:',np.max(var_data_2d))
        print('Tile',itp,',min:',np.min(var_data_2d))

        if itp == 1:
            plt_var=var_data_2d
        else:
            plt_var=np.ma.concatenate((plt_var,var_data_2d),axis=0)
        data_raw.close()

#    plt_var=np.vstack(var_data_all)
    print('Dimension of data set=',plt_var.shape)

    cs_max=np.max(plt_var)
    cs_min=np.min(plt_var)
    print('cs_max=',cs_max)
    print('cs_min=',cs_min)

    cs_cmap='gist_ncar_r'
    cbar_extend='neither'

    for it in range(num_tiles):
        itp=it+1
        glon_tile=np.squeeze(glon[it,:,:])
        if itp == 1:
            glon_tile=(glon_tile+180)%360-180
        glat_tile=np.squeeze(glat[it,:,:])
        var_tile=np.squeeze(plt_var[it,:,:])
        c_glon=np.round(np.mean(glon_tile),decimals=2)
        c_glat=np.round(np.mean(glat_tile),decimals=2)
        print("c_glon, c_glat for tile",str(it+1),"=",c_glon,c_glat)

        out_title=out_title_base+var_nm+'::Tile'+str(itp)
        out_fn=out_fn_base+var_nm+'_tile'+str(itp)

        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(c_glon,c_glat)))
        ax.set_title(out_title, fontsize=6)
        # Call background plot
        back_plot(ax)

        cs=ax.pcolormesh(glon_tile,glat_tile,var_tile,cmap=cs_cmap,
            rasterized=True,vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
        divider=make_axes_locatable(ax)
        ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cbar=plt.colorbar(cs,cax=ax_cb,extend='neither')
        cbar.ax.tick_params(labelsize=6)
        cbar.set_label(var_nm,fontsize=6)

#    for it in range(num_tiles):
#        cs=ax.pcolormesh(glon[it,:,:],glat[it,:,:],plt_var[it,:,:],cmap=cs_cmap,rasterized=True,
#            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
#    divider=make_axes_locatable(ax)
#    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
#    fig.add_axes(ax_cb)
#    cbar=plt.colorbar(cs,cax=ax_cb,extend=cbar_extend)
#    cbar.ax.tick_params(labelsize=6)
#    cbar.set_label(var_nm,fontsize=6)

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

