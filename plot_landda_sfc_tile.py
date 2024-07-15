#!/usr/bin/env python3

###################################################################### CHJ #####
## Name		: plot_landda_sfc_tile.py
## Language	: Python 3.7
## Usage	: Plot sfc output of land-DA workflow analysis task
## Input files  : sfc_data.tile#.nc
## NOAA/EPIC
## History ===============================
## V000: 2024/06/14: Chan-Hoo Jeon : Preliminary version
## V001: 2024/07/15: Chan-Hoo Jeon : add comparison plot
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

# SFC_1 data:
path_sfc1="/scratch2/NAGAPE/epic/Chan-hoo.Jeon/global-workflow_test/comroot/snow_gdas/gdas.20211220/12/model_data/atmos/restart"
fn_sfc1_base='20211220.180000.sfc_data.tile'
fn_sfc1_ext='.nc'
#path_sfc1="/scratch2/NAGAPE/epic/Chan-hoo.Jeon/global-workflow_test/stmp/RUNDIRS/snow_gsi/sfcanl.1388101"
#fn_sfc1_base='fnbgsi.00'
#fn_sfc1_ext=''
#path_sfc1="/scratch2/NAGAPE/epic/Chan-hoo.Jeon/landda_test/ptmp/test/tmp/DATA_SHARE/20000103"
#fn_sfc1_base='20000103.000000.sfc_data.tile'
#fn_sfc1_ext='.nc'

sfc1_out_txt='before'

# SFC_2 data:
#path_sfc2='/scratch2/NAGAPE/epic/Chan-hoo.Jeon/global-workflow_test/comroot/snow_gdas/gdas.20211220/18/analysis/snow'
path_sfc2='/scratch2/NAGAPE/epic/Chan-hoo.Jeon/global-workflow_test/stmp/RUNDIRS/snow_gdas/gdassnowanl_18/anl'
fn_sfc2_base='20211220.180000.sfc_data.tile'
fn_sfc2_ext='.nc'
#path_sfc2=path_sfc1
#fn_sfc2_base='fnbgso.00'
#fn_sfc2_ext=''
#path_sfc2="/scratch2/NAGAPE/epic/Chan-hoo.Jeon/landda_test/ptmp/test/com/landda/v1.2.1/landda.20000103"
#fn_sfc2_base='20000103.000000.sfc_data.tile'
#fn_sfc2_ext='.nc'

sfc2_out_txt='after'

# increment:
opt_analysis='jedi'
path_xainc=path_sfc2
fn_xainc_base='snowinc.20211220.180000.sfc_data.tile'
#fn_xainc_base='20000103.000000.xainc.sfc_data.tile'
fn_xainc_ext='.nc'

# variable
#sfc_var_nm='snwdph'  # land-DA (jedi)
sfc_var_nm='snodl'   # global-wflow (gdas)
#sfc_var_nm='stc'
#sfc_var_nm='smc'
#sfc_var_nm='slc'

# Vertical layer number
zlyr_num=2

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
    sfc1_data=get_sfc(path_sfc1,fn_sfc1_base,fn_sfc1_ext,sfc_var_nm,sfc1_out_txt,0)
    print('data 1: ',sfc1_data.shape)
    # get sfc data from dir2
    sfc2_data=get_sfc(path_sfc2,fn_sfc2_base,fn_sfc2_ext,sfc_var_nm,sfc2_out_txt,0)
    print('data 2: ',sfc2_data.shape)

    (ntile,ny,nx)=sfc1_data.shape
    sfc_xainc_data=np.zeros([ntile,ny,nx])
    print('zero array: ',sfc_xainc_data.shape)
    opt_inc=0
    # get sfc increment data from dir2
    if sfc_var_nm == 'snodl' or sfc_var_nm == 'snwdph':
        if opt_analysis == 'jedi':
            sfc_xainc_data=get_sfc(path_xainc,fn_xainc_base,fn_xainc_ext,sfc_var_nm,'xainc',0)
            print('inc: ',sfc_xainc_data.shape)
            opt_inc=1

    # compare sfc1 and sfc2
    compare_sfc(sfc1_data,sfc2_data,sfc_xainc_data,opt_inc)
        

# geo lon/lat from orography ======================================== CHJ =====
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
def get_sfc(path_sfc,fn_sfc_base,fn_sfc_ext,sfc_var_nm,sfc_out_txt,ref_opt):
# =================================================================== CHJ =====

    print(' ===== sfc files: '+sfc_var_nm+' :: '+sfc_out_txt+' ========================')
    # open the data file
    sfc_data_all=[]
    if ref_opt == 1:
        snowxy_all=[]
        slmsk_all=[]

    for it in range(num_tiles):
        itp=it+1
        fn_sfc=fn_sfc_base+str(itp)+fn_sfc_ext
        fp_sfc=os.path.join(path_sfc,fn_sfc)

        try: sfc=xr.open_dataset(fp_sfc)
        except: raise Exception('Could NOT find the file',fp_sfc)
#        if it == 1:
#            print(sfc)

        # Extract variable
        sfc_data=np.ma.masked_invalid(sfc[sfc_var_nm].data)
        if sfc_out_txt == 'xainc':
            sfc_data2d=np.squeeze(sfc_data,axis=(0,1))
            cmap_range_opt='symmetry'
            cs_cmap='seismic'
            nm_svar='\u0394'+sfc_var_nm
        else:
            if sfc_var_nm == 'stc' or sfc_var_nm == 'smc' or sfc_var_nm == 'slc':
                sfc_data3d=np.squeeze(sfc_data,axis=0)
                sfc_data2d=sfc_data3d[zlyr_num,:,:]
            else:
                sfc_data2d=np.squeeze(sfc_data,axis=0)

            cmap_range_opt='fixed'
            cs_cmap='gist_ncar_r'
            nm_svar=sfc_var_nm

        sfc_data_all.append(sfc_data2d[None,:])

        if ref_opt == 1:
            # Number of snow layers
            snowxy_orig=np.ma.masked_invalid(sfc['snowxy'].data)
            snowxy2d=np.squeeze(snowxy_orig,axis=0)
            snowxy_all.append(snowxy2d[None,:])
            # sea-land mask
            slmsk_orig=np.ma.masked_invalid(sfc['slmsk'].data)
            slmsk2d=np.squeeze(slmsk_orig,axis=0)
            slmsk_all.append(slmsk2d[None,:])

    sfc_var=np.vstack(sfc_data_all)
    if ref_opt == 1:
        snowxy=np.vstack(snowxy_all)
        snowxy_max=np.max(snowxy)
        snowxy_min=np.min(snowxy)
        print('snowxy_max=',snowxy_max)
        print('snowxy_min=',snowxy_min)
        slmsk=np.vstack(slmsk_all)
        slmsk_max=np.max(slmsk)
        slmsk_min=np.min(slmsk)
        print('slmsk_max=',slmsk_max)
        print('slmsk_min=',slmsk_min)

    sfc_var_max=np.max(sfc_var)
    sfc_var_min=np.min(sfc_var)
    print('sfc_var_max=',sfc_var_max)
    print('sfc_var_min=',sfc_var_min)
    sfc_var_min08=sfc_var_min*0.8
    sfc_var_max08=sfc_var_max*0.8
    print('sfc_var_max08=',sfc_var_max08)
    print('sfc_var_min08=',sfc_var_min08)

    if cmap_range_opt=='symmetry':
        n_rnd=0
        tmp_cmp=max(abs(sfc_var_max08),abs(sfc_var_min08))
        cs_min=round(-tmp_cmp,n_rnd)
        cs_max=round(tmp_cmp,n_rnd)
        cbar_extend='both'
    elif cmap_range_opt=='round':
        n_rnd=0
        cs_min=round(sfc_var_min08,n_rnd)
        cs_max=round(sfc_var_max08,n_rnd)
        cbar_extend='both'
    elif cmap_range_opt=='real':
        cs_min=sfc_var_min
        cs_max=sfc_var_max
        cbar_extend='neither'
    elif cmap_range_opt=='fixed':
        cs_min=0.0
        cs_max=150.0
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
   
    if ref_opt == 1:
        for it in range(num_tiles):
            glon_tile=np.squeeze(glon[it,:,:])
            if it == 0:
                glon_tile=(glon_tile+180)%360-180
            glat_tile=np.squeeze(glat[it,:,:])
            c_glon=np.round(np.mean(glon_tile),decimals=2)
            c_glat=np.round(np.mean(glat_tile),decimals=2)
            print(glon_tile.shape)
            print(glat_tile.shape)
            print("c_glon, c_glat for tile",str(it+1),"=",c_glon,c_glat)

            snowxy_tile=np.squeeze(snowxy[it,:,:])
            slmsk_tile=np.squeeze(slmsk[it,:,:])  
            print(snowxy_tile.shape)
            print(slmsk_tile.shape)       
            print("glon:  max=",np.max(glon_tile)," min=",np.min(glon_tile))
            print("glat:  max=",np.max(glat_tile)," min=",np.min(glat_tile))
            print("snowxy:max=",np.max(snowxy_tile)," min=",np.min(snowxy_tile))
            print("slmsk: max=",np.max(slmsk_tile), " min=",np.min(slmsk_tile))
            # plot snow layers
            out_title=out_title_base+'snowxy::Tile '+str(it+1)
            out_fn=out_fn_base+'snowxy_tile'+str(it+1)
            fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(c_glon,c_glat)))

            ax.set_title(out_title, fontsize=6)
            back_plot(ax)
            cs=ax.pcolormesh(glon_tile,glat_tile,snowxy_tile,cmap=plt.cm.get_cmap('Set1',3),
                rasterized=True,vmin=-2,vmax=0,transform=ccrs.PlateCarree())
            divider=make_axes_locatable(ax)
            ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            cbar=plt.colorbar(cs,cax=ax_cb,extend='neither')
            cbar.ax.tick_params(labelsize=6)
            cbar.set_label('snowxy',fontsize=6)
            out_file(out_fn,300)

            # plot sea-land mask
            out_title=out_title_base+'slmsk::Tile '+str(it+1)
            out_fn=out_fn_base+'slmsk_tile'+str(it+1)
            fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(c_glon,c_glat)))
            ax.set_title(out_title, fontsize=6)
            back_plot(ax)
            cs=ax.pcolormesh(glon[it,:,:],glat[it,:,:],slmsk_tile,cmap=plt.cm.get_cmap('Paired',2),
                rasterized=True,vmin=0,vmax=1,transform=ccrs.PlateCarree())
            divider=make_axes_locatable(ax)
            ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            cbar=plt.colorbar(cs,cax=ax_cb,extend='neither')
            cbar.ax.tick_params(labelsize=6)
            cbar.set_label('slmsk',fontsize=6)
            out_file(out_fn,300) 
    return sfc_var


# Compare two data set and plot ===================================== CHJ =====
def compare_sfc(sfc_data1,sfc_data2,inc_data,opt_inc):
# =================================================================== CHJ =====
    print(' ===== compare files =============================================')
    print(' data 1: ',sfc_data1.shape)
    print(' data 2: ',sfc_data2.shape)

    diff_data=sfc_data2-sfc_data1
    print(' diff. data: ',diff_data.shape)

    diff_data_max=np.max(diff_data)
    diff_data_min=np.min(diff_data)
    print('diff_data_max=',diff_data_max)
    print('diff_data_min=',diff_data_min)

    diff_data_max08=diff_data_max*0.8
    diff_data_min08=diff_data_min*0.8
    print('diff_data_max08=',diff_data_max08)
    print('diff_data_min08=',diff_data_min08)
    diff_data_cs_max=max(abs(diff_data_max08),abs(diff_data_min08))
    plot_increment(diff_data,'diff_sfc',diff_data_cs_max)

    if opt_inc == 1:
        diff_inc=diff_data-inc_data
        diff_inc_max=np.max(diff_inc)
        diff_inc_min=np.min(diff_inc)
        print('diff_inc_max=',diff_inc_max)
        print('diff_inc_min=',diff_inc_min)

        diff_inc_max08=diff_inc_max*0.8
        diff_inc_min08=diff_inc_min*0.8
        print('diff_inc_max08=',diff_inc_max08)
        print('diff_inc_min08=',diff_inc_min08)
        diff_inc_cs_max=max(abs(diff_inc_max08),abs(diff_inc_min08))
        plot_increment(diff_inc,'diff_inc',diff_inc_cs_max)


# increment/difference plot ========================================== CHJ =====
def plot_increment(plt_var,plt_var_nm,cs_max):
# ==================================================================== CHJ =====

    cs_cmap='seismic'
    nm_svar='\u0394'+plt_var_nm
    n_rnd=0
    cs_min=cs_max*-1.0
    cbar_extend='both'

    out_title=out_title_base+plt_var_nm
    out_fn=out_fn_base+plt_var_nm

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_title(out_title, fontsize=6)
    # Call background plot
    back_plot(ax)

    for it in range(num_tiles):
        cs=ax.pcolormesh(glon[it,:,:],glat[it,:,:],plt_var[it,:,:],cmap=cs_cmap,rasterized=True,
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

