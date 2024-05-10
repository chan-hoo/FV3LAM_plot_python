#!/usr/bin/env python3

###################################################################### CHJ #####
## Name		: plot_landda_comp_restart.py
## Language	: Python 3.9
## Usage	: Compare two restart NetCDF files
## Input files  : NetCDF(.nc) files
## NOAA/EPIC
## History ===============================
## V000: 2024/05/10: Chan-Hoo Jeon : Preliminary version
###################################################################### CHJ #####

import os, sys
import pygrib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable

# HPC machine ('hera','orion')
machine='orion'

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
    out_fig_dir="/work/noaa/epic/chjeon/FIG_output/"
    mfdt_kwargs={'parallel':False,'combine':'by_coords'}
else:
    sys.exit('ERROR: Required input data are NOT set !!!')

plt.switch_backend('agg')

# Case-dependent input =============================================== CHJ =====

# Path to the directories where the input files are located.
dnm_in1="/work/noaa/epic/chjeon/landda_test/land-DA_workflow/sorc/build/test/"
dnm_in2="/work/noaa/epic/UFS_Land-DA/test_base/mem000/restarts/vector/"

# Input file name
fnm_in1="ufs_land_restart.2019-12-22_00-00-00.nc"
fnm_in2="ufs_land_restart_back.2019-12-22_00-00-00.nc"

# File for lon/lat
if machine=='hera':
    dnm_geo=""
elif machine=='orion':
    dnm_geo="/work/noaa/epic/UFS_Land-DA/inputs/forcing/era5/init"

fnm_geo="ufs-land_C96_init_2010-12-31_23-00-00.nc"

# Domain name
domain_nm='C96'

# Variables
vars_comp=["snow_depth"]

ilvl=1
ilvlm=ilvl-1

n_rnd=3

# Basic forms of output file and title
out_fname_base='landda_comp_'
out_title_base='COMP::'

# Colormap range option ('symmetry','round','real','fixed')
cmap_range_org='real'
cmap_range='symmetry'
cmap_range_err='round'

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'


# Main part (will be called at the end) ======================= CHJ =====
def main():
# ============================================================= CHJ =====
    global compf1,compf2

    print(' ===== INPUT files ==========================================')
    # Load sea-land mask and lon/lat from 6 tiles 
    load_geo()

    # Open the data file
    fpath=os.path.join(dnm_in1,fnm_in1)
    try: compf1=nc.Dataset(fpath)
    except: raise Exception('Could NOT find the file',fname)
    print(compf1)
    fpath=os.path.join(dnm_in2,fnm_in2)
    try: compf2=nc.Dataset(fpath)
    except: raise Exception('Could NOT find the file',fname)
    print(compf2)

    # Variables
    for svar in vars_comp:
        comp_plot(svar)


# ===== load lon/lat from tiles ================================= CHJ =====
def load_geo():
# =============================================================== CHJ =====
    global lon, lat

    # Open obs file for lon/lat
    fpath=os.path.join(dnm_geo,fnm_geo)
    try: mdat=nc.Dataset(fpath)
    except: raise Exception('Could NOT find the file',fpath)
    print(mdat)

    lon_orig=mdat.variables['longitude']
    lat_orig=mdat.variables['latitude']
    lon=np.squeeze(lon_orig)
    lat=np.squeeze(lat_orig)

    print('lon: dimension=', lon.ndim,' size=', lon.shape)
    print('lat: dimension=', lat.ndim,' size=', lat.shape)

    # Longitude 0:360 => -180:180
    lon_max=np.max(lon)
    lon_min=np.min(lon)
    if lon_max>180:
        print('ORIG: lon_min=',lon_min,', lon_max=',lon_max)
        lon=(lon+180)%360-180
        lon_min=np.min(lon)
        lon_max=np.max(lon)
        print('NEW : lon_min=',lon_min,', lon_max=',lon_max)


# ===== plot ==================================================== CHJ =====
def comp_plot(svar):
# =============================================================== CHJ =====

    # Extract data array
    print(' ===== '+svar+' ===== File 1  ===============================')
    fld1=compf1.variables[svar]
    print('File 1: dimension=',fld1.ndim, ', size=',fld1.shape)
    (ntime1,nloc1)=fld1.shape
    if ntime1==1:
        fld1_loc=np.squeeze(fld1[0,:])
    else:
        fld1_loc=np.squeeze(fld1[ilvlm,:])

    print('Field 1: size=', nloc1)

    print(' ===== '+svar+' ===== File 2  ===============================')
    fld2=compf2.variables[svar]
    print('File 2 dimensions=',fld2.ndim, ', size=',fld2.shape)
    (ntime2,nloc2)=fld2.shape
    if ntime2==1:
        fld2_loc=np.squeeze(fld2[0,:])
    else:
        fld2_loc=np.squeeze(fld2[ilvlm,:])

    if nloc1!=nloc2:
        sys.exit('ERROR: array size mismatched!!!')

    print('Field 2: size=', nloc2)

    out_title_fld=out_title_base+svar
    out_comp_fname=out_fname_base+svar

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon)
    lon_max=np.max(lon)
    lat_min=np.min(lat)
    lat_max=np.max(lat)

    print('lon_min=',lon_min,', lon_max=',lon_max)
    print('lat_min=',lat_min,', lat_max=',lat_max)

    # Plot extent
    extent=[lon_min,lon_max,lat_min,lat_max]
    # for CONUS
#    extent=[-125,-66,23,53]
    print(extent)

#    c_lon=np.mean(extent[:2])
    c_lon=-77.0369 # D.C.
    print(' c_lon=',c_lon)

    # Difference
    fcomp=fld2_loc-fld1_loc
    nm_svar='\u0394'+svar

    lb_ext='both'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    scat_sz=1.5

    f1_max=np.max(fld1_loc)
    f1_min=np.min(fld1_loc)
    print(' fld1_max=',f1_max)
    print(' fld1_min=',f1_min)
    f2_max=np.max(fld2_loc)
    f2_min=np.min(fld2_loc)
    print(' fld2_max=',f2_max)
    print(' fld2_min=',f2_min)

# ===== Individual plots ========================================= CHJ =====
    f12_max=max(f1_max,f2_max)
    f12_min=min(f1_min,f2_min)
    cs_cmap_org='gist_ncar'

    # Make the colormap range symmetry
    print(' cmap range_org=',cmap_range_org)
    if cmap_range_org=='symmetry':
        tmp_cmp=max(abs(f12_max),abs(f12_min))
        cs_min_12=round(-tmp_cmp,n_rnd)
        cs_max_12=round(tmp_cmp,n_rnd)
    elif cmap_range_org=='round':
        cs_min_12=round(f12_min,n_rnd)
        cs_max_12=round(f12_max,n_rnd)
    elif cmap_range_org=='real':
        cs_min_12=f12_min
        cs_max_12=f12_max
    elif cmap_range_org=='fixed':
        cs_min_12=-6.0
        cs_max_12=6.0
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    if cs_min_12==cs_max_12:
        cs_min_12=cs_min_12-0.1
        cs_max_12=cs_max_12+0.1

    print(' cs_max_org=',cs_max_12)
    print(' cs_min_org=',cs_min_12)

    # Plot field: DATA 1
    out_title_fld_1=out_title_fld+"::F1"
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    back_plot(ax)
    ax.set_title(out_title_fld_1,fontsize=9)
    cs=ax.scatter(lon,lat,transform=ccrs.PlateCarree(),c=fld1_loc,cmap=cs_cmap_org,
                  vmin=cs_min_12,vmax=cs_max_12,s=scat_sz)
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(svar,fontsize=6)
    # Output figure
    out_comp_fname_1=out_comp_fname+"_f1"
    ndpi=300
    out_file(out_comp_fname_1,ndpi)

    # Plot field: DATA 2
    out_title_fld_2=out_title_fld+"::F2"
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    back_plot(ax)
    ax.set_title(out_title_fld_2,fontsize=9)
    cs=ax.scatter(lon,lat,transform=ccrs.PlateCarree(),c=fld2_loc,cmap=cs_cmap_org,
                  vmin=cs_min_12,vmax=cs_max_12,s=scat_sz)
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(svar,fontsize=6)
    # Output figure
    out_comp_fname_2=out_comp_fname+"_f2"
    ndpi=300
    out_file(out_comp_fname_2,ndpi)


# ===== Difference plots ======================================== CHJ ===== 
    print(' COMP. field=',nm_svar)

    out_title_fld_3=out_title_fld+"::F1-F2"
    cs_cmap='seismic'
    # Max and Min of the field
    fmax=np.max(fcomp)
    fmin=np.min(fcomp)
    print(' fld_comp_max=',fmax)
    print(' fld_comp_min=',fmin)

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
        cs_min=-6.0
        cs_max=6.0
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    if cs_min==cs_max:
        cs_min=cs_min-0.1
        cs_max=cs_max+0.1

    print(' cs_max=',cs_max)
    print(' cs_min=',cs_min)

    # Plot field
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    back_plot(ax)
    ax.set_title(out_title_fld_3,fontsize=9)
    cs=ax.scatter(lon,lat,transform=ccrs.PlateCarree(),c=fcomp,cmap=cs_cmap,
                  vmin=cs_min,vmax=cs_max,s=scat_sz)
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(nm_svar,fontsize=6)

    # Output figure
    ndpi=300
    out_file(out_comp_fname,ndpi)


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

