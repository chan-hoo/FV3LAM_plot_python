###################################################################### CHJ #####
## Name		: plot_fv3lam_comp2f.py
## Language	: Python 3.7
## Usage	: Compare two NetCDF or GRIB2 files for fv3 regional modeling
## Input files  : NetCDF(.nc) or GIRB2(.grib2) files
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/05/07: Chan-Hoo Jeon : Preliminary version
## V001: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
## V002: 2020/10/20: Chan-Hoo Jeon : Add opt. for orography
## V003: 2021/01/20: Chan-Hoo Jeon : Modify opt. for gfs_data
## V004: 2021/03/05: Chan-Hoo Jeon : Simplify the script
## V005: 2021/04/09: Chan-Hoo Jeon : Add individual plots
## V006: 2021/04/20: Chan-Hoo Jeon : Add relative error plot
###################################################################### CHJ #####

import os, sys
import pygrib
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

# Path to the directories where the input files are located.
#dnm_in1="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/expt_dirs/test_community/2020122700/"
dnm_in1="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/debug_ic/new_ic_grb2"
dnm_in2="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/debug_ic/new_ic_netcdf"

# Input file name
fnm_in1='gfs_data_grb2.nc'
fnm_in2='gfs_data_netcdf.nc'

print(fnm_in1[-5:])

# Variables
#vars_comp=["ps","delp"]
#vars_comp=["sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel","ice_aero","liq_aero"]
vars_comp=["u_w","v_w","u_s","v_s"]

if fnm_in1[-2:]=='nc':
    ftype=1
elif fnm_in1[-5:]=='grib2':
    ftype=2
    grb_name='2 metre temperature'
    grb_typlvl='heightAboveGround'
else:
    sys.exit('ERROR: wrong data type !!!')

ilvl=1
ilvlm=ilvl-1

n_rnd=3

# Basic forms of output file and title
out_fname_base='fv3lam_comp_'
out_title_base='COMP::'

# Colormap range option ('symmetry','round','real','fixed')
cmap_range_org='round'
cmap_range='symmetry'
cmap_range_err='round'

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'


# Main part (will be called at the end) ======================= CHJ =====
def main():
# ============================================================= CHJ =====
    global compf1,compf2

    print(' ===== INPUT files ==========================================')
    # open the data file
    if ftype==1:
        fname=os.path.join(dnm_in1,fnm_in1)
        try: compf1=xr.open_mfdataset(fname,**mfdt_kwargs)
        except: raise Exception('Could NOT find the file',fname)
        print(compf1)
        fname=os.path.join(dnm_in2,fnm_in2)
        try: compf2=xr.open_mfdataset(fname,**mfdt_kwargs)
        except: raise Exception('Could NOT find the file',fname)
        print(compf2)

    elif ftype==2:
        fname=os.path.join(dnm_in1,fnm_in1)
        try: compf1=pygrib.open(fname)
        except: raise Exception('Could NOT find the file',fname)
        fname=os.path.join(dnm_in2,fnm_in2)
        try: compf2=pygrib.open(fname)
        except: raise Exception('Could NOT find the file',fname)

    # Variables
    for svar in vars_comp:
        comp_plot(svar)


# ===== plot ==================================================== CHJ =====
def comp_plot(svar):
# =============================================================== CHJ =====

    # Extract data array
    if ftype==1:
        print(' ===== '+svar+' ===== File 1  ===============================')
        sfld1=np.ma.masked_invalid(compf1[svar].data)
        if sfld1.ndim==2:
            (nys1,nxs1)=sfld1.shape
            print(' File 1: 2D: nys=',nys1,' nxs=',nxs1)
            sfld2d1=sfld1
        elif sfld1.ndim==3:
            (nts1,nys1,nxs1)=sfld1.shape
            print(' File 1: time(level)+2D: nts=',nts1,' nys=',nys1,' nxs=',nxs1)
            sfld2d1=np.squeeze(sfld1[ilvlm,:,:])

        print(fnm_in1[0:8])
        if fnm_in1[0:3]=='oro':
            lon=np.ma.masked_invalid(compf1["geolon"].data)
            lat=np.ma.masked_invalid(compf1["geolat"].data)
        elif fnm_in1[0:8]=='gfs_data':
            if svar=='u_w' or svar=='v_w':
                lon=np.ma.masked_invalid(compf1["geolon_w"].data)
                lat=np.ma.masked_invalid(compf1["geolat_w"].data)
            elif svar=='u_s' or svar=='v_s':
                lon=np.ma.masked_invalid(compf1["geolon_s"].data)
                lat=np.ma.masked_invalid(compf1["geolat_s"].data)
            else:
                lon=np.ma.masked_invalid(compf1["geolon"].data)
                lat=np.ma.masked_invalid(compf1["geolat"].data)
        else:
            lon=np.ma.masked_invalid(compf1["lon"].data)
            lat=np.ma.masked_invalid(compf1["lat"].data)
 
        print(' ===== '+svar+' ===== File 2  ===============================')
        sfld2=np.ma.masked_invalid(compf2[svar].data)
        if sfld2.ndim==2:
            (nys2,nxs2)=sfld2.shape
            print(' File 2: 2D: nys=',nys2,' nxs=',nxs2)
            sfld2d2=sfld2
        elif sfld2.ndim==3:
            (nts2,nys2,nxs2)=sfld2.shape
            print(' File 2: time(level)+2D: nts=',nts2,' nys=',nys2,' nxs=',nxs2)
            sfld2d2=np.squeeze(sfld2[ilvlm,:,:])

        if nys1!=nys2 or nxs1!=nxs2:
            sys.exit('ERROR: array size mismatched!!!')

        out_title_fld=out_title_base+svar
        out_comp_fname=out_fname_base+svar

    elif ftype==2:
        print(' ===== '+svar+' ===== File 1  ===============================')
        grbv=compf1.select(name=grb_name,typeOfLevel=grb_typlvl)[ilvlm]
        sfld2d1=grbv.values
        stnm=grbv.shortName
        lat,lon=grbv.latlons()
       
        print(' ===== '+svar+' ===== File 2  ===============================')
        grbv=compf2.select(name=grb_name,typeOfLevel=grb_typlvl)[ilvlm]
        sfld2d2=grbv.values

        out_title_fld=out_title_base+svar
        out_comp_fname=out_fname_base+stnm

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon)
    lon_max=np.max(lon)
    lat_min=np.min(lat)
    lat_max=np.max(lat)

    print(' lon_min=',lon_min,', lon_max=',lon_max)
    print(' lat_min=',lat_min,', lat_max=',lat_max)

    # Plot extent
    extent=[lon_min-5,lon_max+5,lat_min-5,lat_max+3]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])

    # Difference
    svcomp=sfld2d1-sfld2d2
    nm_svar='\u0394'+svar

    lb_ext='both'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3

    f1_max=np.max(sfld2d1)
    f1_min=np.min(sfld2d1)
    print(' fld1_max=',f1_max)
    print(' fld1_min=',f1_min)
    f2_max=np.max(sfld2d2)
    f2_min=np.min(sfld2d2)
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

    out_title_fld_1 = out_title_fld+"::Data 1"
    out_title_fld_2 = out_title_fld+"::Data 2"


    # Plot field: DATA 1
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title_fld_1,fontsize=9)
    cs=ax.pcolormesh(lon,lat,sfld2d1,cmap=cs_cmap_org,rasterized=True,
        vmin=cs_min_12,vmax=cs_max_12,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(svar,fontsize=8)
    # Output figure
    out_comp_fname_1 = out_comp_fname+"_dat1"
    ndpi=300
    out_file(out_comp_fname_1,ndpi)

    # Plot field: DATA 2
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title_fld_2,fontsize=9)
    cs=ax.pcolormesh(lon,lat,sfld2d2,cmap=cs_cmap_org,rasterized=True,
        vmin=cs_min_12,vmax=cs_max_12,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(svar,fontsize=8)
    # Output figure
    out_comp_fname_2 = out_comp_fname+"_dat2"
    ndpi=300
    out_file(out_comp_fname_2,ndpi)



# ===== Difference plots ======================================== CHJ ===== 
    print(' COMP. field=',nm_svar)

    cs_cmap='seismic'
    # Max and Min of the field
    fmax=np.max(svcomp)
    fmin=np.min(svcomp)
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
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title_fld,fontsize=9)
    cs=ax.pcolormesh(lon,lat,svcomp,cmap=cs_cmap,rasterized=True,
        vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(nm_svar,fontsize=8)

    # Output figure
    ndpi=300
    out_file(out_comp_fname,ndpi)


# ===== Relative error plot ========================================= CHJ =====

    err_rel=svcomp/sfld2d1*100
    err_rel[sfld2d1==0.0]=0.0
    err_rel=np.absolute(err_rel)
    out_title_fld_3 = out_title_fld+"::Relative error (%)"
    out_label = nm_svar+": relative error (%)"

    cs_cmap='gist_ncar_r'
    # Max and Min of the field
    fmax_e=np.max(err_rel)
    fmin_e=np.min(err_rel)
    print(' fld_comp_err_max=',fmax_e)
    print(' fld_comp_err_min=',fmin_e)

    # Make the colormap range symmetry
    print(' cmap range for relative error map=',cmap_range_err)
    if cmap_range_err=='symmetry':
        tmp_cmp=max(abs(fmax_e),abs(fmin_e))
        cs_min_e=round(-tmp_cmp,n_rnd)
        cs_max_e=round(tmp_cmp,n_rnd)
    elif cmap_range_err=='round':
        cs_min_e=round(fmin_e,n_rnd)
        cs_max_e=round(fmax_e,n_rnd)
    elif cmap_range_err=='real':
        cs_min_e=fmin_e
        cs_max_e=fmax_e
    elif cmap_range_err=='fixed':
        cs_min_e=-6.0
        cs_max_e=6.0
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')
 
    if cs_min_e==cs_max_e:
        cs_min_e=cs_min_e-0.1
        cs_max_e=cs_max_e+0.1

    print(' cs_max_err=',cs_max_e)
    print(' cs_min_err=',cs_min_e)


    # Plot field
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title_fld_3,fontsize=9)
    cs=ax.pcolormesh(lon,lat,err_rel,cmap=cs_cmap,rasterized=True,
        vmin=cs_min_e,vmax=cs_max_e,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend='max')
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(out_label,fontsize=8)

    # Output figure
    out_comp_fname_3 = out_comp_fname+"_err"
    ndpi=300
    out_file(out_comp_fname_3,ndpi)
 



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

