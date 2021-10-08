###################################################################### CHJ #####
## Name		: plot_fv3lam_restart.py
## Language	: Python 3.7
## Usage	: Plot restart files for fv3 regional modeling
## Input files  : fv_XXX.res.tile1.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2021/10/06: Chan-Hoo Jeon : Preliminary version
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
dnm_data="/scratch2/NCEPDEV/stmp3/Chan-hoo.Jeon/expt_dirs/test_da_no2_gsdhrrr25_1day/2019080112/RESTART/"

# Domain name
domain_nm='GSD_HRRR_25km'

# restart file name
fnm_input_core='fv_core.res.tile1.nc'
fnm_input_swnd='fv_srf_wnd.res.tile1.nc'
fnm_input_trcr='fv_tracer.res.tile1.nc'

# Variables
# fv_core.res.tile1.nc: u, v, W, DZ, T, delp, phis
# fv_srf_wnd.res.tile1.nc: u_srf, v_srf
# fv_tracer.res.tile1.nc: tracers
vars_res=["T","no2","pm25co"]


# Vertical layer number (3d only)
lvl=1

# basic forms of title and file name
out_title_base='FV3LAM::RESTART::'+domain_nm+'::'
out_fname_base='fv3lam_out_restart_'+domain_nm+'_'

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'


# Main part (will be called at the end) ======================= CHJ =====
def main():
# ============================================================= CHJ =====
    global rcore,rswnd,rtrcr
    global lon,lat,lon_w,lat_w,lon_s,lat_s
    global extent,c_lon,c_lat

    print(' ===== INPUT: Grid fron gfs_dat.nc ==================================')
    # open the data file
    dnm_grd=dnm_data+"../INPUT/"
    fnm_grd='gfs_data.nc'
    fname=os.path.join(dnm_grd,fnm_grd)
    try: grdf=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(grdf)

    print(' ===== INPUT: RESTART (core) ========================================')
    # open the data file
    fname=os.path.join(dnm_data,fnm_input_core)
    try: rcore=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(rcore)

    print(' ===== INPUT: RESTART (srf_wnd) =====================================')
    # open the data file
    fname=os.path.join(dnm_data,fnm_input_swnd)
    try: rswnd=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(rswnd)

    print(' ===== INPUT: RESTART (tracer) ======================================')
    # open the data file
    fname=os.path.join(dnm_data,fnm_input_trcr)
    try: rtrcr=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(rtrcr)

    lon=np.ma.masked_invalid(grdf["geolon"].data)
    lat=np.ma.masked_invalid(grdf["geolat"].data)
    lon_w=np.ma.masked_invalid(grdf["geolon_w"].data)
    lat_w=np.ma.masked_invalid(grdf["geolat_w"].data)
    lon_s=np.ma.masked_invalid(grdf["geolon_s"].data)
    lat_s=np.ma.masked_invalid(grdf["geolat_s"].data)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon)
    lon_max=np.max(lon)
    lat_min=np.min(lat)
    lat_max=np.max(lat)

    print(' lon_min=',lon_min,', lon_max=',lon_max)
    print(' lat_min=',lat_min,', lat_max=',lat_max)

    # Plot extent
    esp=1
    extent=[lon_min-esp,lon_max+esp,lat_min-esp,lat_max+esp]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])

    # Variables
    for svar in vars_res:
        res_plot(svar)



# ===== plot ================================================== CHJ =====
def res_plot(svar):
# ============================================================= CHJ =====

    lvlm=lvl-1
    lvls=format(lvl,'03d')

    print(' ===== '+svar+' ===== RESTART ===============================')
    # Extract data array
    if svar=="u" or svar=="v" or svar=="W" or svar=="DZ" or svar=="T" or svar=="delp" or svar=="phis":
        sfld=np.ma.masked_invalid(rcore[svar].data)
    elif svar=="u_srf" or svar=="v_srf":
        sfld=np.ma.masked_invalid(rswnd[svar].data)
    else:
        sfld=np.ma.masked_invalid(rtrcr[svar].data)

    ndim_svar=sfld.ndim

    if ndim_svar==3:
        (nts,nys,nxs)=sfld.shape
        print(' time+2D: nts=',nts,' nys=',nys,' nxs=',nxs)
        sfld2d=np.squeeze(sfld,axis=0)
        out_title_fld=out_title_base+svar
        out_res_fname=out_fname_base+svar
    elif ndim_svar==4:
        (nts,nls,nys,nxs)=sfld.shape
        print(' time+3D: nts=',nts,' nls=',nls,' nys=',nys,' nxs=',nxs)
        sfld2d=sfld[0,lvlm,:,:]
        out_title_fld=out_title_base+svar+'(L='+lvls+')'
        out_res_fname=out_fname_base+svar+'_L'+lvls
    else:
        sys.exit('ERROR: wrong dimension !!!')

    nm_svar=svar
    cs_cmap='jet'
    lb_ext='both'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=2
    cmap_range='round'
    lonf=lon
    latf=lat
 
    if svar=="u":
        lonf=lon_s
        latf=lat_s
        cs_cmap="rainbow"
    elif svar=="v":
        lonf=lon_w
        latf=lat_w
        cmap_range="symmetry"
    elif svar=="T":
        n_rnd=2
    elif svar=="u_srf":
        cmap_range="symmetry"
    elif svar=="v_srf":
        cmap_range="symmetry"
    elif svar=="no2":
        cmap_range="real"
        cs_cmap="gist_ncar_r"
    elif svar=="pm25co":
        cmap_range="round"
        n_rnd=4
        cs_cmap="rainbow"

    print(' RESTART field=',nm_svar)

    # Max and Min of the field
    fmax=np.max(sfld2d)
    fmin=np.min(sfld2d)
    print(' fld_max=',fmax)
    print(' flx_min=',fmin)

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
        cs_min=-10.0
        cs_max=10.0
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_max=',cs_max)
    print(' cs_min=',cs_min)


    # Plot field
    if domain_nm[:7]=='RRFS_NA':
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(
                            central_longitude=-107,central_latitude=53)))
    else:
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
        ax.set_extent(extent, ccrs.PlateCarree())

    back_plot(ax)
    ax.set_title(out_title_fld,fontsize=9)
    cs=ax.pcolormesh(lonf,latf,sfld2d,cmap=cs_cmap,rasterized=True,
        vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(nm_svar,fontsize=6)

    # Output figure
    ndpi=300
    out_file(out_res_fname,ndpi)

  

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

