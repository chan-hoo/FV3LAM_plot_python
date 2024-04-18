###################################################################### CHJ #####
## Name		: plot_fv3lam_dynf.py
## Language	: Python 3.7
## Usage	: Plot an output, dyn, for fv3 regional modeling
## Input files  : dynfXXX.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/05/06: Chan-Hoo Jeon : Preliminary version
## V001: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
## V002: 2021/03/05: Chan-Hoo Jeon : Simplify the script
## V003: 2021/06/24: Chan-Hoo Jeon : Add a projection for RRFS_NA domain
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
dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/ufs_srw_app/srw_dev_test/expt_dirs/grid_RRFS_NA_13km/2019070100"

# Domain name
domain_nm='RRFS_NA_13km'

# grid file name
fnm_hr='f003'

fnm_input='dyn'+fnm_hr+'.nc'

# Variables
#vars_dyn=["cld_amt","clwmr","delz","dnvvelmax","dpres"]
#vars_dyn=["dzdt","grle","hgtsfc","icmr"]
#vars_dyn=["maxvort01"]
#vars_dyn=["o3mr","pressfc","rwmr","snmr"]
vars_dyn=["tmp"]
#vars_dyn=["uhmax03","uhmax25","uhmin03","uhmin25","upvvelmax"]
#vars_dyn=["ugrd","vgrd"]
#vars_dyn=["hgtsfc"]

# Vertical layer number (3d only)
lvl=1

# basic forms of title and file name
out_title_base='FV3LAM::DYN::'+domain_nm+'::'
out_fname_base='fv3lam_out_dyn_'+domain_nm+'_'

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'


# Main part (will be called at the end) ======================= CHJ =====
def main():
# ============================================================= CHJ =====
    global dynf,lon,lat
    global extent,c_lon,c_lat

    print(' ===== OUTPUT: dyn =======================================')
    # open the data file
    fname=os.path.join(dnm_data,fnm_input)
    try: dynf=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(dynf)

    lon=np.ma.masked_invalid(dynf["lon"].data)
    lat=np.ma.masked_invalid(dynf["lat"].data)

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
    for svar in vars_dyn:
        dyn_plot(svar)



# ===== plot ================================================== CHJ =====
def dyn_plot(svar):
# ============================================================= CHJ =====

    lvlm=lvl-1
    lvls=format(lvl,'03d')


    print(' ===== '+svar+' ===== dyn ===============================')
    # Extract data array
    sfld=np.ma.masked_invalid(dynf[svar].data)

    ndim_svar=sfld.ndim

    if ndim_svar==3:
        (nts,nys,nxs)=sfld.shape
        print(' time+2D: nts=',nts,' nys=',nys,' nxs=',nxs)
        sfld2d=np.squeeze(sfld,axis=0)
        out_title_fld=out_title_base+svar+'::'+fnm_hr
        out_dyn_fname=out_fname_base+svar+'_'+fnm_hr
    elif ndim_svar==4:
        (nts,nls,nys,nxs)=sfld.shape
        print(' time+3D: nts=',nts,' nls=',nls,' nys=',nys,' nxs=',nxs)
        sfld2d=sfld[0,lvlm,:,:]
        out_title_fld=out_title_base+svar+'(L='+lvls+')::'+fnm_hr
        out_dyn_fname=out_fname_base+svar+'_L'+lvls+'_'+fnm_hr
    else:
        sys.exit('ERROR: wrong svar !!!')

    nm_svar=svar
    cs_cmap='jet'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=2
    cmap_range='round'
 
    if svar=="cld_amt":
        nm_svar='Cloud amount'
        n_rnd=5
    elif svar=="clwmr":
        nm_svar='Cloud-water mixing ratio'
        n_rnd=5
    elif svar=="grle":
        nm_svar='Graupel mixing ratio'
        n_rnd=5
    elif svar=="hgtsfc":
        nm_svar='Surface height'
        cs_cmap='gist_ncar_r'
        n_rnd=0
    elif svar=="icmr":
        nm_svar='Cloud-ice mixing ratio'
        n_rnd=4
    elif svar=="maxvort01":
        n_rnd=5
    elif svar=="maxvort02":
        n_rnd=5
    elif svar=="maxvorthy1":
        n_rnd=5
    elif svar=="o3mr":
        nm_svar='Ozone mixing ratio'
        n_rnd=7
    elif svar=="pressfc":
        nm_svar='Surface air pressure'
        n_rnd=0
    elif svar=="rwmr":
        nm_svar='Rain-water mixing ratio'
        n_rnd=5
    elif svar=="snmr":
        nm_svar='Snow-water mixing ratio'
        n_rnd=5
    elif svar=="spfh":
        nm_svar='Specific humidity'
        n_rnd=3
    elif svar=="srh01":
        n_rnd=2
    elif svar=="srh03":
        n_rnd=5
    elif svar=="tmp":
        nm_svar='Temperature'
        n_rnd=1
    elif svar=="ustm":
        nm_svar='u-comp. of storm motion'
        n_rnd=2
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="vstm":  
        nm_svar='v-comp. of storm motion'
        n_rnd=2
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'

    print(' DYN field=',nm_svar)

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
    cs=ax.pcolormesh(lon,lat,sfld2d,cmap=cs_cmap,rasterized=True,
        vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(nm_svar,fontsize=6)

    # Output figure
    ndpi=300
    out_file(out_dyn_fname,ndpi)

  

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

