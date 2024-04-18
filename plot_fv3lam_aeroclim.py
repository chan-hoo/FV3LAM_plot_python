###################################################################### CHJ #####
## Name		: plot_fv3lam_aeroclim.py
## Language	: Python 3.7
## Usage	: Plot input MERRA2 aerosol climatology data
## Input files  : aeroclim.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2021/12/08: Chan-Hoo Jeon : Preliminary version
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
dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/ufs_srw_app/srw_dev_test/expt_dirs/grid_RRFS_CONUS_3km_ics_FV3GFS_lbcs_FV3GFS_suite_GFS_v15_thompson_mynn_lam3km/2019070100/"

# month
mnth=1

mnths=format(mnth,'02d')

# vertical level
lvl=1


# grid file name
fnm_input='aeroclim.m'+mnths+'.nc'

# variables
#vars_atm=["AIRDENS","BCPHILIC","BCPHOBIC","OCPHILIC","OCPHOBIC","PS","SO4"]
vars_atm=["PS"]

# basic forms of title and file name: base+static field name
out_title_base='FV3LAM::MERRA2::M'+mnths+'::'
out_fname_base='fv3lam_merra2_'

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'


# Main part (will be called at the end) =================== CHJ =====
def main():
# ========================================================= CHJ =====
    global merra2,lat,lon,extent,c_lon

    print(' ===== INPUT: '+fnm_input+' ================================')
    # open the data file
    fname=os.path.join(dnm_data,fnm_input)
    try: merra2=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(merra2)

    lon1d=np.ma.masked_invalid(merra2['lon'].data)
    lat1d=np.ma.masked_invalid(merra2['lat'].data)

    lon,lat=np.meshgrid(lon1d,lat1d,sparse=False)    

    print(lat.shape)
    print(lon.shape)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon1d)
    lon_max=np.max(lon1d)
    lat_min=np.min(lat1d)
    lat_max=np.max(lat1d)

   # extent=[lon_min-5,lon_max+5,lat_min-5,lat_max+3]
    extent=[lon_min,lon_max,lat_min,lat_max]
    c_lon=np.mean(extent[:2])


    # Variables
    for svar in vars_atm:
        merra2_plot(svar)


# MERRA2 plot ================================================= CHJ =====
def merra2_plot(svar):
# ============================================================= CHJ =====

    print(' ===== '+svar+' ===== MERRA2 aerosol climatology =============')
    # Extract data array
    sfld=np.ma.masked_invalid(merra2[svar].data)

    ndim_svar=sfld.ndim
    lvlm=lvl-1
    lvls=format(lvl,'03d')

    if ndim_svar==3:
        (nts,nys,nxs)=sfld.shape
        print(' time+2D: nts=',nts,' nys=',nys,' nxs=',nxs)
        sfld2d=np.squeeze(sfld,axis=0)
        out_title_fld=out_title_base+svar
        out_fname=out_fname_base+svar
    elif ndim_svar==4:
        (nts,nls,nys,nxs)=sfld.shape
        print(' time+3D: nts=',nts,' nls=',nls,' nys=',nys,' nxs=',nxs)
        sfld2d=sfld[0,lvlm,:,:]
        out_title_fld=out_title_base+svar+'(L='+lvls+')'
        out_fname=out_fname_base+svar+'_L'+lvls
    else:
        sys.exit('ERROR: wrong dimension !!!')

    nm_svar=svar
    cs_cmap='gist_ncar'
    lb_ext='both'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=2
    cmap_range='real'

    print(' MERR2 field=',nm_svar)

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
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title_fld,fontsize=9)
    cs=ax.pcolormesh(lon,lat,sfld2d,cmap=cs_cmap,rasterized=True,
        vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(nm_svar,fontsize=8)

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

