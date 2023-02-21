###################################################################### CHJ #####
## Name		: plot_rrfscmaq_emission.py
## Language	: Python 3.7
## Usage	: Plot an input emission file for rrfs_cmaq
## Input files  : GBBEPx_CRES.emissions.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2021/06/28: Chan-Hoo Jeon : Preliminary version
## V001: 2021/07/01: Chan-Hoo Jeon : Change to scatter plot
## V002: 2021/07/02: Chan-Hoo Jeon : Add grid_spec
## V003: 2021/07/05: Chan-Hoo Jeon : lon/lat check
## V004: 2023/02/21: Chan-Hoo Jeon : Add AQM_NA_13km
###################################################################### CHJ #####

import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from netCDF4 import Dataset
from scipy.io import netcdf
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
    path_fix="/scratch2/NCEPDEV/naqfc/RRFS_CMAQ/nexus/fix/"
elif machine=='orion':
    cartopy.config['data_dir']='/home/chjeon/tools/NaturalEarth'
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
    mfdt_kwargs={'parallel':False,'combine':'by_coords'}
else:
    sys.exit('ERROR: Required input data are NOT set !!!')

plt.switch_backend('agg')

# Case-dependent input =============================================== CHJ =====

# Domain name
domain_nm='AQM_NA_13km'

dnm_in="/scratch2/NCEPDEV/naqfc/RRFS_CMAQ/RAVE_fire/"
fnm_input='Hourly_Emissions_regrid_NA_13km_20230217_t06z_h72.nc'
grid_spec='grid_spec_793.nc'

# Variables
#vars_data=["CO2","CO","SO2","OC","BC","PM2.5","NOx","NH3","MeanFRP"]
vars_data=["MeanFRP"]

# Select time step
ihr=1

# basic forms of title and file name
out_title_base='AQM::Emission::'+domain_nm+'::'
out_fname_base='rrfscmaq_emission_'+domain_nm+'_'

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'



# Main part (will be called at the end) ======================= CHJ =====
def main():
# ============================================================= CHJ =====
    global ds,lon_o,lat_o
    global extent,c_lon,c_lat

    print(' ===== grid_spec =========================================')
    fname=os.path.join(path_fix,grid_spec)
    try: grd=Dataset(fname,'r')
    except: raise Exception('Could NOT find the file',fname)
    print(grd)

    lon_o=grd.variables['grid_lont'][:]
    lat_o=grd.variables['grid_latt'][:]
    lonc_o=grd.variables['grid_lon'][:]
    latc_o=grd.variables['grid_lat'][:]
   

    print(' ===== Input data ========================================')
    # open the data file
    fname=os.path.join(dnm_in,fnm_input)
    # Check the variables in netcdf file
    try: ds=Dataset(fname,'r')
    except: raise Exception('Could NOT find the file',fname)
    print(ds)

# =============================================================================
    
#    lon_max=np.max(lon)
    # Longitude 0:360 => -180:180
#    if lon_max>180:
#        lon=(lon+180)%360-180

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon_o)
    lon_max=np.max(lon_o)
    lat_min=np.min(lat_o)
    lat_max=np.max(lat_o)

    print(' lon_min=',lon_min,', lon_max=',lon_max)
    print(' lat_min=',lat_min,', lat_max=',lat_max)

    # Plot extent
    esp=1
    extent=[lon_min-esp,lon_max+esp,lat_min-esp,lat_max+esp]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])

    # Variables
    for svar in vars_data:
        data_plot(svar)



# ===== plot ================================================== CHJ =====
def data_plot(svar):
# ============================================================= CHJ =====

    print(' ===== '+svar+' ===== data ===============================')
    # Extract data array
    sfld=ds.variables[svar][:]

    ndim_svar=sfld.ndim

    if ndim_svar==2:
        (nys,nxs)=sfld.shape
        print(' 2D: nys=',nys,' nxs=',nxs)
        sfld2d=sfld
    elif ndim_svar==3:
        (nts,nys,nxs)=sfld.shape
        print(' time+2D: nts=',nts,' nys=',nys,' nxs=',nxs)
        sfld2d=sfld[ihr,:,:]
        sfld2d=np.squeeze(sfld2d)

    print(sfld.shape)
    print(sfld2d.shape)

    # Check if dimensions of grid and vars are matched
    if sfld2d.shape==lon_o.shape:
        print('Data and grid are matched !!!')
    else:
        sys.exit('ERROR: Size of data does NOT match with that of grid  !!!')

    # extract non-zero cells
    lon_pts=lon_o[sfld2d>0]
    lat_pts=lat_o[sfld2d>0]
    sfld_pts=sfld2d[sfld2d>0]

    print(lon_pts.shape)
    print(lat_pts.shape)
    print(sfld_pts.shape)

    out_title_fld=out_title_base+svar
    out_fname=out_fname_base+svar

    nm_svar=svar
    cs_cmap='jet'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    scat_sz=2

    cmap_range='round'
    cmap_fix_min=0.0
    cmap_fix_max=10.0

    if svar=="CO2":
#        n_rnd=8
        cmap_range='fixed'
        cmap_fix_min=0.0
        cmap_fix_max=3e-7
    elif svar=="CO":
#        n_rnd=9
        cmap_range='fixed'
        cmap_fix_min=0.0
        cmap_fix_max=1e-8
    elif svar=="SO2":
#        n_rnd=10
        cmap_range='fixed'
        cmap_fix_min=0.0
        cmap_fix_max=3e-9
    elif svar=="OC":
#        n_rnd=9
        cmap_range='fixed'
        cmap_fix_min=0.0
        cmap_fix_max=1e-8
    elif svar=="BC":
#        n_rnd=11
        cmap_range='fixed'
        cmap_fix_min=0.0
        cmap_fix_max=3e-10
    elif svar=="PM2.5":
#        n_rnd=9
        cmap_range='fixed'
        cmap_fix_min=0.0
        cmap_fix_max=6e-9
    elif svar=="NOx":
#        n_rnd=10
        cmap_range='fixed'
        cmap_fix_min=0.0
        cmap_fix_max=1.5e-9
    elif svar=="NH3":
#        n_rnd=10
        cmap_range='fixed'
        cmap_fix_min=0.0
        cmap_fix_max=1.25e-9
    elif svar=="MeanFRP":
#        n_rnd=2
        cmap_range='fixed'
        cmap_fix_min=0.0
        cmap_fix_max=50.0
    else:
        n_rnd=7


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
        cs_min=cmap_fix_min
        cs_max=cmap_fix_max
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
    cs=ax.scatter(lon_pts,lat_pts,transform=ccrs.PlateCarree(),c=sfld_pts,cmap=cs_cmap,
                  vmin=cs_min,vmax=cs_max,s=scat_sz)
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(nm_svar,fontsize=6)

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

