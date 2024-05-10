###################################################################### CHJ #####
## Name         : plot_landda_ioda.py
## Language     : Python 3.7
## Usage        : Plot input snow depth ioda data
## Input files  : ghcn_snwd_ioda.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2024/04/18: Chan-Hoo Jeon : Preliminary version
###################################################################### CHJ #####

import os, sys
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
dnm_data="/scratch2/NAGAPE/epic/Chan-hoo.Jeon/ghcn_ioda/obs/"

# input file date
#date_input='20230501'
date_input='20000103'

# input file name
fnm_input='fake_ghcn_snwd_ioda_'+date_input+'.nc'

# basic forms of title and file name: base+static field name
out_title_base='LAND-DA::Snow depth (ioda)::'+date_input+'::'
out_fname_base='landda_snwd_ioda_'+date_input+'_'

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'


# Main part (will be called at the end) =================== CHJ =====
def main():
# ========================================================= CHJ =====
    global mdat,lat,lon,extent,c_lon

    print(' ===== INPUT: '+fnm_input+' ================================')
    # open the data file
    fpath=os.path.join(dnm_data,fnm_input)
    try: mdat=nc.Dataset(fpath)
    except: raise Exception('Could NOT find the file',fpath)
    print(mdat)
    print(mdat.groups['MetaData'])
    print(mdat.groups['ObsError'])
    print(mdat.groups['ObsValue'])
    print(mdat.groups['PreQC'])

    longitude=mdat.groups['MetaData'].variables['longitude'][:]
    latitude=mdat.groups['MetaData'].variables['latitude'][:]
    datetime=mdat.groups['MetaData'].variables['dateTime'][:]
#    stationElevation=mdat.groups['MetaData'].variables['stationElevation'][:]
    height=mdat.groups['MetaData'].variables['height'][:]
    stationID=mdat.groups['MetaData'].variables['stationIdentification'][:]

#    lon_max=np.max(lon)
    # Longitude 0:360 => -180:180
#    if lon_max>180:
#        lon=(lon+180)%360-180

    lon=longitude
    lat=latitude
#    lon,lat=np.meshgrid(longitude,latitude,sparse=False)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon)
    lon_max=np.max(lon)
    lat_min=np.min(lat)
    lat_max=np.max(lat)

    extent=[lon_min,lon_max,lat_min,lat_max]
    # for CONUS
    extent=[-125,-66,23,53]
    print(extent)

    c_lon=np.mean(extent[:2])
    #c_lon=-77.0369 # D.C.
    print(' c_lon=',c_lon)

    # Variables
#    vars_out=["ObsValue","ObsError","PreQC"]
    vars_out=["ObsValue"]
    for svar in vars_out:
        svar_plot(svar)


# MERRA2 plot ================================================= CHJ =====
def svar_plot(svar):
# ============================================================= CHJ =====

    print(' ===== '+svar+' ===== IODA total snow depth =============')
    # Extract data array
    sfld=mdat.groups[svar].variables['totalSnowDepth'][:]

    out_title_fld=out_title_base+svar
    out_fname=out_fname_base+svar

    cs_cmap='gist_ncar_r'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    scat_sz=1.5
    n_rnd=2
    cmap_range='fixed'

    print(' svar name=',svar)

    # Max and Min of the field
    fmax=np.max(sfld)
    fmin=np.min(sfld)
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
        cs_min=0
        cs_max=2000.0
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_max=',cs_max)
    print(' cs_min=',cs_min)
    print(' extent=',extent)

    # Plot field
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title_fld,fontsize=9)
#    cs=ax.pcolormesh(lon,lat,sfld,cmap=cs_cmap,rasterized=True,
#            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    cs=ax.scatter(lon,lat,transform=ccrs.PlateCarree(),c=sfld,cmap=cs_cmap,
                  vmin=cs_min,vmax=cs_max,s=scat_sz)
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(svar,fontsize=8)

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

