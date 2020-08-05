###################################################################### CHJ #####
## Name		: plot_fv3lam_fixam.py
## Language	: Python 3.7
## Usage	: Plot global static 'grb' files for fv3 regional modeling
## Input files  : [svar].grb
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/03/18: Chan-Hoo Jeon : Preliminary version
## V001: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
###################################################################### CHJ #####

import os, sys
import pygrib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
#import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable


# HPC machine ('hera','orion')
machine='hera'

print(' You are on', machine)

# Path to Natural Earth Data-set for background plot
if machine=='hera':
    cartopy.config['data_dir']='/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/NaturalEarth'
elif machine=='orion':
    cartopy.config['data_dir']='/home/chjeon/tools/NaturalEarth'
else:
    sys.exit('ERROR: path to Natural Earth Data is not set !!!')

plt.switch_backend('agg')

# Global variables ======================================== CHJ =====
# ..... Case-dependent input :: should be changed case-by-case .....
# ******
# INPUT
# ******
# Path to the directory where the input NetCDF files are located.
dnm_data="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/run_C96/"
#dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/regional_workflow/fix/fix_am/"

# Static atmospheric fields
#vars_fixam=["glacier","maxice","snoclim","soilmgldas","seaice","RTGSST"]
vars_fixam=["RTGSST"]

# *******
# OUTPUT
# *******
# Path to directory
if machine=='hera':
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig/"
elif machine=='orion':
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
else:
    sys.exit('ERROR: path to output directory is not set !!!')

# basic forms of title and file name: base+static field name
out_title_base='FIX_AM::'
out_fname_base='fv3_fixam_'

# Resolution of background natural earth data ('50m' or '110m')
back_res='110m'


# Main part (will be called at the end) ==================== CHJ =====
def main():
# ========================================================== CHJ =====

    for svar in vars_fixam:
        fixam_set(svar) 


# Static field set ========================================= CHJ =====
def fixam_set(svar):
# ========================================================== CHJ =====

    global cs_label,cs_cmap,lb_ext,cs_tick

    # file names 
    if svar=='glacier':
        fnm_in='global_glacier.2x2.grb'
    elif svar=='maxice':
        fnm_in='global_maxice.2x2.grb'
    elif svar=='snoclim':
        fnm_in='global_snoclim.1.875.grb'
    elif svar=='soilmgldas':
        fnm_in='global_soilmgldas.t1534.3072.1536.grb'
    elif svar=='seaice':
        fnm_in='CFSR.SEAICE.1982.2012.monthly.clim.grb'
    elif svar=='RTGSST':
        fnm_in='RTGSST.1982.2012.monthly.clim.grb'
    else:
        sys.exit('ERROR: wrong svar: '+svar)

    # open the data file
    fname=os.path.join(dnm_data,fnm_in)
    try: grbs=pygrib.open(fname)
    except: raise Exception('Could NOT find the file',fname)

    print(' ===== '+svar+' ===== fix_am data =======================')


#    for grb in grbs:
#        print(grb.name)
#        print(grb.typeOfLevel)
#        print(grb.level)
#        print(grb.validDate)
#        print(grb.analDate)
#        print(grb.Nx)
#        print(grb.Ny)
#        print(grb.shortName)
#        print(grb.missingValue)
    
    # default 
    cs_cmap='YlGnBu'
    lb_ext='neither'
    n_rnd=0
    
    if svar=='glacier':
        grb1=grbs.select(typeOfLevel='surface')[0]
        out_title_fld=out_title_base+svar
        out_fxam_fname=out_fname_base+svar
        cs_label=svar
        fxam=grb1.values
        lat,lon=grb1.latlons()
        fixam_plot(svar,fxam,lat,lon,out_title_fld,out_fxam_fname,n_rnd)
    elif svar=='maxice':
        grb1=grbs.select(typeOfLevel='surface')[0]
        print(grb1.name)
        out_title_fld=out_title_base+svar+'::Ice Cover (surface)'
        out_fxam_fname=out_fname_base+svar
        cs_label='Ice cover (1=ice, 0=no ice)'
        fxam=grb1.values
        lat,lon=grb1.latlons()
        fixam_plot(svar,fxam,lat,lon,out_title_fld,out_fxam_fname,n_rnd)
    elif svar=='snoclim':
        for ig in range(0,12):   # month
            grbm=grbs.select(name='Snow Fall water equivalent')[ig]
            smon=str(ig+1).zfill(2)
            print(grbm.name)
            print(grbm.analDate)
            fxam=grbm.values
            lat,lon=grbm.latlons()
            out_title_fld=out_title_base+svar+'::Snow Fall Water Equivalent::M'+smon
            out_fxam_fname=out_fname_base+svar+'_m'+smon
            cs_label='Snow fall water equivalent'
            cs_cmap='nipy_spectral_r'
            fixam_plot(svar,fxam,lat,lon,out_title_fld,out_fxam_fname,n_rnd)
    elif svar=='soilmgldas':
        ic=0
        for ig in range(0,12):   # month
            for il in [0,10,40,100]:   # level 
                grbm=grbs.select(typeOfLevel='depthBelowLandLayer')[ic]
                smon=str(ig+1).zfill(2)
                slvl=str(il).zfill(3)
                print('level=',grbm.level)
                print(grbm.analDate)
                fxam=grbm.values
                lat,lon=grbm.latlons()
                out_title_fld=out_title_base+svar+'::M'+smon+'::L'+slvl
                out_fxam_fname=out_fname_base+svar+'_m'+smon+'_L'+slvl
                cs_label=svar
                n_rnd=5
                cs_cmap='nipy_spectral_r'
                fixam_plot(svar,fxam,lat,lon,out_title_fld,out_fxam_fname,n_rnd)
                ic+=1
    elif svar=='seaice':
        for ig in range(0,12):   # month
            grbm=grbs.select(typeOfLevel='meanSea')[ig]
            smon=str(ig+1).zfill(2)
            print(grbm.name)
            print(grbm.typeOfLevel)
            print(grbm.analDate)
            fxam=grbm.values
            lat,lon=grbm.latlons()
            out_title_fld=out_title_base+svar+'::Ice Cover (Mean Sea)::M'+smon
            out_fxam_fname=out_fname_base+svar+'_m'+smon
            cs_label='Ice cover (1=ice, 0=no ice)'
            n_rnd=2
            fixam_plot(svar,fxam,lat,lon,out_title_fld,out_fxam_fname,n_rnd)
    elif svar=='RTGSST':
        for ig in range(0,12):   # month
            grbm=grbs.select(name='Temperature')[ig]
            smon=str(ig+1).zfill(2)
            print(grbm.name)
            print(grbm.typeOfLevel)
            print(grbm.analDate)
            fxam=grbm.values
            lat,lon=grbm.latlons()
            out_title_fld=out_title_base+svar+'::SST::M'+smon
            out_fxam_fname=out_fname_base+svar+'_m'+smon
            cs_label='Temperature'
            n_rnd=2
            cs_cmap='jet'
            fixam_plot(svar,fxam,lat,lon,out_title_fld,out_fxam_fname,n_rnd)
    else:
        sys.exit('ERROR: Wrong svar !!')



# Static field plot ========================================= CHJ =====
def fixam_plot(svar,fxam,lat,lon,out_title_fld,out_fxam_fname,n_rnd):
# ========================================================== CHJ =====

    
    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon)
    lon_max=np.max(lon)
    lat_min=np.min(lat)
    lat_max=np.max(lat)
    
   # extent=[lon_min-5,lon_max+5,lat_min-5,lat_max+3]
    extent=[lon_min,lon_max,lat_min,lat_max]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])
 
    # Max and Min of the field
    fmax=np.nanmax(fxam)
    fmin=np.nanmin(fxam)

    print('fld_min=',fmin)
    print('fld_max=',fmax)
    
    cs_min=round(fmin,n_rnd)
    cs_max=round(fmax,n_rnd)

    print('cs_min=',cs_min)
    print('cs_max=',cs_max)
 
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3

    # Plot field
    fig,ax1=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax1.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax1)
    ax1.set_title(out_title_fld,fontsize=9)
    cs=ax1.pcolormesh(lon,lat,fxam,cmap=cs_cmap,rasterized=True,
            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    # extend(pointed end): 'neither'|'both'|'min'|'max'  
    divider=make_axes_locatable(ax1)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(cs_label,fontsize=8)

    # Output figure
    out_file(out_fxam_fname) 

       
   

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
#    states=cfeature.NaturalEarthFeature('cultural','admin_1_states_provinces',
#                      back_res,edgecolor='black',facecolor='none',
#                      linewidth=fline_wd,linestyle=':',alpha=falpha)
    borders=cfeature.NaturalEarthFeature('cultural','admin_0_countries',
                      back_res,edgecolor='red',facecolor='none',
                      linewidth=fline_wd,alpha=falpha)

#    ax.add_feature(land)
    ax.add_feature(lakes)
#    ax.add_feature(states)
    ax.add_feature(borders)
    ax.add_feature(coastline)




# Output file ============================================= CHJ =====
def out_file(out_file):
# ========================================================= CHJ =====
    # Output figure
    plt.savefig(out_fig_dir+out_file+'.png',dpi=300,bbox_inches='tight')
    plt.close('all')



# Main call ================================================ CHJ =====
if __name__=='__main__':
    main()

