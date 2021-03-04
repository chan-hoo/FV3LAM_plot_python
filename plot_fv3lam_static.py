###################################################################### CHJ #####
## Name		: plot_fv3lam_static.py
## Language	: Python 3.7
## Usage	: Plot static fields for fv3 regional modeling
## Input files  : CXX_[svar].tile7.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/02/28: Chan-Hoo Jeon : Preliminary version
## V001: 2020/03/02: Chan-Hoo Jeon : Plot global (original) data
## V002: 2020/03/05: Chan-Hoo Jeon : Add high-resolution earth data (50m)
## V003: 2020/03/24: Chan-Hoo Jeon : Separate the path to orography
## V004: 2020/03/31: Chan-Hoo Jeon : Modify plot options for temperature
## V005: 2020/04/07: Chan-Hoo Jeon : Add refine ratio to output titles
## V006: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
## V007: 2021/03/04: Chan-Hoo Jeon : Simplify the script
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
# prt_data: global surface climatology data
# out_fig_dir: directory where the output files are created
# mfdt_kwargs: mfdataset argument

if machine=='hera':
    cartopy.config['data_dir']='/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/NaturalEarth'
    prt_data="/scratch1/NCEPDEV/global/glopara/fix/fix_sfc_climo/"
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig/"
    mfdt_kwargs={'parallel':False}
elif machine=='orion':
    cartopy.config['data_dir']='/home/chjeon/tools/NaturalEarth'
    prt_data="/work/noaa/global/glopara/fix/fix_sfc_climo/"
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
    mfdt_kwargs={'parallel':False,'combine':'by_coords'}
else:
    sys.exit('ERROR: Required input Data are NOT set !!!')

plt.switch_backend('agg')

# Case-dependent parameters  ========================================= CHJ =====
# Path to the directory where the input NetCDF files are located.
dnm_data="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/expt_dirs/test_community/sfc_climo/"

# Path to the orography file
oro_path="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/expt_dirs/test_community/2020122700/INPUT/"

# Basic form of Input NetCDF file names: res+static+basic form
fnm_in_base='.tile7.halo4.nc' # for halo4

# Grid resolution ('C96'/'C768'):
res='C403'

# Orography data file (w/ halo0 or halo4) for coordinates
fnm_in_orog='oro_data.tile7.halo4.nc'  # for halo4

# Static fields
# comment out below in case of no time-limit
#vars_static=["facsf","maximum_snow_albedo","slope_type","soil_type",
#             "substrate_temperature","vegetation_type",
#             "vegetation_greenness","snowfree_albedo"]
# use below two separately in case of time-limit < 30 min.
#vars_static=["facsf","maximum_snow_albedo","slope_type","soil_type",
#             "substrate_temperature","vegetation_type",
#             "vegetation_greenness"]
# for snowfree_albedo, split vars_subs in the funtion of static_plot into two.
#vars_static=["snowfree_albedo"]

vars_static=["substrate_temperature","vegetation_greenness"]

# basic forms of title and file name: base+static field name
out_title_base='FV3LAM::Static (FIX)::'
out_fname_base='fv3lam_static_'

out_title_base_g='FV3LAM::FIX::Original::'

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'


# Main part (will be called at the end) =================== CHJ =====
def main():
# ========================================================= CHJ =====

   # Open '_oro_data' file for coordinates
    oro_data()

    # Static fields
    for svar in vars_static:
        static_plot(svar)

 
# Orography data for (lon,lat) ============================= CHJ =====
def oro_data():
# ========================================================== CHJ =====
    global glon,glat,extent,c_lon,c_lat
     # open the orography file
    fname=os.path.join(oro_path,fnm_in_orog)
    try: oro=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)

    print(' ===== OROGRAPHY ================================')
    print(oro)

    # Extract longitudes, and latitudes
    glon=np.ma.masked_invalid(oro['geolon'].data)
    glat=np.ma.masked_invalid(oro['geolat'].data)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(glon)
    lon_max=np.max(glon)
    lat_min=np.min(glat)
    lat_max=np.max(glat)

    # Plot extent
    extent=[lon_min-5,lon_max+5,lat_min-5,lat_max+3]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])
  


# Static field plot ======================================== CHJ =====
def static_plot(svar):
# ========================================================== CHJ =====
    
    fnm_input_var=res+'.'+svar+fnm_in_base

    # open the data file
    fname=os.path.join(dnm_data,fnm_input_var)
    try: sttc=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)

    print(' ===== '+svar+' ===== regional data =======================')
    print(sttc)

    # original files ========================= CHJ =====
    if svar=='facsf':
        svar_g='facsf.1.0.nc'
    elif svar=='maximum_snow_albedo':
        svar_g='maximum_snow_albedo.0.05.nc'
    elif svar=='slope_type':
        svar_g='slope_type.1.0.nc'
    elif svar=='soil_type':
        svar_g='soil_type.statsgo.0.03.nc'
    elif svar=='substrate_temperature':
        svar_g='substrate_temperature.2.6x1.5.nc'
    elif svar=='vegetation_type':
        svar_g='vegetation_type.igbp.0.05.nc'
    elif svar=='vegetation_greenness':
        svar_g='vegetation_greenness.0.144.nc'
    elif svar=='snowfree_albedo':
        svar_g='snowfree_albedo.4comp.0.05.nc'
    else:
        sys.exit('ERROR: wrong svar: '+svar)
    # ========================================= CHJ =====

    # open the global data file
    fnm_g=os.path.join(prt_data,svar_g)
    try: sttc_g=xr.open_mfdataset(fnm_g,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)

    print(' ===== '+svar_g+' ===== original data =====================')
    print(sttc_g)

    # Extract data array
    if svar=="snowfree_albedo":
# in case of no time-limit
        svar_sub=["visible_black_sky_albedo","visible_white_sky_albedo",
                  "near_IR_black_sky_albedo","near_IR_white_sky_albedo"]
# use below one by one in case of time-limit < 30min.
#        svar_sub=["visible_black_sky_albedo","visible_white_sky_albedo"]
#        svar_sub=["near_IR_black_sky_albedo","near_IR_white_sky_albedo"]

        for v_sub in svar_sub:
            sfld=np.ma.masked_invalid(sttc[v_sub].data)
            sfld_g=np.ma.masked_invalid(sttc_g[v_sub].data)

            lat1d=np.ma.masked_invalid(sttc_g['lat'].data)
            lon1d=np.ma.masked_invalid(sttc_g['lon'].data)
            glon_g,glat_g=np.meshgrid(lon1d,lat1d)

            static_plot_sub(v_sub,sfld,sfld_g,glon_g,glat_g)
    else:
        sfld=np.ma.masked_invalid(sttc[svar].data)
        sfld_g=np.ma.masked_invalid(sttc_g[svar].data)

        lat1d=np.ma.masked_invalid(sttc_g['lat'].data)
        lon1d=np.ma.masked_invalid(sttc_g['lon'].data)
        glon_g,glat_g=np.meshgrid(lon1d,lat1d)

        static_plot_sub(svar,sfld,sfld_g,glon_g,glat_g)

      

# Static field plot ======================================== CHJ =====
def static_plot_sub(svar,sfld,sfld_g,glon_g,glat_g):
# ========================================================== CHJ =====
    print('Static field=',svar)
    (nft,nfr,nfc)=sfld.shape
    (ngr,ngc)=glon.shape

    if nfr!=ngr or nfc!=ngc:
        sys.exit('ERROR:: mismatched array size ::'+svar)

    # Max and Min of the field
    sfld_max=np.max(sfld)
    sfld_min=np.min(sfld)
    print('fld_max=',sfld_max)
    print('flx_min=',sfld_min)

# extend(pointed end): 'neither'|'both'|'min'|'max'  
    if svar=='facsf' or svar=='maximum_snow_albedo' or svar=='vegetation_greenness':
        svar_cmap='jet'
        svar_vmin=0
        svar_vmax=1
        lb_extend='neither'
        svar_tick=np.arange(0.0,1.2,0.2)
    elif svar=='slope_type':
        svar_cmap=plt.cm.get_cmap('Paired',9)
        svar_vmin=0.5
        svar_vmax=9.5
        lb_extend='neither'
        svar_tick=np.arange(1,10)
    elif svar=='soil_type':
        svar_cmap=plt.cm.get_cmap('tab20',16)
        svar_vmin=0.5
        svar_vmax=16.5
        lb_extend='neither'
        svar_tick=np.arange(1,17)
    elif svar=='vegetation_type':
        svar_cmap=plt.cm.get_cmap('tab20',19)
        svar_vmin=0.5
        svar_vmax=19.5
        lb_extend='neither'
        svar_tick=np.arange(1,20,2)
    elif svar=='substrate_temperature':
        svar_cmap='jet'
        cs_max=round(sfld_max,1)
        cs_min=round(sfld_min,1)
        cs_del=round((cs_max-cs_min-0.1)/5,2)
        print(' cs_max=',cs_max)
        print(' cs_min=',cs_min)
        svar_vmin=cs_min
        svar_vmax=cs_max
        lb_extend='neither'
        svar_tick=np.arange(cs_min,cs_max,cs_del)
    elif svar=='visible_black_sky_albedo':
        svar_cmap='YlGnBu'
        svar_vmin=0.0
        svar_vmax=0.7
        lb_extend='max'
        svar_tick=np.arange(0.0,0.8,0.1)
    elif svar=='visible_white_sky_albedo':
        svar_cmap='YlOrRd'
        svar_vmin=0.0
        svar_vmax=0.7
        lb_extend='max'
        svar_tick=np.arange(0.0,0.8,0.1)
    elif svar=='near_IR_black_sky_albedo':
        svar_cmap='YlGnBu'
        svar_vmin=0.0
        svar_vmax=0.6
        lb_extend='max'
        svar_tick=np.arange(0.0,0.7,0.1)
    elif svar=='near_IR_white_sky_albedo':
        svar_cmap='YlOrRd'
        svar_vmin=0.0
        svar_vmax=0.6
        lb_extend='max'
        svar_tick=np.arange(0.0,0.7,0.1)
    else:
        sys.exit('ERROR:: wrong field !!!')


    for imn in range(nft):
        if nft==1:
        # Remove the time dimension in the array in case of 1
            sfld2d=np.squeeze(sfld,axis=0)
            sfld2d_g=np.squeeze(sfld_g,axis=0)
            out_title_fld=out_title_base+svar
            out_title_fld_g=out_title_base_g+svar
            out_sfld_fname=out_fname_base+svar
        else:
            print('month=',imn+1)
            sfld2d=sfld[imn,:,:]
            sfld2d_g=sfld_g[imn,:,:]
            mon_nm='M'+format(imn+1,'02d')
            out_title_fld=out_title_base+svar+'::'+mon_nm
            out_title_fld_g=out_title_base_g+svar+'::'+mon_nm
            out_sfld_fname=out_fname_base+svar+'_'+mon_nm

        # Plot field
        fig,(ax1,ax2)=plt.subplots(2,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
        # Plot 1: original global data
        ax1.set_extent(extent, ccrs.PlateCarree())
        # Call background plot
        back_plot(ax1)
        ax1.set_title(out_title_fld_g,fontsize=9)
        cs=ax1.pcolormesh(glon_g,glat_g,sfld2d_g,cmap=svar_cmap,
              rasterized=True,vmin=svar_vmin,vmax=svar_vmax,
              transform=ccrs.PlateCarree())
        # extend(pointed end): 'neither'|'both'|'min'|'max'  
        divider=make_axes_locatable(ax1)
        ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_extend,ticks=svar_tick)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label(svar,fontsize=8)

        # Plot 2: regional data
        ax2.set_extent(extent, ccrs.PlateCarree())
        # Call background plot
        back_plot(ax2)
        ax2.set_title(out_title_fld,fontsize=9)
        cs=ax2.pcolormesh(glon,glat,sfld2d,cmap=svar_cmap,rasterized=True,
              vmin=svar_vmin,vmax=svar_vmax,transform=ccrs.PlateCarree())
        # extend(pointed end): 'neither'|'both'|'min'|'max'  
        divider=make_axes_locatable(ax2)
        ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_extend,ticks=svar_tick)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label(svar,fontsize=8)


        # Output figure
        out_file(out_sfld_fname) 

       
   

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
def out_file(out_file):
# ========================================================= CHJ =====
    # Output figure
    plt.savefig(out_fig_dir+out_file+'.png',dpi=300,bbox_inches='tight')
    plt.close('all')



# Main call ================================================ CHJ =====
if __name__=='__main__':
    main()

