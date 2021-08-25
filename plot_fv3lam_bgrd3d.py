###################################################################### CHJ #####
## Name		: plot_fv3lam_bgrd3d.py
## Language	: Python 3.7
## Usage	: Plot output 'BGRD3D.tmXX(*.lev.grib2)' files for fv3 regional modeling
## Input files  : BGRD3DXX.tmXX(*.lev.grib2)
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/05/12: Chan-Hoo Jeon : Preliminary version
## V001: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
## V002: 2020/07/20: Chan-Hoo Jeon : Add new cmap and background image option
## V003: 2021/03/05: Chan-Hoo Jeon : Simplify the script
## V004: 2021/06/24: Chan-Hoo Jeon : Add a projection for RRFS_NA domain
## V005: 2021/08/25; Chan-Hoo Jeon : Add vars for RRFS-CMAQ
###################################################################### CHJ #####

import os, sys
import pygrib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
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

if machine=='hera':
    cartopy.config['data_dir']='/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/NaturalEarth'
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig/"
elif machine=='orion':
    cartopy.config['data_dir']='/home/chjeon/tools/NaturalEarth'
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
else:
    sys.exit('ERROR: Required input data are NOT set !!!')

plt.switch_backend('agg')

# Case-dependent input =============================================== CHJ =====

# Path to the directory where the input NetCDF files are located.
dnm_data="/scratch2/NCEPDEV/stmp3/Chan-hoo.Jeon/expt_dirs/test_cmaq_conus13_1day/2019080112/postprd/"
#dnm_data="/scratch2/NCEPDEV/naqfc/Chan-hoo.Jeon/test_wcoss/"

# Domain name
domain_nm='RRFS_CONUS_13km'

# Input file hour number
fnm_hr='f003'

# Input file (BGRD3D / NATLEV)
#fnm_in='rrfs.t00z.bgrd3d'+fnm_hr+'.tm00.grib2'
fnm_in='rrfs.t12z.natlev'+fnm_hr+'.tm00.grib2'

# Output fields
##### hybrid levels #####
#vars_grb2=["Temperature on model surface"]
#vars_grb2=["Specific humidity on model surface"]
#vars_grb2=["Cloud mixing ratio on model surface"]
#vars_grb2=["Rain mixing ratio on model surface"]
#vars_grb2=["Snow mixing ratio on model surface"]
#vars_grb2=["Turbulent kinetic energy on model surface"]
#vars_grb2=["U component wind on model surface"]
#vars_grb2=["V component wind on model surface"]

##### pressure levels: 250, 500, 700, 850 #####
#vars_grb2=["Temperature on pressure surface"]
#vars_grb2=["Specific humidity on pressure surface"]
#vars_grb2=["U component wind on pressure surface"]
#vars_grb2=["V component wind on pressure surface"]

##### flight levels: 305, 457, 610, 914, 1524, 1829, 2134, 2743, 3658, 4572 #####
#vars_grb2=["U wind at flight levels"]
#vars_grb2=["V wind at flight levels"]

##### Boundary layers: 3000, 6000, 9000, 12000, 15000, 18000 #####
#vars_grb2=["Temperature in boundary layer"]
#vars_grb2=["Dew point temperature in boundary layer"]
#vars_grb2=["Specific humidity in boundary layer"]
#vars_grb2=["U wind in boundary layer"]
#vars_grb2=["V wind in boundary layer"]

##### RRFS-CMAQ #####
vars_grb2=["pmtf","ozcon"]


# Total number of vertical layers in the case
nlvl=64

# Layer number to be plotted
ilvl=1

# basic forms of title and file name
out_title_base='FV3LAM::BGRD3D::'+domain_nm+'::'+fnm_hr+'::'
out_fname_base='fv3lam_out_bgrd3d_'+domain_nm+'_'+fnm_hr+'_'

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'
# high-resolution background image ('on', 'off')
back_img='off'


# Main part (will be called at the end) ==================== CHJ =====
def main():
# ========================================================== CHJ =====
    global grbs

    # open the data file
    fname=os.path.join(dnm_data,fnm_in)
    try: grbs=pygrib.open(fname)
    except: raise Exception('Could NOT find the file',fname)

    for grb in grbs:
        print(grb.name)
        print(grb.typeOfLevel)
        print(grb.level)
#        print(grb.validDate)
#        print(grb.analDate)
#        print(grb.Nx)
#        print(grb.Ny)
        print(grb.shortName)
#        print(grb.missingValue)

    for svar in vars_grb2:
        plot_grb2(svar) 



# Output field plot ======================================== CHJ =====
def plot_grb2(svar):
# ========================================================== CHJ =====

    ilvlm=ilvl-1

    print(' ===== '+svar+' ===== lev ===============================')

    nm_svar=svar
    cs_cmap='jet'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=2
    cmap_range='round'
    df_min=-10
    df_max=10

    if svar=="Temperature on model surface":
        grbv=grbs.select(name="Temperature",typeOfLevel="hybrid")[ilvlm]
    elif svar=="Temperature on pressure surface":     
        grbv=grbs.select(name="Temperature",typeOfLevel="isobaricInhPa")[ilvlm]
    elif svar=="Temperature in boundary layer":       
        grbv=grbs.select(name="Temperature",typeOfLevel="pressureFromGroundLayer")[ilvlm]
    elif svar=="Dew point temperature in boundary layer":       
        grbv=grbs.select(name="Dew point temperature",typeOfLevel="pressureFromGroundLayer")[ilvlm]
    elif svar=="Specific humidity on model surface":
        grbv=grbs.select(name="Specific humidity",typeOfLevel="hybrid")[ilvlm]
        n_rnd=3 
    elif svar=="Specific humidity on pressure surface":       
        grbv=grbs.select(name="Specific humidity",typeOfLevel="isobaricInhPa")[ilvlm]
        n_rnd=5
    elif svar=="Specific humidity in boundary layer":       
        grbv=grbs.select(name="Specific humidity",typeOfLevel="pressureFromGroundLayer")[ilvlm]
        n_rnd=3
    elif svar=="Cloud mixing ratio on model surface":
        grbv=grbs.select(name="Cloud mixing ratio",typeOfLevel="hybrid")[ilvlm]
        n_rnd=5
    elif svar=="Rain mixing ratio on model surface":
        grbv=grbs.select(name="Rain mixing ratio",typeOfLevel="hybrid")[ilvlm]
        n_rnd=5
    elif svar=="Snow mixing ratio on model surface":
        grbv=grbs.select(name="Snow mixing ratio",typeOfLevel="hybrid")[ilvlm]
        n_rnd=5
    elif svar=="Turbulent kinetic energy on model surface":
        grbv=grbs.select(name="Turbulent kinetic energy",typeOfLevel="hybrid")[ilvlm]
        n_rnd=5
    elif svar=="U component wind on model surface":
        grbv=grbs.select(name="U component of wind",typeOfLevel="hybrid")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="U component wind on pressure surface":  
        grbv=grbs.select(name="U component of wind",typeOfLevel="isobaricInhPa")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="U wind at flight levels":    
        grbv=grbs.select(name="U component of wind",typeOfLevel="heightAboveSea")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="U wind in boundary layer":       
        grbv=grbs.select(name="U component of wind",typeOfLevel="pressureFromGroundLayer")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="V component wind on model surface":
        grbv=grbs.select(name="V component of wind",typeOfLevel="hybrid")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="V component wind on pressure surface":    
        grbv=grbs.select(name="V component of wind",typeOfLevel="isobaricInhPa")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="V wind at flight levels":    
        grbv=grbs.select(name="V component of wind",typeOfLevel="heightAboveSea")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="V wind in boundary layer":       
        grbv=grbs.select(name="V component of wind",typeOfLevel="pressureFromGroundLayer")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="ozcon":
        grbv=grbs.select(name="Ozone Concentration (PPB)",shortName="ozcon")[ilvlm]
    elif svar=="pmtf":
        grbv=grbs.select(name="Particulate matter (fine)",shortName="pmtf")[ilvlm]
        cs_cmap="gist_ncar_r"        
    else:
        sys.exit('ERROR: Wrong svar or Not set up yet !!! ::'+svar)

    sval=grbv.values
    sname=grbv.name
    if sname=='Temperature' or sname=='Dew point temperature':
        # Kelvin to Fahrenheit
        sval=(sval-273.15)*1.8+32.0
        nm_svar=svar+' ('+chr(176)+'F)'
    print(' Name=',sname)   
    stnm=grbv.shortName
    print(' ShortName=',stnm)
    tlxtr=str(grbv.level)
    print(' Level=',tlxtr)
    nmxtr=grbv.typeOfLevel
    print(' Type of level=',nmxtr)
    lat,lon=grbv.latlons()
    slvl=format(ilvlm+1,'03d')

    out_title=out_title_base+svar+'::L='+tlxtr
    out_fname=out_fname_base+stnm+'_'+nmxtr+'_L'+slvl


    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon)
    lon_max=np.max(lon)
    lat_min=np.min(lat)
    lat_max=np.max(lat)

    print(' lon_min=',lon_min,', lon_max=',lon_max)
    print(' lat_min=',lat_min,', lat_max=',lat_max)

    esp=1    
    extent=[lon_min-esp,lon_max+esp,lat_min-esp,lat_max+esp]
#    extent=[lon_min,lon_max,lat_min,lat_max]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])
 
    # Max and Min of the field
    fmax=np.nanmax(sval)
    fmin=np.nanmin(sval)

    print(' fld_min=',fmin)
    print(' fld_max=',fmax)

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
        cs_min=df_min
        cs_max=df_max
    elif cmap_range=='designed':
        cs_min=df_min
        cs_max=None
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_min=',cs_min)
    print(' cs_max=',cs_max)
 
    # Plot field
    if domain_nm[:7]=='RRFS_NA':
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(
                            central_longitude=-107,central_latitude=53)))
    else:
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
        ax.set_extent(extent, ccrs.PlateCarree())

    back_plot(ax)
    ax.set_title(out_title,fontsize=9)
    cs=ax.pcolormesh(lon,lat,sval,cmap=cs_cmap,rasterized=True,
            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(nm_svar,fontsize=6)

    # Output figure
    ndpi=300
    out_file(out_fname,ndpi) 



# new colormap option ====================================== CHJ =====
def new_cmap():
# ========================================================== CHJ =====
    c_lvls=np.linspace(5,70,14)
    c_list=['turquoise','dodgerblue','mediumblue','lime','limegreen','green','#EEEE00','#EEC900','darkorange','red','firebrick','darkred','fuchsia']
    new_cmap=colors.ListedColormap(c_list)
    new_norm=colors.BoundaryNorm(c_lvls,new_cmap.N)

    return new_cmap,new_norm



# Background plot ========================================== CHJ =====
def back_plot(ax):
# ========================================================== CHJ =====
    fline_wd=0.5  # line width
    falpha=0.3 # transparency

    # natural_earth
    land=cfeature.NaturalEarthFeature('physical','land',back_res,
                      edgecolor='face',facecolor=cfeature.COLORS['land'],
                      alpha=falpha)
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

    # high-resoultion background image
    if back_img=='on':
        ax.background_img(name='NE', resolution='high')

    # On/off features
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

