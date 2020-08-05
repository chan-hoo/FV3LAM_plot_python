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

# Path to Natural Earth Data-set for background plot
if machine=='hera':
    path_NE='/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/NaturalEarth'
elif machine=='orion':
    path_NE='/home/chjeon/tools/NaturalEarth'
else:
    sys.exit('ERROR: path to Natural Earth Data is not set !!!')

cartopy.config['data_dir']=path_NE
os.environ["CARTOPY_USER_BACKGROUNDS"]=path_NE+'/raster_files'

plt.switch_backend('agg')

# Global variables ======================================== CHJ =====
# ..... Case-dependent input :: should be changed case-by-case .....
# ******
# INPUT
# ******
# Path to the directory where the input NetCDF files are located.
dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/expt_dirs/test_community_hrrr25/2020061800/postprd/"

# Input file (BGRD3DXX.tmXX)
#fnm_in='fv3sar.tz.conus.natlev.f01.grib2'
fnm_in='BGRD3D_2017000000600'

# Domain name:
domain='HRRR'

# Grid resolution ('C96'/'C768'):
res='C401'

# Grid type ('ESG'/'GFDL')
gtype='ESG'

# GFDL grid-refinement ratio (for ESG grid, refine=0)
if gtype=='ESG':
    refine=0
elif gtype=='GFDL':
    refine=3

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
vars_grb2=["V wind in boundary layer"]


# Total number of vertical layers in the case
nlvl=64

# Layer number to be plotted
ilvl=1

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

# basic forms of title and file name
if gtype=='ESG':
    out_title_base='FV3-LAM::'+domain+'(ESG)::'+res+'::'
    out_fname_base='fv3_out_lev_'+domain+'_esg_'+res+'_'
elif gtype=='GFDL':
    out_title_base='FV3-LAM::'+domain+'(GFDL)::'+res+'(x'+str(refine)+')'+'::'
    out_fname_base='fv3_out_lev_'+domain+'_gfdl_'+res+'_'

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'
# high-resolution background image ('on', 'off')
back_img='on'



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
        cs_min=-10.0
        cs_max=10.0
    elif cmap_range=='designed':
        cs_min=5
        cs_max=None
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_min=',cs_min)
    print(' cs_max=',cs_max)
 
    # Plot field
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title,fontsize=9)
    cs=ax.pcolormesh(lon,lat,sval,cmap=cs_cmap,rasterized=True,
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

