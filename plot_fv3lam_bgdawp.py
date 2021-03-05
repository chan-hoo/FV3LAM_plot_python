###################################################################### CHJ #####
## Name		: plot_fv3lam_bgdawp.py
## Language	: Python 3.7
## Usage	: Plot output 'BGDAWP(prs.grib2)' files for fv3 regional modeling
## Input files  : BGDAWPXX.tmXX(*.prs.grib2)
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/05/18: Chan-Hoo Jeon : Preliminary version
## V001: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
## V002: 2020/07/17: Chan-Hoo Jeon : Add new cmap and background image
## V003: 2021/03/05: Chan-Hoo Jeon : Simplify the script
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
# Path to the directory where the input NetCDF files are located.
dnm_data="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/expt_dirs/test_community/2020122700/postprd/"

# Input file hour number
fnm_hr='f003'

# Input file (BGDAWPXX.tmXX)
fnm_in='rrfs.t00z.bgdawp'+fnm_hr+'.tm00.grib2'

# Output fields
##### pressure levels: 
##### 2,5,7,10,20,30,50,70,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,
##### 500,525,550,575,600,625,650,675,700,725,750,775,800,825,850,875,900,925,950,975,1000
#vars_grb2=["Temperature on pressure surface"]
#vars_grb2=["Relative humidity on pressure surface"]
#vars_grb2=["Dew point temperature on pressure surface"]
#vars_grb2=["Specific humidity on pressure surface"]
#vars_grb2=["U component of wind on pressure surface"]
#vars_grb2=["V component of wind on pressure surface"]
#vars_grb2=["Vertical velocity on pressure surface"]
#vars_grb2=["Absolute vorticity on pressure surface"]
#vars_grb2=["Turbulent kinetic energy on pressure surface"]
#vars_grb2=["Cloud mixing ratio on pressure surface"]
#vars_grb2=["Ice water mixing ratio on pressure surface"]
#vars_grb2=["Rain mixing ratio on pressure surface"]
#vars_grb2=["Snow mixing ratio on pressure surface"]

##### Specific height: 30,50,80,100,305,457,610,914,1524,1829,2134,2743,3658,4572
#vars_grb2=["Temperature at specific height"]
#vars_grb2=["Specific humidity at specific height"]
#vars_grb2=["U component of wind at specific height"]
#vars_grb2=["V component of wind at specific height"]

##### Boundary layers: 3000, 6000, 9000, 12000, 15000, 18000 #####
#vars_grb2=["Temperature in boundary layer"]
#vars_grb2=["Dew point temperature in boundary layer"]
#vars_grb2=["Specific humidity in boundary layer"]
#vars_grb2=["U component of wind in boundary layer"]
#vars_grb2=["V component of wind in boundary layer"]
#vars_grb2=["Vertical velocity in boundary layer"]

##### Single layer #####
#vars_grb2=["U component of storm motion"]
#vars_grb2=["V component of storm motion"]
#vars_grb2=["2 meter temperature"]
#vars_grb2=["2 meter dew point temperature"]
#vars_grb2=["2 meter relative humidity"]
#vars_grb2=["10 meter U wind component"]
#vars_grb2=["10 meter V wind component"]
#vars_grb2=["10 meter wind speed"]
#vars_grb2=["Total cloud cover"]
#vars_grb2=["Total precipitation"]
#vars_grb2=["Storm surface runoff"]
#vars_grb2=["Sea surface temperature"]
#vars_grb2=["Composite radar reflectivity"]
#vars_grb2=["Orography"]

vars_grb2=["Total precipitation","Composite radar reflectivity"]

# Layer number to be plotted
ilvl=1

# basic forms of title and file name
out_title_base='FV3LAM::BGDAWP::'+fnm_hr+'::'
out_fname_base='fv3lam_out_bgdawp_'+fnm_hr+'_'

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

    print(' ===== '+svar+' ===== prs ===============================')

    nm_svar=svar
    cs_cmap='jet'
    cs_norm=None
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=2
    cmap_range='round'

    if svar=="Temperature on pressure surface":
        grbv=grbs.select(name="Temperature",typeOfLevel="isobaricInhPa")[ilvlm]
    elif svar=="Temperature at specific height":
        grbv=grbs.select(name="Temperature",typeOfLevel="heightAboveGround")[ilvlm]
    elif svar=="Temperature in boundary layer":
        grbv=grbs.select(name="Temperature",typeOfLevel="pressureFromGroundLayer")[ilvlm]
    elif svar=="Relative humidity on pressure surface":
        grbv=grbs.select(name="Relative humidity",typeOfLevel="isobaricInhPa")[ilvlm]
        n_rnd=5
    elif svar=="Dew point temperature on pressure surface":
        grbv=grbs.select(name="Dew point temperature",typeOfLevel="isobaricInhPa")[ilvlm]
    elif svar=="Specific humidity on pressure surface":
        grbv=grbs.select(name="Specific humidity",typeOfLevel="isobaricInhPa")[ilvlm]
        n_rnd=6
    elif svar=="Specific humidity at specific height":
        grbv=grbs.select(name="Specific humidity",typeOfLevel="heightAboveGround")[ilvlm]
        n_rnd=4
    elif svar=="Specific humidity in boundary layer":
        grbv=grbs.select(name="Specific humidity",typeOfLevel="pressureFromGroundLayer")[ilvlm]
        n_rnd=4
    elif svar=="U component of wind on pressure surface":
        grbv=grbs.select(name="U component of wind",typeOfLevel="isobaricInhPa")[ilvlm]
        cmap_range='symmetry'
    elif svar=="U component of wind at specific height":
        grbv=grbs.select(name="U component of wind",typeOfLevel="heightAboveGround")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="U component of wind in boundary layer":
        grbv=grbs.select(name="U component of wind",typeOfLevel="pressureFromGroundLayer")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="U component of storm motion":
        grbv=grbs.select(name="U-component storm motion",typeOfLevel="heightAboveGroundLayer")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="V component of wind on pressure surface":
        grbv=grbs.select(name="V component of wind",typeOfLevel="isobaricInhPa")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="V component of wind at specific height":
        grbv=grbs.select(name="V component of wind",typeOfLevel="heightAboveGround")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="V component of wind in boundary layer":
        grbv=grbs.select(name="V component of wind",typeOfLevel="pressureFromGroundLayer")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="V component of storm motion":
        grbv=grbs.select(name="V-component storm motion",typeOfLevel="heightAboveGroundLayer")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="Vertical velocity in boundary layer":
        grbv=grbs.select(name="Vertical velocity",typeOfLevel="pressureFromGroundLayer")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="Absolute vorticity on pressure surface":
        grbv=grbs.select(name="Absolute vorticity",typeOfLevel="isobaricInhPa")[ilvlm]
        n_rnd=5
    elif svar=="Turbulent kinetic energy on pressure surface":
        grbv=grbs.select(name="Turbulent kinetic energy",typeOfLevel="isobaricInhPa")[ilvlm]
        n_rnd=5
    elif svar=="Cloud mixing ratio on pressure surface":
        grbv=grbs.select(name="Cloud mixing ratio",typeOfLevel="isobaricInhPa")[ilvlm]
        n_rnd=4
    elif svar=="Ice water mixing ratio on pressure surface":
        grbv=grbs.select(name="Ice water mixing ratio",typeOfLevel="isobaricInhPa")[ilvlm]
        n_rnd=4
    elif svar=="Rain mixing ratio on pressure surface":
        grbv=grbs.select(name="Rain mixing ratio",typeOfLevel="isobaricInhPa")[ilvlm]
        n_rnd=4
    elif svar=="Snow mixing ratio on pressure surface":
        grbv=grbs.select(name="Snow mixing ratio",typeOfLevel="isobaricInhPa")[ilvlm]
        n_rnd=4
    elif svar=="2 meter temperature":
        ilvlm=0
        grbv=grbs.select(name="2 metre temperature",typeOfLevel="heightAboveGround")[ilvlm]     
    elif svar=="2 meter dew point temperature":
        ilvlm=0
        grbv=grbs.select(name="2 metre dewpoint temperature",typeOfLevel="heightAboveGround")[ilvlm]
    elif svar=="2 meter relative humidity":
        ilvlm=0
        grbv=grbs.select(name="2 metre relative humidity",typeOfLevel="heightAboveGround")[ilvlm]
        n_rnd=3
    elif svar=="10 meter U wind component":
        ilvlm=0
        grbv=grbs.select(name="10 metre U wind component",typeOfLevel="heightAboveGround")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="10 meter V wind component":
        ilvlm=0
        grbv=grbs.select(name="10 metre V wind component",typeOfLevel="heightAboveGround")[ilvlm]
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="10 meter wind speed":
        ilvlm=0
        grbv=grbs.select(name="10 metre wind speed",typeOfLevel="heightAboveGround")[ilvlm]
    elif svar=="Total precipitation":
        ilvlm=0
        grbv=grbs.select(name="Total Precipitation",typeOfLevel="surface")[ilvlm]
        lb_ext='max'
        cs_cmap,cs_norm=new_cmap()
        cmap_range='designed'
    elif svar=="Storm surface runoff":
        ilvlm=0
        grbv=grbs.select(name="Storm surface runoff",typeOfLevel="surface")[ilvlm]
    elif svar=="Sea surface temperature":
        ilvlm=0
        grbv=grbs.select(name="Sea surface temperature",typeOfLevel="surface")[ilvlm]
    elif svar=="Total cloud cover":
        ilvlm=0
        grbv=grbs.select(name="total cloud cover")[ilvlm]
        cs_cmap="cubehelix_r"
    elif svar=="Composite radar reflectivity":
        ilvlm=0
        grbv=grbs.select(name="Maximum/Composite radar reflectivity")[ilvlm]
        cs_cmap,cs_norm=new_cmap()
        nm_svar=svar+' (dBZ)'
        cmap_range='designed'
        lb_ext='max'
    elif svar=="Orography":
        ilvlm=0
        grbv=grbs.select(name="Orography")[ilvlm]
        cs_cmap="gist_ncar_r"
        n_rnd=0
    else:
        sys.exit('ERROR: Wrong svar or Not set up yet !!! ::'+svar)

    sval=grbv.values
    sname=grbv.name
    print(' Name=',sname)
    stnm=grbv.shortName
    print(' ShortName=',stnm)
    tlxtr=str(grbv.level)
    print(' Level=',tlxtr)
    nmxtr=grbv.typeOfLevel
    print(' Type of level=',nmxtr)
 
    if stnm=='t' or stnm=='dpt' or stnm=='2t' or stnm=='2d' or stnm=='sst':
        # Kelvin to Fahrenheit
        sval=(sval-273.15)*1.8+32.0
        nm_svar=svar+' ('+chr(176)+'F)'
    elif stnm=='u' or stnm=='v' or stnm=='ustm' or stnm=='vstm' or stnm=='w' or stnm=='10u' or stnm=='10v' or stnm=='10si':
        # Meter to Knot (kn)
        sval=sval*1.94384
        nm_svar=svar+' (kn)'


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
            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree(),norm=cs_norm)
    if cmap_range=='designed':
        cs.cmap.set_under('white',alpha=0.)
        cs.cmap.set_over('black')
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

