###################################################################### CHJ #####
## Name		: plot_fv3lam_ani_comref.py
## Language	: Python 3.7
## Usage	: Create animation of hourly comparison of reflectivity by fv3 and mrms radar
## Input files  : BGDAWPXX.tmXX(*.prs.grib2) and QCComposite.XXX.grib2
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/08/03: Chan-Hoo Jeon : Preliminary version
## V001: 2020/08/05: Chan-Hoo Jeon : Fix the issue of colorbar duplication
###################################################################### CHJ #####

import os, sys
import pygrib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER,LATITUDE_FORMATTER
import matplotlib.animation as animation

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
# Domain name:
domain='HRRR'
# Resolution:
res='25km'

# date
s_date='20200618'
# plotting hours for animation
plot_hrs=['01','02','03','04','05','06']

icmax=len(plot_hrs)

# Path to the directory where the input NetCDF files are located.
dnm_data_mdl="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/expt_dirs/test_community_hrrr25/2020061800/postprd/"
dnm_data_dat="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/TMP/mrms_"+s_date+"/new_grib2/"


# Output fields
svar="Reflectivity"


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

# Title and file name
out_title_sub1='FV3-LAM'
out_title_sub2='MRMS Radar Data'

out_fname='fv3_out_ani_refl_comp_'+domain

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'
# high-resolution background image ('on', 'off')
back_img='off'



# Main part (will be called at the end) ==================== CHJ =====
def main():
# ========================================================== CHJ =====

    read_data()
    plot_comp()



# Read data ================================================ CHJ =====
def read_data():
# ========================================================== CHJ =====
  
    global lon1,lat1,lon2,lat2,extent,sval_mdl,sval_dat

    for ic in range(0,icmax):
        shr=plot_hrs[ic]
        print(ic,shr)
        # Input file (BGDAWPXX.tmXX)
        fnm_in_mdl=domain+'.t00z.bgdawp'+shr+'.tm00'
        fnm_in_dat='QCComposite_00.50_'+s_date+shr+'00_new.grib2'

        print(' ===== Read result file ===== '+plot_hrs[ic]+' ==========================')

        # open the data file
        fname=os.path.join(dnm_data_mdl,fnm_in_mdl)
        try: grbs=pygrib.open(fname)
        except: raise Exception('Could NOT find the file',fname)

        grbv1=grbs.select(name="Maximum/Composite radar reflectivity")[0]
        sval=grbv1.values
        sval[sval<0]=0
        if ic==0:
            Nx1=grbv1.Nx
            Ny1=grbv1.Ny
            sval_mdl=np.zeros((Ny1,Nx1,icmax))
        sval_mdl[:,:,ic]=sval

        print(' ===== Read radar data ===== '+plot_hrs[ic]+' ==============================')
        # open the data file
        fname=os.path.join(dnm_data_dat,fnm_in_dat)
        try: grbs=pygrib.open(fname)
        except: raise Exception('Could NOT find the file',fname)

        grbv2=grbs.select(typeOfLevel="heightAboveSea",level=500)[0]
        sval=grbv2.values
        sval[sval<0]=0
        if ic==0:
            Nx2=grbv2.Nx
            Ny2=grbv2.Ny
            sval_dat=np.zeros((Ny2,Nx2,icmax))
        sval_dat[:,:,ic]=sval

        if ic==0:
            sname=grbv1.name
            print(' Name=',sname)
            stnm=grbv1.shortName
            print(' ShortName=',stnm)
            lat1,lon1=grbv1.latlons()
            # Highest and lowest longitudes and latitudes for plot extent
            lon_min=np.min(lon1)
            lon_max=np.max(lon1)
            lat_min=np.min(lat1)
            lat_max=np.max(lat1)
            print(' lon_min=',lon_min,', lon_max=',lon_max)
            print(' lat_min=',lat_min,', lat_max=',lat_max)
            esp=1
            extent=[lon_min-esp,lon_max+esp,lat_min-esp,lat_max+esp]
            lat2,lon2=grbv2.latlons()



# Plot comparison ========================================== CHJ =====
def plot_comp():
# ========================================================== CHJ =====
    global fig,spec_fig,ax1,ax2,lb_fnt,cs_cmap,cs_norm,cs_min,cs_max

    ani_fps=1
    lb_fnt=3

    cs_cmap,cs_norm=new_cmap()
    nm_svar='Reflectivity (dBZ)'
    lb_ext='max'
    cs_min=5
    cs_max=None

    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])

    # Plot field
    fig=plt.figure(figsize=(3,2.8))   #(width,height)
    spec_fig=fig.add_gridspec(2,1)
    out_title_super='Composite Reflectivity::'+domain+'('+res+')::'+s_date+'/  '
    fig.suptitle(out_title_super,fontsize=lb_fnt+2)
    ax1=fig.add_subplot(spec_fig[0,0],projection=ccrs.PlateCarree(c_lon))
    ax2=fig.add_subplot(spec_fig[1,0],projection=ccrs.PlateCarree(c_lon))

    cbar_plot(fig,ax1,nm_svar,lb_ext,lb_fnt)
    cbar_plot(fig,ax2,nm_svar,lb_ext,lb_fnt)


    print(' Animation working ...')
    anim=animation.FuncAnimation(fig,animate,init_func=init,frames=icmax,interval=1/ani_fps,blit=False)
    print(' Animation complete ... and saving ...')
    anim.save(out_fig_dir+out_fname+'.gif',dpi=300,fps=ani_fps,writer='imagemagick')




# Animation ================================================ CHJ =====
def animate(ic):
# ========================================================== CHJ =====
    shr=plot_hrs[ic]
    out_title_super='Composite Reflectivity::'+domain+'('+res+')::'+s_date+'/'+shr
    fig.suptitle(out_title_super,fontsize=lb_fnt+2)
    # subplot 1: composite reflectivity + target point lines
    ax1.clear()
    ax1.set_extent(extent,ccrs.PlateCarree())
    back_plot(ax1)
    ax1.set_title(out_title_sub1,fontsize=lb_fnt+1,color='red')
    cs1=ax1.pcolormesh(lon1,lat1,sval_mdl[:,:,ic],cmap=cs_cmap,rasterized=True,
            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree(),norm=cs_norm)
    cs1.cmap.set_under('white',alpha=0.)
    cs1.cmap.set_over('black')
    gridline_plot(ax1,lb_fnt)

    # subplot 2: radar data
    ax2.clear()
    ax2.set_extent(extent,ccrs.PlateCarree())
    back_plot(ax2)
    ax2.set_title(out_title_sub2,fontsize=lb_fnt+1,color='red')
    cs2=ax2.pcolormesh(lon2,lat2,sval_dat[:,:,ic],cmap=cs_cmap,rasterized=True,
            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree(),norm=cs_norm)
    cs2.cmap.set_under('white',alpha=0.)
    cs2.cmap.set_over('black')
    gridline_plot(ax2,lb_fnt)


# Initialization for animation ============================= CHJ =====
def init():
# ========================================================== CHJ =====
    return 


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



# gridline plot ========================================== CHJ =====
def gridline_plot(ax,lb_sz):
# ========================================================== CHJ ===== 
    gl=ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                    linewidth=0.2,color='chocolate',alpha=0.3,linestyle='-')
    gl.xlabels_top=False
    gl.ylabels_right=False
    gl.xlocator=mticker.FixedLocator([-130,-120,-110,-100,-90,-80,-70,-60])
    gl.ylocator=mticker.FixedLocator([20,30,40,50,60])
    gl.xformatter=LONGITUDE_FORMATTER
    gl.yformatter=LATITUDE_FORMATTER
    gl.xlabel_style={'size':lb_sz,'color':'black'}
    gl.ylabel_style={'size':lb_sz,'color':'black'}



# color bar plot ========================================== CHJ =====
def cbar_plot(fig,ax,nm_svar,lb_ext,lb_sz):
# ========================================================== CHJ =====
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=matplotlib.colorbar.ColorbarBase(ax_cb,cmap=cs_cmap,norm=cs_norm,orientation='vertical',extend=lb_ext)
    cbar.ax.tick_params(labelsize=lb_sz)
    cbar.set_label(nm_svar,fontsize=lb_sz)



# Output file ============================================= CHJ =====
def out_file(out_file,ndpi):
# ========================================================= CHJ =====
    # Output figure
    plt.savefig(out_fig_dir+out_file+'.png',dpi=ndpi,bbox_inches='tight')
    plt.close('all')



# Main call ================================================ CHJ =====
if __name__=='__main__':
    main()

