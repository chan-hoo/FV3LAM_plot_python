###################################################################### CHJ #####
## Name		: plot_fv3lam_co2his.py
## Language	: Python 3.7
## Usage	: Plot historical co2 data files for fv3 regional modeling
## Input files  : co2historicaldata_20XX.txt
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/07/14: Chan-Hoo Jeon : Preliminary version
## V001: 2021/03/05: Chan-Hoo Jeon : Simplify the script
###################################################################### CHJ #####

import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import pandas as pd
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
elif machine=='orion':
    cartopy.config['data_dir']='/home/chjeon/tools/NaturalEarth'
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
else:
    sys.exit('ERROR: Required input data are NOT set !!!')

plt.switch_backend('agg')

# Case-dependent input =============================================== CHJ =====
# Path to the directory where the input file is located.
dnm_data="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/run_C96/"

# input file name
fnm_in='co2historicaldata_2020.txt'

year=fnm_in[-8:-4]
print(' year=',year)

# basic forms of title and file name
out_title_base='FV3LAM::Monthly CO2 in '+year+'::'
out_fname_base='fv3lam_co2his_'

# Resolution of background natural earth data ('50m' or '110m')
back_res='110m'


# Main part (will be called at the end) ==================== CHJ =====
def main():
# ========================================================== CHJ =====

    global lon,lat,extent,c_lon,tmax,tmin

    # open the data file
    fname=os.path.join(dnm_data,fnm_in)
    try: co2h=pd.read_csv(fname,sep='\s+',header=None,skiprows=1,na_values=[-99.99])
    except: raise Exception('Could NOT find the file',fname)

    print(' ===== CO2 history data =======================')
    print(co2h)  

    co2h.shape

    tmax=np.max(np.max(co2h))
    tmin=np.min(np.min(co2h))
    print(' Total max=',tmax)
    print(' Total min=',tmin)

#    lon1d=np.linspace(7.5,352.5,24)
    lon1d=np.linspace(0,345,24)
    lat1d=np.linspace(82.5,-82.5,12)
    print(' lon=',lon1d)
    print(' lat=',lat1d)

    lon,lat=np.meshgrid(lon1d,lat1d,sparse=False)
    #print(lon)
    #print(lat)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon1d)
    lon_max=np.max(lon1d)
    lat_min=np.min(lat1d)
    lat_max=np.max(lat1d)
    
   # extent=[lon_min-5,lon_max+5,lat_min-5,lat_max+3]
    extent=[lon_min,lon_max,lat_min,lat_max]
    c_lon=np.mean(extent[:2])
   # c_lat=np.mean(extent[2:])


    for im in range(0,12):
        imp1=im+1
        print(' month=',imp1)
        im_s=12*im
        co2h_mn=co2h.loc[im_s:im_s+11,:]
        print(co2h_mn.shape)

        data_plot(co2h_mn,imp1)



# Plot data ================================================ CHJ =====
def data_plot(co2h_mn,imp1):
# ========================================================== CHJ =====

    cs_cmap='YlOrBr'
    lb_ext='both'
    n_rnd=2
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3

    mon_nm=format(imp1,'02d')
    out_title_fld=out_title_base+'M'+mon_nm
    out_fld_fname=out_fname_base+'m'+mon_nm
   
    cs_label='Monthly averaged CO2'
    cmap_range='fixed'

    # Max and Min of the field
    fmax=np.max(np.max(co2h_mn))
    fmin=np.min(np.min(co2h_mn))
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
        cs_min=tmin
        cs_max=tmax
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_max=',cs_max)
    print(' cs_min=',cs_min)


    # Plot field
    fig,ax1=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax1.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax1)
    ax1.set_title(out_title_fld,fontsize=9)
    cs=ax1.pcolormesh(lon,lat,co2h_mn,cmap=cs_cmap,rasterized=True,
            vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
#    cs=ax1.contourf(lon,lat,co2h_mn,cmap=cs_cmap,vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())

    # extend(pointed end): 'neither'|'both'|'min'|'max'  
    divider=make_axes_locatable(ax1)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(cs_label,fontsize=8)

    # Output figure
    out_file(out_fld_fname) 

       
   

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

