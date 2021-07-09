###################################################################### CHJ #####
## Name		: plot_fv3lam_phyf.py
## Language	: Python 3.7
## Usage	: Plot an output, phy, for fv3 regional modeling
## Input files  : phyfXXX.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/05/07: Chan-Hoo Jeon : Preliminary version
## V001: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
## V002: 2021/03/05: Chan-Hoo Jeon : Simplify the script
## V003: 2021/06/24: Chan-Hoo Jeon : Add a projection for RRFS_NA domain
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
dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/ufs_srw_app/srw_dev_test/expt_dirs/inline_post/2019070100/"
#dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/ufs_srw_app/srw_dev_test/expt_dirs/grid_RRFS_NA_13km/2019070100/"

# Domain name
domain_nm='RRFS_CONUS_25km'
#domain_nm='RRFS_NA_13km'

# grid file name
fnm_hr='f003'

fnm_input='phy'+fnm_hr+'.nc'

# Variables
#vars_phy=["acond","albdo_ave","alnsf","alnwf","alvsf","alvwf"]
#vars_phy=["cduvb_ave","cnwat","cpofp","cprat_ave","cpratb_ave"]
#vars_phy=["csdlf","csdsf","csulf","csulftoa","csusf","csusftoa"]
#vars_phy=["csusf","csusftoa","cwork_aveclm","dlwrf","dlwrf_ave"]
#vars_phy=["dswrf","dswrf_avetoa","duvb_ave","evbs_ave","evcw_ave"]
#vars_phy=["f10m","facsf","facwf","ffhh","ffmm","fldcp","fricv"]
#vars_phy=["gflux","gflux_ave","hgt_hyblev1","hpbl","icec"]
#vars_phy=["icetk","land","lhtfl","lhtfl_ave","nbdsf_ave"]
#vars_phy=["nddsf_ave","orog","pevpr","pevpr_ave","prateb_ave","presssfc"]
#vars_phy=["pwatclm","rh02max","rh02min","sbsno_ave","sfcr","sfexc","shdmax"]
#vars_phy=["shdmin","shtfl","shtfl_ave","sltyp","snoalb"]
#vars_phy=["snod","snohf","snowc_ave","soill1","soill2"]
#vars_phy=["soill3","soill4","soilm","soilt1","soilt2"]
#vars_phy=["soilt3","soilt4","soilw1","soilw2","soilw3"]
#vars_phy=["soilw4","sotyp","spd10max","spfh2m","spfh_hyblev1"]
#vars_phy=["spfhmax_max2m","spfhmin_min2m","ssrun_acc","sunsd_acc"]
#vars_phy=["t02max","t02min","tg3","tisfc"]
#vars_phy=["tmax_max2m","tmin_min2m","tmp2m","tmpsfc","tprcp","trans_ave"]
#vars_phy=["u-gwd_ave","u10max","uflx","uflx_ave","ugrd10m"]
#vars_phy=["ugrd_hyblev1","ulwrf","ulwrf_ave","ulwrf_avetoa","uswrf"]
#vars_phy=["uswrf_ave","uswrf_avetoa","v-gwd_ave","v10max","vbdsf_ave"]
#vars_phy=["vddsf_ave","veg","vflx","vflx_ave","vgrd10m"]
#vars_phy=["vgrd_hyblev1","vtype","watr_acc","weasd","wilt"]
vars_phy=["tmpsfc"]

# basic forms of title and file name
out_title_base='FV3LAM::PHY::'+domain_nm+'::'
out_fname_base='fv3lam_out_phy_'+domain_nm+'_'

# grid plot
plot_grid='on'

# Grid point plot (every 'n_skip' rows/columns) 
n_skip=5

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'



# Main part (will be called at the end) ======================= CHJ =====
def main():
# ============================================================= CHJ =====
    global phyf,lon,lat
    global extent,c_lon,c_lat
    global lon_min,lat_min

    print(' ===== DATA: phy =======================================')
    # open the data file
    fname=os.path.join(dnm_data,fnm_input)
    try: phyf=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(phyf)

    lon=np.ma.masked_invalid(phyf["lon"].data)
    lat=np.ma.masked_invalid(phyf["lat"].data)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(lon)
    lon_max=np.max(lon)
    lat_min=np.min(lat)
    lat_max=np.max(lat)

    print(' lon_min=',lon_min,', lon_max=',lon_max)
    print(' lat_min=',lat_min,', lat_max=',lat_max)

    # Plot extent
    esp=1
    extent=[lon_min-esp,lon_max+esp,lat_min-esp,lat_max+esp]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])

    # Grid plot
    if plot_grid=='on':
        grid_plot()

    # Variables
    for svar in vars_phy:
        phy_plot(svar)


# ===== grid plot ============================================= CHJ =====
def grid_plot():
# ============================================================= CHJ =====

    # grid points (every 'n_skip' rows/columns from 2nd row/col)
    grdx_slc=lon[::n_skip,::n_skip]
    grdy_slc=lat[::n_skip,::n_skip]

    if domain_nm[:7]=='RRFS_NA':
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(
                            central_longitude=-107,central_latitude=53)))
        ref_lon=-133.5
        ref_lat=lat_min-5.5
        lgd_loc='lower left'
    else:
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
        ax.set_extent(extent, ccrs.PlateCarree())
        ref_lon=lon_min
        ref_lat=lat_min-3
        lgd_loc='lower right'

    out_grd_title=out_title_base+'Grid'
    out_grd_fname=out_fname_base+'grid' 

    ax.set_title(out_grd_title, fontsize=9)

    # Call background plot
    back_plot(ax)

    # Scatter plot (zorder: lowest-plot on bottom, highest-plot on top)
    sp_scale=2
    s1=ax.scatter(grdx_slc,grdy_slc,transform=ccrs.PlateCarree(),marker='o',facecolors="None",edgecolors='red',linewidth=0.3,s=sp_scale,zorder=3)

    ref_txt='Grid: every '+str(n_skip)+'th col/row'
    plt.text(ref_lon,ref_lat,ref_txt,transform=ccrs.Geodetic(),fontsize=8)
#    plt.legend((s1),('super-grid'),scatterpoints=1,loc=lgd_loc,ncol=1,fontsize=8)

    # Output figure
    ndpi=300
    out_file(out_grd_fname,ndpi)




# ===== plot ================================================== CHJ =====
def phy_plot(svar):
# ============================================================= CHJ =====

    print(' ===== '+svar+' ===== phy ===============================')
    # Extract data array
    sfld=np.ma.masked_invalid(phyf[svar].data)

    (nts,nys,nxs)=sfld.shape
    print(' time+2D: nts=',nts,' nys=',nys,' nxs=',nxs)
    sfld2d=np.squeeze(sfld,axis=0)
    out_title_fld=out_title_base+svar+'::'+fnm_hr
    out_phy_fname=out_fname_base+svar+'_'+fnm_hr

    nm_svar=svar
    cs_cmap='jet'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=2
    cmap_range='round'
 
    if svar=="sltyp":
        nm_svar='Surface slope type'
        cs_cmap=plt.cm.get_cmap('Paired',10)
        n_rnd=0
    elif svar=="snoalb":
        nm_svar='Max. snow albedo over land'
        cs_cmap='nipy_spectral_r'
        n_rnd=2
    elif svar=="sotyp":
        nm_svar='Soil type'
        cs_cmap=plt.cm.get_cmap('tab20',17)
        n_rnd=0
    elif svar=="tprcp":
        nm_svar='Precipitation rate'
        n_rnd=4
    elif svar=="ugrd10m":
        nm_svar='u-comp. of 10m wind'
        n_rnd=2
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="ugrd_hyblev1":
        nm_svar='u-comp. of current of hybrid level 1'
        n_rnd=2
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="vgrd10m":
        nm_svar='v-comp. of 10m wind'
        n_rnd=2
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="vgrd_hyblev1":
        nm_svar='v-comp. of current of hybrid level 1'
        n_rnd=2
        cmap_range='symmetry'
        lb_ext='both'
        cs_cmap='seismic'
    elif svar=="vtype":
        nm_svar='Vegetation type'
        cs_cmap=plt.cm.get_cmap('tab20',20)
        n_rnd=0


    if svar=="tmp2m":
        sfld2d=(sfld2d-273.15)*1.8+32.0


    print(' PHY field=',nm_svar)

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
    if domain_nm[:7]=='RRFS_NA':
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(
                            central_longitude=-107,central_latitude=53)))
    else:
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
        ax.set_extent(extent, ccrs.PlateCarree())

    back_plot(ax)
    ax.set_title(out_title_fld,fontsize=9)
    cs=ax.pcolormesh(lon,lat,sfld2d,cmap=cs_cmap,rasterized=True,
        vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(nm_svar,fontsize=6)

    # Output figure
    ndpi=300
    out_file(out_phy_fname,ndpi)

  

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

