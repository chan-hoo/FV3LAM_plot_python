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
dnm_out="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/expt_dirs/test_community/2020122700/"

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
out_title_base='FV3LAM::PHY::'
out_fname_base='fv3lam_out_phy_'

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'



# Main part (will be called at the end) ======================= CHJ =====
def main():
# ============================================================= CHJ =====
    global phyf,lon,lat
    global extent,c_lon,c_lat

    print(' ===== OUTPUT: phy =======================================')
    # open the data file
    fname=os.path.join(dnm_out,fnm_input)
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

    # Variables
    for svar in vars_phy:
        phy_plot(svar)



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
    out_dyn_fname=out_fname_base+svar+'_'+fnm_hr

    nm_svar=svar
    cs_cmap='jet'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=5
    cmap_range='round'
 
    if svar=="acond":
        nm_svar='Aerodynamic conductance'
        n_rnd=3
    elif svar=="albdo_ave":
        nm_svar='Albedo:average'
        n_rnd=2
    elif svar=="alnsf":
        nm_svar='Near-IR black sky albedo'
        cs_cmap='YlGnBu'
        n_rnd=3
    elif svar=="alnwf":
        nm_svar='Near-IR white sky albedo'
        cs_cmap='YlOrRd'
        n_rnd=3
    elif svar=="alvsf":
        nm_svar='Visible black sky albedo'
        cs_cmap='YlGnBu'
        n_rnd=3
    elif svar=="alvwf":
        nm_svar='Visible white sky albedo'
        cs_cmap='YlOrRd'
        n_rnd=3
    elif svar=="cduvb_ave":
        nm_svar='Clear sky UV-B solar flux'
        n_rnd=2
    elif svar=="cnwat":
        nm_svar='Canopy moisture content'
        n_rnd=2
    elif svar=="cpofp":
        nm_svar='Percent of frozen precipitation'
        n_rnd=2
    elif svar=="cprat_ave":
        nm_svar='Convective precipitatin'
        n_rnd=5
    elif svar=="cpratb_ave":
        n_rnd=5
    elif svar=="csdlf":
        nm_svar='Clear sky downward long wave flux'
        n_rnd=2
    elif svar=="csdsf":
        nm_svar='Clear sky downward solar flux'
        n_rnd=2
    elif svar=="csulf":
        nm_svar='Clear sky upward long wave flux'
        n_rnd=2
    elif svar=="csulftoa":
        nm_svar='csulf: top of atmos.'
        n_rnd=2
    elif svar=="csusf":
        nm_svar='Clear sky upward solar flux'
        n_rnd=2
    elif svar=="csusftoa":
        nm_svar='csusf: top of atmos.'
        n_rnd=2
    elif svar=="cwork_aveclm":
        nm_svar='Cloud work function'
        n_rnd=2
    elif svar=="dlwrf":
        nm_svar='Downward long wave rad. flux'
        n_rnd=2
    elif svar=="dlwrf_ave":
        nm_svar='dlwrf: average'
        n_rnd=2
    elif svar=="dswrf":
        nm_svar='Downward short wave rad. flux'
        n_rnd=2
    elif svar=="dswrf_ave":
        nm_svar='dswrf: average'
        n_rnd=2
    elif svar=="dswrf_avetoa":
        nm_svar='dswrf:average:top of atmos'
        n_rnd=2
    elif svar=="duvb_ave":
        nm_svar='UV-B downward solar flux'
        n_rnd=2
    elif svar=="evbs_ave":
        nm_svar='Direct evaporation from bair soil'
        n_rnd=2
    elif svar=="evcw_ave":
        nm_svar='Canopy water evaporation'
        n_rnd=2
    elif svar=="f10m":
        nm_svar='10m wind speed over lowest value'
        n_rnd=2
    elif svar=="facsf":
        nm_svar='Fract. coverage w/ strong cosz'
        n_rnd=2
    elif svar=="facwf":  
        nm_svar='fract. coverage w/ weak cosz'
        n_rnd=2
    elif svar=="ffhh":
        nm_svar='Surface exchange coeff. for heat'
        n_rnd=2
    elif svar=="ffmm":
        nm_svar='Surface exchange coeff. for momentum'
        n_rnd=2
    elif svar=="fldcp":
        nm_svar='Field capacity fraction'
        n_rnd=2
    elif svar=="fricv":
        nm_svar='Friction velocity'
        n_rnd=3
    elif svar=="gflux":
        nm_svar='Ground heat flux'
        n_rnd=2
    elif svar=="gflux_ave":
        nm_svar='Ground heat flud: average'
        n_rnd=2
    elif svar=="hgt_hyblev1":
        nm_svar='Height of hybrid level 1'
        n_rnd=2
    elif svar=="hpbl":
        nm_svar='Planetary boundary layer height'
        n_rnd=2
    elif svar=="icec":
        nm_svar='Sea-ice fraction'
        n_rnd=1
    elif svar=="icetk":
        nm_svar='Sea-ice depth'
        n_rnd=2
    elif svar=="land":
        nm_svar='Land cover'
        cs_cmap=plt.cm.get_cmap('Paired',3)
        n_rnd=0
    elif svar=="lhtfl":
        nm_svar='Latent heat net flux'
        n_rnd=2
    elif svar=="lhtfl_ave":
        nm_svar='Latent heat net flux:avg'
        n_rnd=2
    elif svar=="nbdsf_ave":
        nm_svar='Near IR diffuse dw. solar flux'
        n_rnd=2
    elif svar=="nddsf_ave":
        n_rnd=2
    elif svar=="orog":
        nm_svar='Surface altitude'
        cs_cmap='terrain_r'
        n_rnd=2
    elif svar=="pevpr":
        nm_svar='Potential evaporation rate'
        n_rnd=2
    elif svar=="pevpr_ave":
        nm_svar='Potential evaporation rate:avg'
        n_rnd=2
    elif svar=="prate_ave":
        n_rnd=4
    elif svar=="prateb_ave":
        n_rnd=5
    elif svar=="pressfc":
        nm_svar='Surface pressure'
        n_rnd=2
    elif svar=="pwatclm":
        nm_svar='Precipitable water'
        n_rnd=2
    elif svar=="rh02max":
        n_rnd=3
    elif svar=="rh02min":
        n_rnd=3
    elif svar=="sbsno_ave":
        nm_svar='Sublimation'
        n_rnd=2
    elif svar=="sfcr":
        nm_svar='Roughness length'
        cs_cmap='gist_earth_r'
        n_rnd=2
    elif svar=="sfexc":
        nm_svar='Exchange coefficient'
        n_rnd=2
    elif svar=="shdmax":
        nm_svar='Maximum areal fractional coverage'
        n_rnd=2
    elif svar=="shdmin":
        nm_svar='Minimum areal fractional coverage'
        n_rnd=2
    elif svar=="shtfl":
        nm_svar='Sensible heat flux'
        n_rnd=2
    elif svar=="shtfl_ave":
        nm_svar='Sensible heat flux: average'
        n_rnd=2
    elif svar=="sltyp":
        nm_svar='Surface slope type'
        cs_cmap=plt.cm.get_cmap('Paired',10)
        n_rnd=0
    elif svar=="snoalb":
        nm_svar='Max. snow albedo over land'
        cs_cmap='nipy_spectral_r'
        n_rnd=2
    elif svar=="snod":
        nm_svar='Physical snow depth'
        n_rnd=2
    elif svar=="snohf":
        nm_svar='Snow phase-change heat flux'
        n_rnd=2
    elif svar=="snowc_ave":
        nm_svar='Snow cover: average'
        n_rnd=1
    elif svar=="soill1":
        nm_svar='Vol. fract. of unfrozen soil moist. 1'
        n_rnd=3
    elif svar=="soill2":
        nm_svar='Vol. fract. of unfrozen soil moist. 2'
        n_rnd=3
    elif svar=="soill3":
        nm_svar='Vol. fract. of unfrozen soil moist. 3'
        n_rnd=3
    elif svar=="soill4":
        nm_svar='Vol. fract. of unfrozen soil moist. 4'
        n_rnd=3
    elif svar=="soilm":
        nm_svar='Soil moisture content'
        n_rnd=2
    elif svar=="soilt1":
        nm_svar='Soil column temperature 1'
        n_rnd=2
    elif svar=="soilt2":
        nm_svar='Soil column temperature 2'
        n_rnd=2
    elif svar=="soilt3":
        nm_svar='Soil column temperature 3'
        n_rnd=2
    elif svar=="soilt4":
        nm_svar='Soil column temperature 4'
        n_rnd=2
    elif svar=="soilw1":
        nm_svar='Total volumetric soil moisture 1'
        n_rnd=2
    elif svar=="soilw2":
        nm_svar='Total volumetric soil moisture 2'
        n_rnd=2
    elif svar=="soilw3":
        nm_svar='Total volumetric soil moisture 3'
        n_rnd=2
    elif svar=="soilw4":
        nm_svar='Total volumetric soil moisture 4'
        n_rnd=2
    elif svar=="sotyp":
        nm_svar='Soil type'
        cs_cmap=plt.cm.get_cmap('tab20',17)
        n_rnd=0
    elif svar=="spd10max":
        nm_svar='10m current speed'
        n_rnd=2
    elif svar=="spfh2m":
        nm_svar='2m specific humidity'
        n_rnd=3
    elif svar=="spfh_hyblev1":
        nm_svar='Specific humidity of hybrid level 1'
        n_rnd=3
    elif svar=="spfhmax_max2m":
        n_rnd=3
    elif svar=="spfhmin_min2m":
        n_rnd=3
    elif svar=="ssrun_acc":
        nm_svar='Storm surface runoff'
        n_rnd=2
    elif svar=="sunsd_acc":
        nm_svar='Sunshine duration'
        n_rnd=1
    elif svar=="t02max":
        nm_svar='Max. of 2m temperature'
        n_rnd=2
    elif svar=="t02min":
        nm_svar='Min. of 2m temperature'
        n_rnd=2
    elif svar=="tg3":
        nm_svar='Deep soil temperature'
        n_rnd=2
    elif svar=="tisfc":
        nm_svar='Sea-ice surface temperature'
        n_rnd=2
    elif svar=="tmax_max2m":
        n_rnd=2
    elif svar=="tmin_min2m":
        n_rnd=2
    elif svar=="tmp2m":
        nm_svar='2m temperatue'
        n_rnd=2
    elif svar=="tmpsfc":
        nm_svar='Sea surface temperature'
        n_rnd=2
    elif svar=="tprcp":
        nm_svar='Precipitation rate'
        n_rnd=4
    elif svar=="trans_ave":
        nm_svar='Transpiration'
        n_rnd=2
    elif svar=="u-gwd_ave":
        nm_svar='Zonal flux of gravity wave stress'
        n_rnd=2
    elif svar=="u10max":
        nm_svar='Max. of u-comp. of 10m wind'
        n_rnd=2
    elif svar=="uflx":
        nm_svar='u-comp. of momentum flux'
        n_rnd=2
    elif svar=="uflx_ave":
        nm_svar='u-comp. of momentum flux:avg'
        n_rnd=2
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
    elif svar=="ulwrf":
        nm_svar='Upward long wave rad. flux'
        n_rnd=2
    elif svar=="ulwrf_ave":
        nm_svar='Upward long wave rad. flux:avg'
        n_rnd=2
    elif svar=="ulwrf_avetoa":
        nm_svar='ulwrf:avg at top of atmos.'
        n_rnd=2
    elif svar=="uswrf":
        nm_svar='Upward short wave rad. flux'
        n_rnd=2
    elif svar=="uswrf_ave":
        nm_svar='Upward short wave rad. flux:avg'
        n_rnd=2
    elif svar=="uswrf_avetoa":
        nm_svar='uswrf:avg at top of atmos.'
        n_rnd=2
    elif svar=="v-gwd_ave":
        nm_svar='Meridional flux of gravity wave stress'
        n_rnd=2
    elif svar=="v10max":
        nm_svar='Max. of v-comp. of 10m wind'
        n_rnd=2
    elif svar=="vbdsf_ave":
        nm_svar='Visible beam downward solar flux'
        n_rnd=2
    elif svar=="vddsf_ave":
        nm_svar='Visible diffuse downward solar flux'
        n_rnd=2
    elif svar=="veg":
        nm_svar='Vegetation'
        n_rnd=2
    elif svar=="vflx":
        nm_svar='v-comp. of momentum flux'
        n_rnd=2
    elif svar=="vflx_ave":
        nm_svar='v-comp. of momentum flux:avg'
        n_rnd=2
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
    elif svar=="watr_acc":
        nm_svar='Water runoff'
        n_rnd=2
    elif svar=="weasd":
        nm_svar='Snow-liquid equivalent'
        n_rnd=1
    elif svar=="wilt":
        nm_svar='Wilting point fraction'
        n_rnd=2
    else:
        sys.exit('ERROR: wrong svar !!!')

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
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title_fld,fontsize=9)
    cs=ax.pcolormesh(lon,lat,sfld2d,cmap=cs_cmap,rasterized=True,
        vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(nm_svar,fontsize=8)

    # Output figure
    ndpi=300
    out_file(out_dyn_fname,ndpi)

  

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

