###################################################################### CHJ #####
## Name		: plot_rrfscmaq_bio.py
## Language	: Python 3.7
## Usage	: Plot an input bio file for rrfs_cmaq
## Input files  : BEIS_CRES.ncf
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2021/06/29: Chan-Hoo Jeon : Preliminary version
###################################################################### CHJ #####

import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from netCDF4 import Dataset
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
dnm_in="/scratch2/NCEPDEV/naqfc/RRFS_CMAQ/aqm/bio/"

# Domain name
domain_nm="GSD_HRRR_25km"
#domain_nm="RRFS_CONUS_13km"

# grid file name
if domain_nm=="GSD_HRRR_25km":
    fnm_input="BEIS_SARC401.ncf"
elif domain_nm=="RRFS_CONUS_13km":
    fnm_input="BEIS_RRFScmaq_C775.ncf"
else:
    sys.exit('ERROR: Wrong domain name !!!')

# Variables
vars_data=["AVG_NOAG_GROW","AVG_NOAG_NONGROW","AVG_NONONAG","AVG_ACETS","AVG_ACTACS",
           "AVG_ACTALS","AVG_APINS","AVG_ATERPS","AVG_ATHUS","AVG_BPHES","AVG_BPINS",
           "AVG_BUTES","AVG_BUTOS","AVG_CAMPHS","AVG_COS","AVG_D3CARS","AVG_DLIMS",
           "AVG_ETHAS","AVG_ETHES","AVG_ETHOS","AVG_FORACS","AVG_FORMS","AVG_GTERPS",
           "AVG_HEXAS","AVG_HEXES","AVG_HEXYS","AVG_ISOPS","AVG_MBOS","AVG_METHS",
           "AVG_MYRCS","AVG_OCIMS","AVG_ORVOCS","AVG_PCYMS","AVG_PROPES","AVG_SABIS",
           "AVG_SESQTS","AVG_TRPOS","LAI_ISOPS","LAI_MBOS","LAI_METHS","AVG_ACETW",
           "AVG_ACTACW","AVG_ACTALW","AVG_APINW","AVG_ATERPW","AVG_ATHUW","AVG_BPHEW",
           "AVG_BPINW","AVG_BUTEW","AVG_BUTOW","AVG_CAMPHW","AVG_COW","AVG_D3CARW",
           "AVG_DLIMW","AVG_ETHAW","AVG_ETHEW","AVG_ETHOW","AVG_FORACW","AVG_FORMW",
           "AVG_GTERPW","AVG_HEXAW","AVG_HEXEW","AVG_HEXYW","AVG_ISOPW","AVG_MBOW",
           "AVG_METHW","AVG_MYRCW","AVG_OCIMW","AVG_ORVOCW","AVG_PCYMW","AVG_PROPEW",
           "AVG_SABIW","AVG_SESQTW","AVG_TRPOW","LAI_ISOPW","LAI_MBOW","LAI_METHW"]

# basic forms of title and file name
out_title_base='RRFS_CMAQ::Bio::'+domain_nm+'::'
out_fname_base='rrfscmaq_bio_'+domain_nm+'_'

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'



# Main part (will be called at the end) ======================= CHJ =====
def main():
# ============================================================= CHJ =====
    global ds,lon,lat
    global extent,c_lon,c_lat

    print(' ===== Input data ========================================')
    # open the data file
    fname=os.path.join(dnm_in,fnm_input)
    try: ds=Dataset(fname,'r')
    except: raise Exception('Could NOT find the file',fname)
    print(ds)
 
    lon=ds.variables['XLONG']
    lat=ds.variables['XLAT']

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
    for svar in vars_data:
        data_plot(svar)



# ===== plot ================================================== CHJ =====
def data_plot(svar):
# ============================================================= CHJ =====

    print(' ===== '+svar+' ===== data ===============================')
    # Extract data array
    sfld=ds.variables[svar][:]

    ndim_svar=sfld.ndim

    if ndim_svar==2:
        (nys,nxs)=sfld.shape
        print(' 2D: nys=',nys,' nxs=',nxs)
        sfld2d=sfld
    elif ndim_svar==3:
        (nts,nys,nxs)=sfld.shape
        print(' time+2D: nts=',nts,' nys=',nys,' nxs=',nxs)
        sfld2d=np.squeeze(sfld,axis=0)

    out_title_fld=out_title_base+svar
    out_fname=out_fname_base+svar

    nm_svar=svar
    cs_cmap='nipy_spectral_r'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    n_rnd=0
    cmap_range='round'
 
    if svar=="AVG_NOAG_GROW":
        n_rnd=0


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

