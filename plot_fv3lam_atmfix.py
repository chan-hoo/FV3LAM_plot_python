###################################################################### CHJ #####
## Name		: plot_fv3lam_atmstat.py
## Language	: Python 3.7
## Usage	: Plot an output, atmos_static, for fv3 regional modeling
## Input files  : atmos_static.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/04/27: Chan-Hoo Jeon : Preliminary version
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
fnm_input='atmos_static.nc'

# Input orography file for reference grid
dnm_ref="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/expt_dirs/test_community/2020122700/INPUT/"
fnm_ref_oro='oro_data.nc'

# variables
vars_atm=["pk","bk","hyam","hybm","zsurf"]

# basic forms of title and file name: base+static field name
out_title_base='FV3LAM::ATMOS_STATIC::'
out_fname_base='fv3lam_out_atmfix_'

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'


# Main part (will be called at the end) =================== CHJ =====
def main():
# ========================================================= CHJ =====
    global atmfx,olon,olat
    global extent,c_lon,c_lat

    print(' ===== OUTPUT: atmos_static ================================')
    # open the data file
    fname=os.path.join(dnm_out,fnm_input)
    try: atmfx=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(atmfx)

    print(' ===== REFERENCE: Orography-grid (halo0) ==================')
    # open the data file
    fname=os.path.join(dnm_ref,fnm_ref_oro)
    try: oro=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(oro)
    # Extract grid info.
    olon=np.ma.masked_invalid(oro["geolon"].data)
    olat=np.ma.masked_invalid(oro["geolat"].data)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(olon)
    lon_max=np.max(olon)
    lat_min=np.min(olat)
    lat_max=np.max(olat)

    # Plot extent
    extent=[lon_min-5,lon_max+5,lat_min-5,lat_max+3]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])


    # Variables
    for svar in vars_atm:
        if svar=='zsurf':
            atmfx_plot_2d(svar)
        else:
            atmfx_plot_1d(svar)



# 1D Static field plot ======================================== CHJ =====
def atmfx_plot_1d(svar):
# ============================================================= CHJ =====

    print(' ===== '+svar+' ===== atmos. static data: 1D =============')
    # Extract data array
    sfld=np.ma.masked_invalid(atmfx[svar].data)
    (nxs,)=sfld.shape
    print(' nxs=',nxs)

    x1=np.linspace(1,nxs,nxs)

    if svar=='pk':
        x_label="Interface number"
        y_label="Pressure at interface"
        n_rnd=1
    elif svar=='bk':
        x_label="Interface number"
        y_label="Vertical coordinate sigma"
        n_rnd=1
    elif svar=='hyam':
        x_label="Interface number"
        y_label="Hybrid coefficient A"
        n_rnd=2
    elif svar=='hybm':
        x_label="Interface number"
        y_label="Hybrid coefficient B"
        n_rnd=1
    else:
        sys.exit('ERROR:: wrong svar ::'+svar)


    # Max and Min of the field
    fmax=np.max(sfld)
    fmin=np.min(sfld)
    print(' fld_max=',fmax)
    print(' flx_min=',fmin)
    cs_max=round(fmax,n_rnd)
    cs_min=round(fmin,n_rnd)
    print(' cs_max=',cs_max)
    print(' cs_min=',cs_min)

    mk_size=2
    mk_face='red'
    mk_edge_c='black'
    mk_edge_lw=0.3
    txt_fnt=7

    print(' SFC field=',y_label)

    out_title_fld=out_title_base+svar
    out_atmfx_fname=out_fname_base+svar

    plt.figure(figsize=(3,3))
    plt.plot(x1,sfld,'o',color='red',markersize=mk_size,
             markerfacecolor=mk_face,markeredgecolor=mk_edge_c,
             markeredgewidth=mk_edge_lw) 
    plt.xlabel(x_label,fontsize=txt_fnt)
    plt.ylabel(y_label,fontsize=txt_fnt)
    plt.xlim(1,nxs)
    plt.ylim(cs_min,cs_max)
    plt.title(out_title_fld,fontsize=txt_fnt+1)
    plt.tick_params(axis='both',which='major',labelsize=txt_fnt-1)
    plt.grid(True,linestyle=':',linewidth=0.5)

    # Output figure
    ndpi=300
    out_file(out_atmfx_fname,ndpi)



# 2D Static field plot ======================================== CHJ =====
def atmfx_plot_2d(svar):
# ============================================================= CHJ =====

    print(' ===== '+svar+' ===== atmos. static data: 2D =============')
    # Extract data array
    sfld=np.ma.masked_invalid(atmfx[svar].data)

    (nys,nxs)=sfld.shape
    print(' nys=',nys,' nxs=',nxs)

    (nyo,nxo)=olon.shape
    if nys==nyo and nxs==nxo:
        print(' Array sizes matched: oro/svar !!!!')
    else:
        sys.exit('ERROR:: mismatched array size with oro. ::'+svar)

    # Max and Min of the field
    fmax=np.max(sfld)
    fmin=np.min(sfld)
    print(' fld_max=',fmax)
    print(' flx_min=',fmin)
    cs_max=round(fmax,2)
    cs_min=round(fmin,2)
    print(' cs_max=',cs_max)
    print(' cs_min=',cs_min)

    nm_svar="Orography (m)"
    cs_cmap='terrain_r'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3


    print(' SFC field=',nm_svar)

    out_title_fld=out_title_base+svar
    out_atmfx_fname=out_fname_base+svar

    # Plot field
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_title_fld,fontsize=9)
    cs=ax.pcolormesh(olon,olat,sfld,cmap=cs_cmap,rasterized=True,
        vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(nm_svar,fontsize=8)

    # Output figure
    ndpi=300
    out_file(out_atmfx_fname,ndpi)

  

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

