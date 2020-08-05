###################################################################### CHJ #####
## Name		: plot_fv3lam_icsfc.py
## Language	: Python 3.7
## Usage	: Plot initial static fields for fv3 regional modeling
## Input files  : sfc_data.tile7.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/03/26: Chan-Hoo Jeon : Preliminary version
## V001: 2020/04/07: Chan-Hoo Jeon : Add refine ratio to output titles
## V002: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
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

# Path to Natural Earth Data-set for background plot
if machine=='hera':
    cartopy.config['data_dir']='/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/NaturalEarth'
elif machine=='orion':
    cartopy.config['data_dir']='/home/chjeon/tools/NaturalEarth'
else:
    sys.exit('ERROR: path to Natural Earth Data is not set !!!')

plt.switch_backend('agg')

# Global variables ======================================== CHJ =====
# ..... Case-dependent input :: should be changed case-by-case .....
# ******
# INPUT
# ******
# Path to the directory where the input NetCDF file is located.
dnm_data="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/C96/"
#dnm_data="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/brz_C768/"

# Path to the orography file
oro_path=dnm_data

# Input NetCDF file names
fnm_input='sfc_data.tile7.nc'

# Domain name:
domain='CONUS'

# Grid resolution ('C96'/'C768'):
res='C96'

# Grid type ('ESG'/'GFDL')
gtype='GFDL'

# GFDL grid-refinement ratio (for ESG grid, refine=0)
if gtype=='ESG':
    refine=0
elif gtype=='GFDL':
    refine=3

# Orography data file (w/ halo0) for coordinates
fnm_in_orog=res+'_oro_data.tile7.halo0.nc'

# Static fields
#vars_isfc=["alvsf","alvwf","alnsf","alnwf"]
vars_isfc=["slmsk","snoalb","zorl"]
#vars_isfc=["t2m","tg3","tisfc","tsea"]
#vars_isfc=["slope","stype","vtype"]
#vars_isfc=["stc"]

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

# basic forms of title and file name: base+static field name
if gtype=='ESG':
    out_title_base='IC(SFC)::'+domain+'(ESG)::'+res+'::'
    out_fname_base='fv3_isfc_'+domain+'_esg_'+res+'_'
elif gtype=='GFDL':
    out_title_base='IC(SFC)::'+domain+'(GFDL)::'+res+'(x'+str(refine)+')::'
    out_fname_base='fv3_isfc_'+domain+'_gfdl_'+res+'_'

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'

# Machine-specific mfdataset arguments
if machine=='hera':
    mfdt_kwargs={'parallel':False}
elif machine=='orion':
    mfdt_kwargs={'parallel':False,'combine':'by_coords'}
else:
    mfdt_kwargs={'parallel':False}


# Main part (will be called at the end) =================== CHJ =====
def main():
# ========================================================= CHJ =====
    global sttc

    # Open '_oro_data' file for coordinates
    oro_data()

    # open the data file
    fname=os.path.join(dnm_data,fnm_input)
    try: sttc=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)

    print(sttc)

    # Surface climatology fields
    for svar in vars_isfc:
        ic_sfc_plot(svar)

 
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


# Initial Static field plot ======================================== CHJ =====
def ic_sfc_plot(svar):
# ========================================================== CHJ =====

    print(' ===== '+svar+' ===== sfc data =======================')
    # Extract data array
    sfld=np.ma.masked_invalid(sttc[svar].data)

    if svar=='stc' or svar=='smc' or svar=='slc':
        (nts,nzs,nys,nxs)=sfld.shape
    else:    
        (nts,nys,nxs)=sfld.shape
        nzs=1

    print(' nts=',nts,' nzs=',nzs,' nys=',nys,' nxs=',nxs)
    (nyo,nxo)=glon.shape
    if nys==nyo and nxs==nxo:
        print(' Array sizes matched: oro/svar !!!!')
    else:  
        sys.exit('ERROR:: mismatched array size with oro. ::'+svar)

    if svar=='slmsk' or svar=='snoalb' or svar=='snoalb' or svar=='slope' or svar=='stype' or svar=='vtype':
        n_rnd=0
    else:
        n_rnd=2

    # Max and Min of the field
    fmax=np.max(sfld)
    fmin=np.min(sfld)
    print(' fld_max=',fmax)
    print(' flx_min=',fmin)
    cs_max=round(fmax,n_rnd)
    cs_min=round(fmin,n_rnd)
    print(' cs_max=',cs_max)
    print(' cs_min=',cs_min)

    # options for each var. =========================== CHJ =====
    # Default
    cs_cmap='jet'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3

    if svar=='alvsf':
        nm_svar='Visible black sky albedo' 
        cs_cmap='YlGnBu'
    elif svar=='alvwf':
        nm_svar='Visible white sky albedo'
        cs_cmap='YlOrRd'
    elif svar=='alnsf':
        nm_svar='Near-IR black sky albedo' 
        cs_cmap='YlGnBu'
    elif svar=='alnwf':
        nm_svar='Near-IR white sky albedo'
        cs_cmap='YlOrRd'
    elif svar=='slmsk':
        nm_svar='Sea(0),land(1),ice(2)'
        cs_min=0
        cs_max=2
        cs_cmap=plt.cm.get_cmap('Paired',3)
    elif svar=='snoalb':
        nm_svar='Maximum snow albedo'
        cs_map='nipy_spectral_r'
    elif svar=='zorl':
        nm_svar='Surface roughness'
        cs_cmap='gist_earth_r'
    elif svar=='t2m':
        nm_svar='2m temperature'
    elif svar=='tg3':        
        nm_svar='Deep soil temperature'
    elif svar=='tisfc':
        nm_svar='Sea-ice surface temperature'
    elif svar=='tsea':
        nm_svar='Sea-surface temperature'
    elif svar=='slope':
        nm_svar='Surface slope type'
        cs_cmap=plt.cm.get_cmap('Paired',10)
        cs_min=-0.5
        cs_max=9.5
    elif svar=='stype':
        nm_svar='Soil type'
        cs_cmap=plt.cm.get_cmap('tab20',17)
        cs_min=-0.5
        cs_max=16.5
    elif svar=='vtype':
        nm_svar='Vegetation type'
        cs_cmap=plt.cm.get_cmap('tab20',20)
        cs_min=-0.5
        cs_max=19.5
    elif svar=='stc':
        nm_svar='Soil temperature'
    else:
        sys.exit('ERROR: wrong svar or not set up yet: '+svar)
    # ================================================= CHJ =====

    print(' SFC field=',nm_svar) 

    for imn in range(nzs):
        if nzs==1:
        # Remove the time dimension in the array in case of 1
            sfld2d=np.squeeze(sfld,axis=0)
            out_title_fld=out_title_base+nm_svar+'('+svar+')'
            out_sfld_fname=out_fname_base+svar
        else:
            print('zaxis=',imn+1)
            sfld2d=sfld[0,imn,:,:]
            z_nm='z'+format(imn+1,'01d')
            out_title_fld=out_title_base+nm_svar+'('+svar+')'+'::'+z_nm
            out_sfld_fname=out_fname_base+svar+'_'+z_nm

        # Plot field
        fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
        ax.set_extent(extent, ccrs.PlateCarree())
        # Call background plot
        back_plot(ax)
        ax.set_title(out_title_fld,fontsize=9)
        cs=ax.pcolormesh(glon,glat,sfld2d,cmap=cs_cmap,rasterized=True,
              vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
        divider=make_axes_locatable(ax)
        ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        if svar=='vtype':
            cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext,ticks=np.arange(0,20,3))
        elif svar=='slmsk':
            cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext,ticks=[0,1,2])
        else:
            cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label(nm_svar,fontsize=8)

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

