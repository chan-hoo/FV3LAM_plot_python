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
## V003: 2021/03/05: Chan-Hoo Jeon : Simplify the script
## V004: 2021/08/16: Chan-Hoo Jeon : Remove oro file and add smc, slc
## V005: 2021/08/16: Chan-Hoo Jeon : Add original data plot
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
    #out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig/"
    out_fig_dir="/scratch2/NCEPDEV/naqfc/Chan-hoo.Jeon/test/fig/"
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
dnm_data="/scratch2/NCEPDEV/naqfc/Chan-hoo.Jeon/test/expt_dirs/test_CONUS25_netcdf/2021081500/INPUT/"

# Input NetCDF file names
fnm_input='sfc_data.nc'

# Path to the original data
dnm_org='/scratch2/NCEPDEV/naqfc/Chan-hoo.Jeon/test/expt_dirs/test_CONUS25_netcdf/2021081500/FV3GFS/for_ICS/'
# File name of the original data
fnm_org='gfs.t00z.sfcanl.nc'

domain_nm='RRFS_CONUS_25'

# Static fields
#vars_isfc=["alvsf","alvwf","alnsf","alnwf"]
#vars_isfc=["slmsk","snoalb"]
#vars_isfc=["slope","stype","vtype"]
vars_isfc=["stc","smc","slc"]
#vars_isfc=["smc"]

# basic forms of title and file name: base+static field name
out_title_base='FV3LAM::IC(SFC)::'+domain_nm+'::'
out_fname_base='fv3lam_isfc_'+domain_nm+'_'

# Flag for comparison with the original data ('on' or 'off')
flg_org='on'

if flg_org=='on':
    out_title_g_base='FV3GFS::IC(SFC)::'

# Resolution of background natural earth data ('50m' or '110m')
back_res='50m'


# Main part (will be called at the end) =================== CHJ =====
def main():
# ========================================================= CHJ =====
    global sttc
    global glon,glat,extent,c_lon,c_lat

    if flg_org=='on':
        global org_data,olon,olat,nxg,nyg

    print(' ===== input sfc data =============================')
    # open the data file
    fname=os.path.join(dnm_data,fnm_input)
    try: sttc=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)

    print(sttc)

    # Extract longitudes, and latitudes
    glon=np.ma.masked_invalid(sttc['geolon'].data)
    glat=np.ma.masked_invalid(sttc['geolat'].data)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(glon)
    lon_max=np.max(glon)
    lat_min=np.min(glat)
    lat_max=np.max(glat)

    # Plot extent
    extent=[lon_min-5,lon_max+5,lat_min-5,lat_max+3]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])

    # Original data plot
    if flg_org=='on':
        print(' ===== original sfc data =============================')
        # open the original data file
        fname=os.path.join(dnm_org,fnm_org)
        try: org_data=xr.open_mfdataset(fname,**mfdt_kwargs)
        except: raise Exception('Could NOT find the file',fname)
        print(org_data)

        # lon/lat of the original data
        olon=np.ma.masked_invalid(org_data['lon'].data)
        olat=np.ma.masked_invalid(org_data['lat'].data)
        (nyg,nxg)=olon.shape

    # Surface climatology fields
    for svar in vars_isfc:
        ic_sfc_plot(svar)



# Initial Static field plot ================================ CHJ =====
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

    svar_org=[]

    if svar=='alvsf':
        nm_svar='Visible black sky albedo' 
        cs_cmap='YlGnBu'
        svar_org=np.array(['alvsf'])
    elif svar=='alvwf':
        nm_svar='Visible white sky albedo'
        cs_cmap='YlOrRd'
        svar_org=np.array(['alvwf'])
    elif svar=='alnsf':
        nm_svar='Near-IR black sky albedo' 
        cs_cmap='YlGnBu'
        svar_org=np.array(['alnsf'])
    elif svar=='alnwf':
        nm_svar='Near-IR white sky albedo'
        cs_cmap='YlOrRd'
        svar_org=np.array(['alnwf'])
    elif svar=='slmsk':
        nm_svar='Sea(0),land(1),ice(2)'
        cs_min=0
        cs_max=2
        cs_cmap=plt.cm.get_cmap('Paired',3)
        svar_org=np.array(['land'])
    elif svar=='snoalb':
        nm_svar='Maximum snow albedo'
        cs_map='nipy_spectral_r'
        svar_org=np.array(['snoalb'])
    elif svar=='zorl':
        nm_svar='Surface roughness'
        cs_cmap='gist_earth_r'
        svar_org=np.array(['sfcr'])
    elif svar=='t2m':
        nm_svar='2m temperature'
        svar_org=np.array(['tmp2m'])
    elif svar=='tg3':        
        nm_svar='Deep soil temperature'
        svar_org=np.array(['tg3'])
    elif svar=='tisfc':
        nm_svar='Sea-ice surface temperature'
        svar_org=np.array(['tisfc'])
    elif svar=='tsea':
        nm_svar='Sea-surface temperature'
        svar_org=np.array(['tmpsfc'])
    elif svar=='slope':
        nm_svar='Surface slope type'
        cs_cmap=plt.cm.get_cmap('Paired',10)
        cs_min=-0.5
        cs_max=9.5
        svar_org=np.array(['sltyp'])
    elif svar=='stype':
        nm_svar='Soil type'
        cs_cmap=plt.cm.get_cmap('tab20',17)
        cs_min=-0.5
        cs_max=16.5
        svar_org=np.array(['sotyp'])
    elif svar=='vtype':
        nm_svar='Vegetation type'
        cs_cmap=plt.cm.get_cmap('tab20',20)
        cs_min=-0.5
        cs_max=19.5
        svar_org=np.array(['vtype'])
    elif svar=='stc':
        nm_svar='Soil temperature (K)'
        svar_org=np.array(['soilt1','soilt2','soilt3','soilt4'])
    elif svar=='slc':
        nm_svar='Unfrozen soil moisture (vol.frac)'
        cs_min=0
        cs_max=0.5
        cs_cmap='nipy_spectral'
        svar_org=np.array(['soill1','soill2','soill3','soill4'])
    elif svar=='smc':
        nm_svar='Soil moisture (vol.frac)'
        cs_min=0
        cs_max=0.5
        cs_cmap='nipy_spectral'
        svar_org=np.array(['soilw1','soilw2','soilw3','soilw4'])
    else:
        sys.exit('ERROR: wrong svar or not set up yet: '+svar)
    # ================================================= CHJ =====

    if flg_org=='on':
        if svar=='stc' or svar=='slc' or svar=='smc':
            sfld2d_org=np.empty((nzs,nyg,nxg))
            sfld2d_org[:]=np.NaN
            for iz in range(nzs):
                sfld_tmp=np.ma.masked_invalid(org_data[svar_org[iz]].data)
                sfld2d_org[iz,:,:]=sfld_tmp[0,:,:]
            sfld2d_org[sfld2d_org==0]='nan'
        else:
            sfld2d_org=np.ma.masked_invalid(org_data[svar_org[0]].data)        

    print(' SFC field=',nm_svar) 

    for imn in range(nzs):
        if nzs==1:
        # Remove the time dimension in the array in case of 1
            sfld2d=np.squeeze(sfld,axis=0)
            out_title_fld=out_title_base+svar
            out_sfld_fname=out_fname_base+svar
            out_title_fld_g=out_title_g_base+svar_org[0]
        else:
            print('zaxis=',imn+1)
            sfld2d=sfld[0,imn,:,:]
            z_nm='z'+format(imn+1,'01d')
            out_title_fld=out_title_base+svar+'::'+z_nm
            out_sfld_fname=out_fname_base+svar+'_'+z_nm
            out_title_fld_g=out_title_g_base+svar_org[imn]

        # Plot field
        if flg_org=='on':
            if domain_nm[:7]=='RRFS_NA':
                fig,(ax1,ax2)=plt.subplots(2,1,subplot_kw=dict(projection=ccrs.Orthographic(
                            central_longitude=-107,central_latitude=53)))
            else:
                fig,(ax1,ax2)=plt.subplots(2,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
                ax1.set_extent(extent, ccrs.PlateCarree())
                ax2.set_extent(extent, ccrs.PlateCarree())

            # Plot 1: original data
            back_plot(ax1)
            ax1.set_title(out_title_fld_g,fontsize=7)
            cs=ax1.pcolormesh(olon,olat,sfld2d_org[imn,:,:],cmap=cs_cmap,
                rasterized=True,vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
            # extend(pointed end): 'neither'|'both'|'min'|'max'  
            divider=make_axes_locatable(ax1)
            ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            if svar=='vtype':
                cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext,ticks=np.arange(0,20,3))
            elif svar=='slmsk':
                cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext,ticks=[0,1,2])
            else:
                cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
            cbar.ax.tick_params(labelsize=5)
            cbar.set_label(nm_svar,fontsize=5)

            # Plot 2: input data
            back_plot(ax2)
            ax2.set_title(out_title_fld,fontsize=7)
            cs=ax2.pcolormesh(glon,glat,sfld2d,cmap=cs_cmap,rasterized=True,
                vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
            # extend(pointed end): 'neither'|'both'|'min'|'max'  
            divider=make_axes_locatable(ax2)
            ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            if svar=='vtype':
                cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext,ticks=np.arange(0,20,3))
            elif svar=='slmsk':
                cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext,ticks=[0,1,2])
            else:
                cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
            cbar.ax.tick_params(labelsize=5)
            cbar.set_label(nm_svar,fontsize=5)

        else:
            if domain_nm[:7]=='RRFS_NA':
                fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Orthographic(
                            central_longitude=-107,central_latitude=53)))
            else:
                fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
                ax.set_extent(extent, ccrs.PlateCarree())
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

