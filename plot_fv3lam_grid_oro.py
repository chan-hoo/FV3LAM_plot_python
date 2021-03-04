###################################################################### CHJ #####
## Name		: plot_fv3lam_grid_oro.py
## Language	: Python 3.7
## Usage	: Plot regional FV3 super-grid and orograpy on the map
## Input files  : CXX_[grid,oro_data].tile7.haloX.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/02/24: Chan-Hoo Jeon : Preliminary version
## V001: 2020/02/27: Chan-Hoo Jeon : Separate Supergrid/Oro-coord.
## V002: 2020/03/05: Chan-Hoo Jeon : Add high-resolution earth data (50m)
## V003: 2020/04/02: Chan-Hoo Jeon : Add more plot options for orography 
## V004: 2020/04/07: Chan-Hoo Jeon : Add refine ratio to output titles
## V005: 2020/04/20: Chan-Hoo Jeon : Print out some ref. for model_configure/input_nml
## V006: 2020/04/23: Chan-Hoo Jeon : Add grid-uniformness plot
## V007: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
## V008: 2020/06/29: Chan-Hoo Jeon : Modify grid_dxy_plot to use area for GFDL/ESG
## V009: 2021/03/04: Chan-Hoo Jeon : Simplify the script
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
import math

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
# Path to the directory where the grid file is located.
dnm_data="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/expt_dirs/test_community/2020122700/INPUT/"

# path to the orography file
dnm_orog=dnm_data

# Grid file
fnm_in_grid='grid.tile7.halo4.nc'

# Orography file
fnm_in_orog='oro_data.tile7.halo4.nc'

# Grid point plot (every 'n_skip' rows/columns)
n_skip=30

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

# Bases of output title and file names
# Grid
out_grd_title='FV3LAM::grid'
out_grd_fname='fv3lam_grid'
# Orography
out_orog_title_base='FV3LAM::orography::'
out_orog_fname_base='fv3lam_orog'

# orography variables:
#orog_vars=["slmsk","land_frac","orog_raw","orog_filt","stddev","convexity",
#           "oa1","oa2","oa3","oa4","ol1","ol2","ol3","ol4",
#           "theta","gamma","sigma","elvmax"]
orog_vars=["orog_filt"]

# Color-map range option flag ('symmetry','round','real','fixed')
cmap_range_grd='round'

# Resolution of background natural earth data ('10m' or '50m' or '110m')
back_res='50m'

# Machine-specific mfdataset arguments
if machine=='hera':
    mfdt_kwargs={'parallel':False}
elif machine=='orion':
    mfdt_kwargs={'parallel':False,'combine':'by_coords'}
else:
    mfdt_kwargs={'parallel':False}


# Main part (will be called at the end) ============================== CHJ =====
def main():
# ==================================================================== CHJ =====
    orog_plot()
    grid_plot()


 
# Orography plot ===================================================== CHJ =====
def orog_plot():
# ==================================================================== CHJ =====
    global oro_x,oro_y,extent,c_lon,c_lat
    global npx,npy
     # open the orography file
    fname=os.path.join(dnm_orog,fnm_in_orog)
    try: oro=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)

    print(' ===== OROGRAPHY ================================')
    print(oro)

    # Extract longitudes, and latitudes
    oro_x=np.ma.masked_invalid(oro['geolon'].data)
    oro_y=np.ma.masked_invalid(oro['geolat'].data)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(oro_x)
    lon_max=np.max(oro_x)
    lat_min=np.min(oro_y)
    lat_max=np.max(oro_y)

    print(' ***** Ref. for lon1/lat1/lon2/lat2 in model_configure *****')
    print(' oro:lon-min(lon1)=',lon_min-360)
    print(' oro:lon-max(lon2)=',lon_max-360)
    print(' oro:lat-min(lat1)=',lat_min)
    print(' oro:lat-max(lat2)=',lat_max)

    print(' ***** npx/npy in input.nml/fv_core_nml *****')
    hcond=fnm_in_orog[-8:-3]
    npy,npx=oro_x.shape
    if hcond=='halo0':
        print(' npx=',npx+1)
        print(' npy=',npy+1)
    elif hcond=='halo4':
        print(' npx=',npx-7)
        print(' npy=',npy-7)
    else:
        sys.exit('ERROR: wrong fnm_in_base !!!!!')

    # Plot extent
    extent=[lon_min-5,lon_max+5,lat_min-5,lat_max+3]
    c_lon=np.mean(extent[:2])
    c_lat=np.mean(extent[2:])

    for ovar in orog_vars:
        oro_f=np.ma.masked_invalid(oro[ovar].data)
        orog_var_plot(oro_f,ovar)
        

   
# Orography plot ===================================================== CHJ =====
def orog_var_plot(oro_f,ovar):
# ==================================================================== CHJ =====

    # Orography: output title and file name
    out_orog_title=out_orog_title_base+ovar
    out_orog_fname=out_orog_fname_base+"_"+ovar


    print(' ===== '+ovar+' ===== orography =======================')

    # Max and Min of the field
    fmax=np.max(oro_f)
    fmin=np.min(oro_f)
    print(' fld_max=',fmax)
    print(' flx_min=',fmin)
    cs_max=round(fmax,2)
    cs_min=round(fmin,2)
    print(' cs_max=',cs_max)
    print(' cs_min=',cs_min)

    # options for each var. =========================== CHJ =====
    # Default
    cs_cmap='jet'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
    ndpi=300

    if ovar=='slmsk':
        nm_svar='Mask: sea(0),land(1)'
        cs_cmap=plt.cm.get_cmap('Paired',3)
        cs_min=0
        cs_max=2
    elif ovar=='land_frac':
        nm_svar='Land fraction' 
    elif ovar=='orog_raw' or ovar=='orog_filt':
        nm_svar='Orography'
        cs_cmap='terrain_r'
    elif ovar=='stddev':
        nm_svar='stddev'
        cs_cmap='gist_ncar_r'
    elif ovar=='convexity':
        nm_svar='Convexity'
        cs_cmap='gist_ncar_r'
    elif ovar=='oa1' or ovar=='oa2' or ovar=='oa3' or ovar=='oa4':
        nm_svar='Asymmetry parameter: '+ovar.upper()
        cs_cmap='seismic'
        cs_min=-0.6
        cs_max=0.6
        ndpi=150
    elif ovar=='ol1' or ovar=='ol2' or ovar=='ol3' or ovar=='ol4':
        nm_svar='Convexity parameter: '+ovar.upper()
        cs_cmap='nipy_spectral_r'
        cs_min=0
        cs_max=1
        ndpi=150
    elif ovar=='theta':
        nm_svar='theta'
        cs_cmap='seismic'
        cs_min=-100
        cs_max=100
    elif ovar=='gamma':
        nm_svar='gamma'
        cs_cmap='gist_ncar_r'
    elif ovar=='sigma':
        nm_svar='sigma'
        cs_cmap='gist_ncar_r'
    elif ovar=='elvmax':
        nm_svar='Maximum elevation'
        cs_cmap='gist_earth_r'
    else:
        sys.exit('ERROR: wrong ovar !!!!!')


    # Plot field
    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    # Call background plot
    back_plot(ax)
    ax.set_title(out_orog_title,fontsize=9)
    cs=ax.pcolormesh(oro_x,oro_y,oro_f,cmap=cs_cmap,rasterized=True,
              vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())
    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    if ovar=='slmsk':
        cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext,ticks=[0,1,2])
    else:
        cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(nm_svar,fontsize=8)


    # Output figure
    out_file(out_orog_fname,ndpi) 




# Grid plot ========================================================== CHJ =====
def grid_plot():
# ==================================================================== CHJ =====
    global grd_x,grd_y
     # open the grid file
    fname=os.path.join(dnm_data,fnm_in_grid)
    try: grd=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)

    print(' ===== GRID ====================================')
    print(grd)

    # Extract longitudes, latitudes, and others
    grd_x=np.ma.masked_invalid(grd['x'].data)
    grd_y=np.ma.masked_invalid(grd['y'].data)
    grd_area=np.ma.masked_invalid(grd['area'].data)

    # array size
    (nyp,nxp)=grd_x.shape
    print('grid array size=',grd_x.shape)
    (ny,nx)=grd_area.shape
    print('area array size=',grd_area.shape)

    # Hightest/Lowest longitudes and latitudes for text
    lon_min=np.min(grd_x)
    lon_max=np.max(grd_x)
    lat_min=np.min(grd_y)
    lat_max=np.min(grd_y)

    # Plot grids (super- and oro-)
    grid_super_oro_plot(lon_min,lat_min)

    # Plot grid sizes
    grid_dxy_plot(grd_area,lon_min,lat_min) 

    # Plot boundary
    grid_bndr_plot(nxp,nyp)

    # Plot around four corners
    c_nms=['R1C1','RxC1','R1Cx','RxCx']
    for c_nm in c_nms:
        grid_crnr_plot(c_nm)
    


# Grid plot: super/oro =============================================== CHJ =====
def grid_super_oro_plot(lon_min,lat_min):
# ==================================================================== CHJ =====
    print(' ===== super-/oro- ===== GRID ====================================')

    # grid points (every 'n_skip' rows/columns from 2nd row/col)
    grdx_slc=grd_x[1:-1:n_skip,1:-1:n_skip]
    grdy_slc=grd_y[1:-1:n_skip,1:-1:n_skip]

    # grid points (every 'n_skip' rows/columns)
    orox_slc=oro_x[::n_skip,::n_skip]
    oroy_slc=oro_y[::n_skip,::n_skip]

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.set_title(out_grd_title, fontsize=9)

    # Call background plot
    back_plot(ax)

    # Scatter plot (zorder: lowest-plot on bottom, highest-plot on top)
    sp_scale=2
    s1=ax.scatter(grdx_slc,grdy_slc,transform=ccrs.PlateCarree(),marker='o',facecolors="None",edgecolors='red',linewidth=0.3,s=sp_scale,zorder=3)

    s2=ax.scatter(orox_slc,oroy_slc,transform=ccrs.PlateCarree(),marker='*',facecolors='green',edgecolors='green',linewidth=0.3,s=sp_scale,zorder=4)

    ref_txt='Super-: every '+str(n_skip)+' (i,j)s from (2,2), Oro-: every '+str(n_skip)+' (i,j)s'
    plt.text(lon_min-3,lat_min-3,ref_txt,transform=ccrs.Geodetic(),fontsize=8)
    plt.legend((s1,s2),('super-grid','oro-grid'),scatterpoints=1,loc='lower right',ncol=1,fontsize=8) 

    # Output figure
    ndpi=300
    out_file(out_grd_fname,ndpi)



# Grid plot: dx/dy =================================================== CHJ =====
def grid_dxy_plot(grd_area,lon_min,lat_min):
# ==================================================================== CHJ =====
    print(' ===== cell size ===== GRID =========================================')

    cell_area=np.zeros((npy,npx))
    cell_dxy=np.zeros((npy,npx))

    for iy in range(npy):
        for jx in range(npx):
            iy2=2*iy
            iy2p1=2*iy+1
            jx2=2*jx
            jx2p1=2*jx+1

            cell_area[iy,jx]=grd_area[iy2,jx2]+grd_area[iy2,jx2p1]+grd_area[iy2p1,jx2]+grd_area[iy2p1,jx2p1]
            cell_dxy[iy,jx]=math.sqrt(cell_area[iy,jx])/1000

    # Max and Min of the field
    fmax=np.max(cell_dxy)
    fmin=np.min(cell_dxy)
    favg=np.average(cell_dxy)
    print(' fld_max=',fmax)
    print(' flx_min=',fmin)
    print(' flx_avg=',favg)

    # Set the colormap range
#    cmap_range_grd='round'
    n_rnd=2
    print(' cmap range=',cmap_range_grd)
    if cmap_range_grd=='symmetry':
        tmp_cmp=max(abs(fmax),abs(fmin))
        cs_min=round(-tmp_cmp,n_rnd)
        cs_max=round(tmp_cmp,n_rnd)
        cs_avg=round(favg,n_rnd)
    elif cmap_range_grd=='round':
        cs_min=round(fmin,n_rnd)
        cs_max=round(fmax,n_rnd)
        cs_avg=round(favg,n_rnd)
    elif cmap_range_grd=='real':
        cs_min=fmin
        cs_max=fmax
        cs_avg=favg
    elif cmap_range_grd=='fixed':
        cs_min=2.7
        cs_max=3.2
        cs_avg=round(favg,n_rnd)
    else:
        sys.exit('ERROR: wrong colormap-range flag !!!')

    print(' cs_min=',cs_min)
    print(' cs_max=',cs_max)
    print(' cs_avg=',cs_avg)


    nm_svar='Cell size (km)'
    cs_cmap='nipy_spectral_r'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())

    out_grd_dx_title=out_grd_title+'::Cell Size'
    ax.set_title(out_grd_dx_title, fontsize=9)

    # Call background plot
    back_plot(ax)

    cs=ax.pcolormesh(oro_x,oro_y,cell_dxy,cmap=cs_cmap,rasterized=True,
              vmin=cs_min,vmax=cs_max,transform=ccrs.PlateCarree())

    ref_txt='Max='+str(round(fmax,2))+', Min='+str(round(fmin,2))+', Avg='+str(round(favg,2))
    plt.text(lon_min-2,lat_min-2,ref_txt,horizontalalignment='left',
             transform=ccrs.Geodetic(),fontsize=9)

    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(nm_svar,fontsize=8)


    # Output figure
    out_grd_dx_fname=out_grd_fname+'_dxy'
    ndpi=300
    out_file(out_grd_dx_fname,ndpi)



# Grid boundary plot ================================================= CHJ =====
def grid_bndr_plot(nxp,nyp):
# ==================================================================== CHJ =====
    print(' ===== boundary ===== GRID ====================================')

    # Boundary: 1C (1st column of the array) 
    grd_B1C_lon=grd_x[:,0]
    grd_B1C_lat=grd_y[:,0]
    # Boundary: 1R (1st row of the array)
    grd_B1R_lon=grd_x[0,:]
    grd_B1R_lat=grd_y[0,:]
    # Boundary: xC (last column of the array)
    grd_BxC_lon=grd_x[:,-1]
    grd_BxC_lat=grd_y[:,-1]
    # Boundary: xR (last row of the array)
    grd_BxR_lon=grd_x[-1,:]
    grd_BxR_lat=grd_y[-1,:]

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(c_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
 
    out_grd_bndr_title=out_grd_title+'::Boundary'
    ax.set_title(out_grd_bndr_title, fontsize=9)

    # Call background plot
    back_plot(ax)

    # Scatter plot (zorder: lowest-plot on bottom, highest-plot on top)
    sp_scale=0.1
    ax.scatter(grd_B1C_lon,grd_B1C_lat,transform=ccrs.PlateCarree(),c='red',s=sp_scale,label='B1C',zorder=3)
    ax.scatter(grd_B1R_lon,grd_B1R_lat,transform=ccrs.PlateCarree(),c='blue',s=sp_scale,label='B1R',zorder=3)
    ax.scatter(grd_BxC_lon,grd_BxC_lat,transform=ccrs.PlateCarree(),c='purple',s=sp_scale,label='BxC',zorder=3)
    ax.scatter(grd_BxR_lon,grd_BxR_lat,transform=ccrs.PlateCarree(),c='green',s=sp_scale,label='BxR',zorder=3)
    # Add text to each boundary
    tsize=9
    ntxt=int(nyp/2)
    txt_x=grd_B1C_lon[ntxt]+1
    txt_y=grd_B1C_lat[ntxt]-1
    ax.text(txt_x,txt_y,'B1C',color='red',fontsize=tsize,transform=ccrs.PlateCarree())
    txt_x=grd_BxC_lon[ntxt]+1
    txt_y=grd_BxC_lat[ntxt]+1
    ax.text(txt_x,txt_y,'BxC',color='purple',fontsize=tsize,transform=ccrs.PlateCarree())
    ntxt=int(nxp/2)
    txt_x=grd_B1R_lon[ntxt]
    txt_y=grd_B1R_lat[ntxt]+1
    ax.text(txt_x,txt_y,'B1R',color='blue',fontsize=tsize,transform=ccrs.PlateCarree())
    txt_x=grd_BxR_lon[ntxt]
    txt_y=grd_BxR_lat[ntxt]+1
    ax.text(txt_x,txt_y,'BxR',color='green',fontsize=tsize,transform=ccrs.PlateCarree())
 
    # File name
    out_grd_bndr_fname=out_grd_fname+'_bndr'
    # Output figure
    ndpi=300
    out_file(out_grd_bndr_fname,ndpi)



# Grid corner plot =================================================== CHJ =====
def grid_crnr_plot(c_nm):
# ==================================================================== CHJ =====
    print(' ===== corner:',c_nm,'===== GRID ====================================')
    
    # Corner:
    N_crn=30
    N_crn_oro=int(N_crn/2)
    if c_nm=='R1C1':
        grdx_crnr=grd_x[0:N_crn,0:N_crn]
        grdy_crnr=grd_y[0:N_crn,0:N_crn]
        orox_crnr=oro_x[0:N_crn_oro,0:N_crn_oro]
        oroy_crnr=oro_y[0:N_crn_oro,0:N_crn_oro]
    elif c_nm=='RxC1':
        grdx_crnr=grd_x[-N_crn:,0:N_crn]
        grdy_crnr=grd_y[-N_crn:,0:N_crn]
        orox_crnr=oro_x[-N_crn_oro:,0:N_crn_oro]
        oroy_crnr=oro_y[-N_crn_oro:,0:N_crn_oro]
    elif c_nm=='R1Cx':
        grdx_crnr=grd_x[0:N_crn,-N_crn:]
        grdy_crnr=grd_y[0:N_crn,-N_crn:]
        orox_crnr=oro_x[0:N_crn_oro,-N_crn_oro:]
        oroy_crnr=oro_y[0:N_crn_oro,-N_crn_oro:]
    elif c_nm=='RxCx':
        grdx_crnr=grd_x[-N_crn:,-N_crn:]
        grdy_crnr=grd_y[-N_crn:,-N_crn:]
        orox_crnr=oro_x[-N_crn_oro:,-N_crn_oro:]
        oroy_crnr=oro_y[-N_crn_oro:,-N_crn_oro:]
    else:
        sys.exit('ERROR: wrong c_nm !!!')

    #print(grdx_crnr.shape)
    #print(orox_crnr.shape)

    # Highest and lowest longitudes and latitudes for plot extent
    lon_min=np.min(grdx_crnr)
    lon_max=np.max(grdx_crnr)
    lat_min=np.min(grdy_crnr)
    lat_max=np.max(grdy_crnr)

    # Plot extent
    extent=[lon_min,lon_max,lat_min,lat_max]
    central_lon=np.mean(extent[:2])
    central_lat=np.mean(extent[2:])

    fig,ax=plt.subplots(1,1,subplot_kw=dict(projection=ccrs.Robinson(central_lon)))
    ax.set_extent(extent, ccrs.PlateCarree())
 
    out_grd_crnr_title=out_grd_title+'::Corner::'+c_nm
    ax.set_title(out_grd_crnr_title, fontsize=9)

    # Call background plot
    back_plot(ax)

    # Scatter plot (zorder: lowest-plot on bottom, highest-plot on top)
    sp_scale=5
    s1=ax.scatter(grdx_crnr,grdy_crnr,transform=ccrs.PlateCarree(),marker='o',facecolors="None",edgecolors='red',linewidth=0.5,s=sp_scale,zorder=3)

    s2=ax.scatter(orox_crnr,oroy_crnr,transform=ccrs.PlateCarree(),marker='*',facecolors='green',edgecolors='green',linewidth=0.5,s=sp_scale,zorder=4)

    ref_txt=str(N_crn)+' rows/columns from the corner '+c_nm
    plt.text(lon_min,lat_min,ref_txt,transform=ccrs.Geodetic(),fontsize=8,color='purple')
    plt.legend((s1,s2),('super-grid','oro-grid'),scatterpoints=1,loc='upper right',ncol=1,fontsize=8) 

    out_grd_crnr_fname=out_grd_fname+'_crnr_'+c_nm

    # Output figure
    ndpi=300
    out_file(out_grd_crnr_fname,ndpi)



# Background plot ========================================== CHJ =====
def back_plot(ax):
# ========================================================== CHJ =====
    fline_wd=0.5  # line width
    falpha=0.3 # transparency

    # natural_earth
#    land=cfeature.NaturalEarthFeature('physical','land',back_res,
#                      edgecolor='face',facecolor=cfeature.COLORS['land'],
#          	      alpha=falpha)
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

