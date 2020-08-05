###################################################################### CHJ #####
## Name		: plot_fv3sar_icbc.py
## Language	: Python 3.7
## Usage	: Plot time-dependent IC/LBC fields for fv3 regional modeling
## Input files  : gfs_bndy.tile7.XXX.nc and gfs_data.tile7.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/03/06: Chan-Hoo Jeon : Preliminary version
## V001: 2020/03/10: Chan-Hoo Jeon : Add velocities to plot options
## V002: 2020/04/03: Chan-Hoo Jeon : Read GFS data and index once
## V003: 2020/04/07: Chan-Hoo Jeon : Add refine ratio to output titles
## V004: 2020/04/21: Chan-Hoo Jeon : Print out max/min every time step
## V005: 2020/05/26: Chan-Hoo Jeon : Modified to include blending layers
## V006: 2020/06/22: Chan-Hoo Jeon : Add opt. for machine-specific arguments
###################################################################### CHJ #####

import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.lines as mlines
import numpy as np
import xarray as xr
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.switch_backend('agg')

# Global variables =================================================== CHJ =====
# ..... Case-dependent input :: should be changed case-by-case .....
# ******
# INPUT
# ******
# Path to the directory where the input NetCDF files are located.
dnm_data="/scratch2/NCEPDEV/stmp1/Chan-hoo.Jeon/run_C96/INPUT_blend10/"
#dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/regional_workflow/fix/fix_sar/pr/"
#dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/TB_work/work_FV3_regional_C96_2019101500_blend/INPUT_blend10/"

# File names of initial and boundary conditions 
fnm_in_bndr_base='gfs_bndy.tile7.'
#bndr_step=["000","003","006","009","012"]   # time steps of input netcdf files (3 digits)
bndr_step=["000"]

# Variables 
#vars_bc=["ps","w","zh","t","sphum","liq_wat","o3mr","ice_wat",
#         "rainwat","snowwat","graupel","u_w","v_w","u_s","v_s"]
#vars_bc=["sphum","liq_wat","o3mr","ice_wat"]
#vars_bc=["rainwat","snowwat","graupel"]
#vars_bc=["u_s","v_s","u_w","v_w"]
vars_bc=["t","u_s"]

# number of the target vertical level (only for 3-D fields)
ilvl=1

# Domain name:
domain='CONUS'

# Grid resolution:
res='C96'

# Grid type ('ESG'/'GFDL')
gtype='GFDL'

# GFDL grid-refinement ratio (for ESG grid, refine=0)
if gtype=='ESG':
    refine=0
elif gtype=='GFDL':
    refine=3

# HPC machine ('hera','orion')
machine='hera'

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

# basic forms of title and file name
if gtype=='ESG':
    out_title_base='IC/LBC::'+domain+'(ESG)::'+res
    out_fname_base='fv3_icbc_'+domain+'_esg_'+res+'_'
elif gtype=='GFDL':
    out_title_base='IC/LBC::'+domain+'(GFDL)::'+res+'(x'+str(refine)+')'
    out_fname_base='fv3_icbc_'+domain+'_gfdl_'+res+'_'

# Number of additional boundary arrays (halo)
n_halo=4

# Number of blending layers (nrows_blend)
n_blend=10

# Total number of halos
n_halo_all=n_halo+n_blend

bc_lvl=format(ilvl,'03d')
sblend=format(n_blend,'02d') 

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

    print(' ===== Read GFS DATA ================================================')
    # Initial GFS Data of the main domain
    fnm_input_var='gfs_data.tile7.nc' 
    fname=os.path.join(dnm_data,fnm_input_var)
    try: tvar=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(tvar)

    print(' ===== IC/BC :: i,j indexes =========================================')
    fnm_input_var=fnm_in_bndr_base+'000.nc'
    # open the data file
    fname=os.path.join(dnm_data,fnm_input_var)
    try: ibcf=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(ibcf)
 
    # variables
    for svar in vars_bc:    # variables
        bndr_index(svar,ibcf)
        # Plot initial gfs data inside the regular domain for comparison
        plot_gfs(tvar,svar)

        for ibc in bndr_step:   # time steps
            # main plot
            plot_bc(svar,ibc)



# Read index of boundary fields ====================================== CHJ =====
def bndr_index(svar,ibcf):
# ==================================================================== CHJ =====
    global i_bottom,j_bottom,i_top,j_top,i_left,j_left,i_right,j_right
    global i_b,i_b_sym,i_bx_sym,j_l_sym
    global i_b_xtk,i_b_xlb,j_b_ytk,j_t_ytk,j_l_ytk,j_l_ylb
    global nd_xj,nd_yi,i_bx_xtk,i_bx_xlb

    print(' ===== BC Field Index :: '+svar+' ==================================')

    if svar=='u_w' or svar=='v_w':
        # top & bottom: (halo,domain i+2*halo+1) : i(x-direction), j(y-direction)
        i_bottom=np.ma.masked_invalid(ibcf['i_w_bottom'].data)
        j_bottom=np.ma.masked_invalid(ibcf['j_w_bottom'].data)
        i_top=np.ma.masked_invalid(ibcf['i_w_top'].data)
        j_top=np.ma.masked_invalid(ibcf['j_w_top'].data)
        # left & right: (domain j,halo)
        i_left=np.ma.masked_invalid(ibcf['i_w_left'].data)
        j_left=np.ma.masked_invalid(ibcf['j_w_left'].data)
        i_right=np.ma.masked_invalid(ibcf['i_w_right'].data)
        j_right=np.ma.masked_invalid(ibcf['j_w_right'].data)
        n_jb=j_bottom.shape[0]
        if n_jb!=n_halo_all:
            sys.exit('ERROR: n_halo+n_blend is not the same as BC files !!!')
    elif svar=='u_s' or svar=='v_s':
        # top & bottom: (halo,domain i+2*halo) : i(x-direction), j(y-direction)
        i_bottom=np.ma.masked_invalid(ibcf['i_s_bottom'].data)
        j_bottom=np.ma.masked_invalid(ibcf['j_s_bottom'].data)
        i_top=np.ma.masked_invalid(ibcf['i_s_top'].data)
        j_top=np.ma.masked_invalid(ibcf['j_s_top'].data)
        # left & right: (domain j,halo+1)
        i_left=np.ma.masked_invalid(ibcf['i_s_left'].data)
        j_left=np.ma.masked_invalid(ibcf['j_s_left'].data)
        i_right=np.ma.masked_invalid(ibcf['i_s_right'].data)
        j_right=np.ma.masked_invalid(ibcf['j_s_right'].data)
        n_jb=i_right.shape[0]
        if n_jb!=n_halo_all:
            sys.exit('ERROR: n_halo+n_blend is not the same as BC files !!!')
    else:
        # top & bottom: (halo,domain i+2*halo) : i(x-direction), j(y-direction)
        i_bottom=np.ma.masked_invalid(ibcf['i_bottom'].data)
        j_bottom=np.ma.masked_invalid(ibcf['j_bottom'].data)
        i_top=np.ma.masked_invalid(ibcf['i_top'].data)
        j_top=np.ma.masked_invalid(ibcf['j_top'].data)
        # left & right: (domain j,halo)
        i_left=np.ma.masked_invalid(ibcf['i_left'].data)
        j_left=np.ma.masked_invalid(ibcf['j_left'].data)
        i_right=np.ma.masked_invalid(ibcf['i_right'].data)
        j_right=np.ma.masked_invalid(ibcf['j_right'].data)
        n_jb=j_bottom.shape[0]
        if n_jb!=n_halo_all:
            sys.exit('ERROR: n_halo+n_blend is not the same as BC files !!!')


    # new coordinates for expanding axis-scales from mid
    # nd_xj: x-dimensional length of the real domain
    nd_xj=len(i_bottom)-2*n_halo
    nd_yi=len(j_left)

    print(' Domain size (nlat,nlon)=(',nd_yi,',',nd_xj,')')

    # bottom
    i_b_sym=axis_contract(i_bottom)
    nib=len(i_b_sym)
    i_b_xtk,i_b_xlb=tick_label(nib,i_b_sym,i_bottom)

    nskp=int(len(j_bottom)/4)
    if nskp<1:
        nskp=1
    j_b_ytk=j_bottom[::nskp]

    # top
    nskp=int(len(j_top)/4)
    if nskp<1:
        nskp=1
    j_t_ytk=j_top[::nskp]

    # left
    j_l_sym=axis_contract(j_left)
    nib=len(j_l_sym)
    j_l_ytk,j_l_ylb=tick_label(nib,j_l_sym,j_left)   

    # x-axis only for domain
    i_b=i_bottom[n_halo:-n_halo]
#    print(i_b)
    i_bx_sym=axis_contract(i_b)
    nib=len(i_b)
    i_bx_xtk,i_bx_xlb=tick_label(nib,i_bx_sym,i_b)

#    print('i_bottom::',i_bottom)
#    print('j_bottom::',j_bottom)
#    print('i_top::',i_top)
#    print('j_top::',j_top)
#    print('i_left::',i_left)
#    print('j_left::',j_left)
#    print('i_right::',i_right)
#    print('j_right::',j_right)
    

 
# Plot each variable ================================================= CHJ =====
def plot_bc(svar,ibc):
# ==================================================================== CHJ =====
    print(' ===== BC Field:: '+svar+' :: ==== at '+ibc+' ==============')

    fnm_input_var=fnm_in_bndr_base+ibc+'.nc'
    # open the data file
    fname=os.path.join(dnm_data,fnm_input_var)
    print(fname)
    try: ibcf=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)

    # Extract data of svar at ibcf
    bc_bottom=np.ma.masked_invalid(ibcf[svar+'_bottom'].data)
    bc_top=np.ma.masked_invalid(ibcf[svar+'_top'].data)
    bc_right=np.ma.masked_invalid(ibcf[svar+'_right'].data)
    bc_left=np.ma.masked_invalid(ibcf[svar+'_left'].data)

    bc1d_bottom=bc_bottom.ravel()
    bc1d_top=bc_top.ravel()
    bc1d_right=bc_right.ravel()
    bc1d_left=bc_left.ravel()
    bc1d_all=np.concatenate([bc1d_bottom,bc1d_top,bc1d_right,bc1d_left])
   
    ndim_bc=bc_bottom.ndim

    if ndim_bc==2: # 2d array
        bc2d_bottom=bc_bottom
        bc2d_top=bc_top
        bc2d_right=bc_right
        bc2d_left=bc_left
        out_bc_fname=out_fname_base+svar+'_t'+ibc+'_B'+sblend
    elif ndim_bc==3: # 3d array (vertical level,:,:)
        bc2d_bottom=bc_bottom[ilvl-1,:,:]
        bc2d_top=bc_top[ilvl-1,:,:]
        bc2d_right=bc_right[ilvl-1,:,:]
        bc2d_left=bc_left[ilvl-1,:,:]
        out_bc_fname=out_fname_base+svar+'_L'+bc_lvl+'_t'+ibc+'_B'+sblend


    var_min=[]
    var_max=[]
    # Max/min: overall
    var_min.append(np.min(bc2d_bottom))
    var_min.append(np.min(bc2d_top))
    var_min.append(np.min(bc2d_left))
    var_min.append(np.min(bc2d_right))
    var_min_all=np.min(var_min)
    var_max.append(np.max(bc2d_bottom))
    var_max.append(np.max(bc2d_top))
    var_max.append(np.max(bc2d_left))
    var_max.append(np.max(bc2d_right))
    var_max_all=np.max(var_max)

    print("MIN=",var_min_all)
    print("MAX=",var_max_all)

    # Plot field
    fig=plt.figure(figsize=(3,3))
    grid=plt.GridSpec(4,4,wspace=0.1,hspace=0.1)
    # center
    ax_c=fig.add_subplot(grid[1:-1,1:-1])
    # left
    ax_l=fig.add_subplot(grid[1:-1,0])
    # right
    ax_r=fig.add_subplot(grid[1:-1,3])
    # top
    ax_t=fig.add_subplot(grid[0,:])
    # bottom
    ax_b=fig.add_subplot(grid[3,:])
    
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=2.5
    
    # bottom
    xj_b,yi_b=np.meshgrid(i_b_sym,j_bottom)
    ax_b.pcolormesh(xj_b,yi_b,bc2d_bottom,cmap=cs_cmap,vmin=cs_min,vmax=cs_max)
    ax_b.tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz)
    ax_b.set_xticks(i_b_xtk)
    ax_b.set_xticklabels(i_b_xlb)
    ax_b.set_yticks(j_bottom)
    ax_b.set_yticklabels(j_bottom)

    if n_blend>0:    
        x1=(i_b_sym[n_halo-1]+i_b_sym[n_halo])*0.5
        x2=(i_b_sym[-n_halo-1]+i_b_sym[-n_halo])*0.5
        y1=0.5
        ln1=mlines.Line2D([x1,x2],[y1,y1],color='w',lw=0.7,ls='--')
        ax_b.add_artist(ln1) 

        x1=(i_b_sym[n_halo-1]+i_b_sym[n_halo])*0.5
        y1=(j_bottom[n_halo-1]+j_bottom[n_halo])*0.5
        y2=j_bottom[-1]
        ln1=mlines.Line2D([x1,x1],[y1,y2],color='w',lw=0.7,ls='--')
        ax_b.add_artist(ln1) 

        x1=(i_b_sym[-n_halo-1]+i_b_sym[-n_halo])*0.5
        y1=(j_bottom[n_halo-1]+j_bottom[n_halo])*0.5
        y2=j_bottom[-1]
        ln1=mlines.Line2D([x1,x1],[y1,y2],color='w',lw=0.7,ls='--')
        ax_b.add_artist(ln1) 

    # top
    xj_t,yi_t=np.meshgrid(i_b_sym,j_top)
    ax_t.pcolormesh(xj_t,yi_t,bc2d_top,cmap=cs_cmap,vmin=cs_min,vmax=cs_max)
    ax_t.tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz,
                     bottom=False,labelbottom=False,top=True,labeltop=False)
    ax_t.set_xticks(i_b_xtk)
    ax_t.set_xticklabels(i_b_xlb)
    ax_t.set_yticks(j_top)
    ax_t.set_yticklabels(j_top)

    if n_blend>0:
        x1=(i_b_sym[n_halo-1]+i_b_sym[n_halo])*0.5
        x2=(i_b_sym[-n_halo-1]+i_b_sym[-n_halo])*0.5
        y1=j_top[-n_halo-1]+0.5
        ln1=mlines.Line2D([x1,x2],[y1,y1],color='w',lw=0.7,ls='--')
        ax_t.add_artist(ln1) 

        x1=(i_b_sym[n_halo-1]+i_b_sym[n_halo])*0.5
        y2=(j_top[-n_halo-1]+j_top[-n_halo])*0.5
        y1=j_top[0]
        ln1=mlines.Line2D([x1,x1],[y1,y2],color='w',lw=0.7,ls='--')
        ax_t.add_artist(ln1) 

        x1=(i_b_sym[-n_halo-1]+i_b_sym[-n_halo])*0.5
        y2=(j_top[-n_halo-1]+j_top[-n_halo])*0.5
        y1=j_top[0]
        ln1=mlines.Line2D([x1,x1],[y1,y2],color='w',lw=0.7,ls='--')
        ax_t.add_artist(ln1) 

    # left
    xj_l,yi_l=np.meshgrid(i_left,j_l_sym)
    ax_l.pcolormesh(xj_l,yi_l,bc2d_left,cmap=cs_cmap,vmin=cs_min,vmax=cs_max)
    ax_l.tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz,
                     labelbottom=False)
    ax_l.set_yticks(j_l_ytk)
    ax_l.set_yticklabels(j_l_ylb)
    ax_l.set_xticks(i_left)

    if n_blend>0:
        x1=i_left[n_halo-1]+0.5
        y1=j_l_sym[0]
        y2=j_l_sym[-1]
        ln1=mlines.Line2D([x1,x1],[y1,y2],color='w',lw=0.7,ls='--')
        ax_l.add_artist(ln1) 

    # right
    xj_r,yi_r=np.meshgrid(i_right,j_l_sym)
    ax_r.pcolormesh(xj_r,yi_r,bc2d_right,cmap=cs_cmap,vmin=cs_min,vmax=cs_max)
    ax_r.tick_params(direction='out',length=tick_ln,width=tick_wd,
                     left=False,labelleft=False,labelbottom=False,right=True)
    ax_r.set_yticks(j_l_ytk)
    ax_r.set_yticklabels(j_l_ylb)
    ax_r.set_xticks(i_right)

    if n_blend>0:
        x1=i_right[-n_halo]-0.5
        y1=j_l_sym[0]
        y2=j_l_sym[-1]
        ln1=mlines.Line2D([x1,x1],[y1,y2],color='w',lw=0.7,ls='--')
        ax_r.add_artist(ln1) 

    # center
    ax_c.axis([0,10,0,10])
    ax_c.axis('off')

    ax_c.text(5,7,out_title_base,fontsize=5,color='black',
              horizontalalignment='center',verticalalignment='center')
    ax_c.text(5,6,'Variable='+svar,fontsize=5,color='black',
              horizontalalignment='center',verticalalignment='center')
    ax_c.text(5,5,'Layer #='+bc_lvl,fontsize=5,color='black',
              horizontalalignment='center',verticalalignment='center')
    ax_c.text(5,4,'Time='+ibc,fontsize=5,color='black',
              horizontalalignment='center',verticalalignment='center')
    ax_c.text(5,2,'N_blend='+str(n_blend),fontsize=4,color='black',
              horizontalalignment='center',verticalalignment='center')   

    ax_c.text(5,9.95,'Top',fontsize=tlb_sz,horizontalalignment='center',
              verticalalignment='center')
    ax_c.text(5,0.05,'Bottom',fontsize=tlb_sz,horizontalalignment='center',
              verticalalignment='center')
    ax_c.text(0.05,5,'Left',fontsize=tlb_sz,horizontalalignment='center',
              verticalalignment='center',rotation='vertical')
    ax_c.text(9.95,5,'Right',fontsize=tlb_sz,horizontalalignment='center',
              verticalalignment='center',rotation='vertical')

    # Output figure
    out_file(out_bc_fname)



# Plot gfs data ====================================================== CHJ =====
def plot_gfs(tvar,svar):
# ==================================================================== CHJ =====
    global cs_cmap,cs_min,cs_max

    print(' ===== GFS Data::'+svar+' ===============================')
  
    sgfs=np.ma.masked_invalid(tvar[svar].data)
    ndim_sgfs=sgfs.ndim

   
    if ndim_sgfs==2: # 2D array
        nlvl=1
        nyi,nxj=sgfs.shape
        sgfs2d=sgfs
        out_gfs_fname=out_fname_base+svar+'_gfs'
        out_title_gfs_o=out_title_base+'::GFS::'+svar+'t000::Regular'
        out_title_gfs_n=out_title_base+'::GFS::'+svar+'t000::Contracted'
    elif ndim_sgfs==3: # 3D array
        nlvl,nyi,nxj=sgfs.shape
        sgfs2d=sgfs[ilvl-1,:,:]
        out_gfs_fname=out_fname_base+svar+'_L'+bc_lvl+'_gfs'
        out_title_gfs_o=out_title_base+'::GFS::'+svar+'::L'+bc_lvl+'t000::Regular'
        out_title_gfs_n=out_title_base+'::GFS::'+svar+'::L'+bc_lvl+'t000::Contracted'
    else:
        sys.exit('ERROR: Wrong dimension of GFS !!!')


    fmin=np.min(sgfs2d)
    fmax=np.max(sgfs2d)

    cmap_range='round'
    n_rnd=2

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

    print(' cs_min=',cs_min)
    print(' cs_max=',cs_max)


    # Automate a range of color bars
    cs_cmap='jet'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=3
 
    # Plot field
    if svar=='u_s' or svar=='v_s':
        jlext_len=len(j_left)+2
        j_l_ext=np.arange(1,jlext_len+1,1)
        j_l_ext_sym=axis_contract(j_l_ext)
        j_l_axs=j_l_ext
        j_l_axs_sym=j_l_ext_sym
    else:
        j_l_axs=j_left
        j_l_axs_sym=j_l_sym

    fig,(ax1,ax2)=plt.subplots(2,1,figsize=(3.5,2.85))
    # Plot 1: regular coordinate
    ax1.set_title(out_title_gfs_o,fontsize=tlb_sz+1)
    xj_o,yi_o=np.meshgrid(i_b,j_l_axs)
    cs=ax1.pcolormesh(xj_o,yi_o,sgfs2d,cmap=cs_cmap,rasterized=True,vmin=cs_min,vmax=cs_max)
    ax1.tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz) 
#    ax1.set_xticks(i_b_xtk)
#    ax1.set_yticks(j_l_ytk)
    # extend(pointed end): 'neither'|'both'|'min'|'max'  
    divider=make_axes_locatable(ax1)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend='both')
    cbar.ax.tick_params(labelsize=tlb_sz)
    cbar.set_label(svar,fontsize=tlb_sz)

    # Plot 2: Contracted coordinate
    ax2.set_title(out_title_gfs_n,fontsize=tlb_sz+1)
    xj_n,yi_n=np.meshgrid(i_bx_sym,j_l_axs_sym)
    cs=ax2.pcolormesh(xj_n,yi_n,sgfs2d,cmap=cs_cmap,rasterized=True,vmin=cs_min,vmax=cs_max)
    ax2.tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz) 
    ax2.set_xticks(i_bx_xtk)
    ax2.set_xticklabels(i_bx_xlb)
    ax2.set_yticks(j_l_ytk)
    ax2.set_yticklabels(j_l_ylb)

    # extend(pointed end): 'neither'|'both'|'min'|'max'  
    divider=make_axes_locatable(ax2)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend='both')
    cbar.ax.tick_params(labelsize=tlb_sz)
    cbar.set_label(svar,fontsize=tlb_sz)


    fig.tight_layout()
    # Output figure
    out_file(out_gfs_fname)



# Contracted coordinate ============================================== CHJ =====
def axis_contract(ij_xy):
# ==================================================================== CHJ =====
    imid=(np.min(ij_xy)+np.max(ij_xy))/2.0 # find mid. pt.
    itmp=np.absolute(ij_xy-imid)
#    print(np.min(ij_xy),np.max(ij_xy),imid,itmp)

    if domain=='CONUS' and res=='C96' and n_blend==10:
        iexp=3
    elif domain=='CONUS' and res=='C96' and n_blend==6:
        iexp=5
    elif domain=='CONUS' and res=='C96' and n_blend==0:
        iexp=10
    elif domain=='BRZ' and res=='C768':
        iexp=16
    else:
        iexp=16

    itmp2=np.power(itmp,iexp)     # expand array from mid
    itmp=itmp2/np.max(itmp2)   # normalize
    ix_xy_cont=np.cumsum(itmp)    # new coord.

    return ix_xy_cont


 
# Ticks and lables for a contracted coordinate ======================= CHJ =====
def tick_label(n_ary,n_crd,o_crd):
# ==================================================================== CHJ =====
    # tick loc.
#    tck_lbl=[2,5,10,19,36,n_ary-39,n_ary-22,n_ary-13,n_ary-8,n_ary-5]
#    tck_lbl=[0,4,12,28,n_ary-29,n_ary-13,n_ary-5,n_ary-1]
    tck_lbl=[0,3,9,13,50,n_ary-51,n_ary-14,n_ary-10,n_ary-4,n_ary-1]
    # new ticks
    n_xtk=[n_crd[ir] for ir in tck_lbl]
    # new labels
    n_xlb=[o_crd[ir] for ir in tck_lbl]
    return n_xtk,n_xlb 



# Output file ======================================================== CHJ =====
def out_file(out_file):
# ==================================================================== CHJ =====
    # Output figure
    plt.savefig(out_fig_dir+out_file+'.png',dpi=300,bbox_inches='tight')
    plt.close('all')



# Main call ========================================================== CHJ =====
if __name__=='__main__':
    main()

