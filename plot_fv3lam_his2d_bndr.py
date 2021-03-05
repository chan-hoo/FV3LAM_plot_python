###################################################################### CHJ #####
## Name		: plot_fv3lam_his2d_bndr.py
## Language	: Python 3.7
## Usage	: Plot boundary of domain from fv3_history2d.nc (output)
## Input files  : fv3_history2d.nc
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2020/12/24: Chan-Hoo Jeon : Preliminary version
## V001: 2021/03/05: Chan-Hoo Jeon : Simplify the script
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

# HPC machine ('hera','orion')
machine='hera'

print(' You are on', machine)

plt.switch_backend('agg')

# Case-dependent input =============================================== CHJ =====
# Path to the directory where the input NetCDF file is located.
dnm_in="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/test/"

# grid file name
fnm_input='fv3_history2d.nc'

# Variables
vars_his=["tmpsfc"]

# Time level of plotting ( 0 < prt_tlvl < maximum-1)
prt_tlvl=1

# Length of boundary (unit: number of grid points)
len_bndr=20

# tick interval for short sides
dtick_s=2
# tick interval for long sides
dtick_l=200

# Path to directory
if machine=='hera':
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig/"
    mfdt_kwargs={'parallel':False}
elif machine=='orion':
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
    mfdt_kwargs={'parallel':False,'combine':'by_coords'}
else:
    sys.exit('ERROR: path to output directory is not set !!!')

# basic forms of title and file name
out_title_base='FV3LAM::HIS2D::'
out_fname_base='fv3lam_out_his2d_'


# Main part (will be called at the end) ======================= CHJ =====
def main():
# ============================================================= CHJ =====
    global his2d

    print(' ===== OUTPUT: history2d =======================================')
    # open the data file
    fname=os.path.join(dnm_in,fnm_input)
    try: his2d=xr.open_mfdataset(fname,**mfdt_kwargs)
    except: raise Exception('Could NOT find the file',fname)
    print(his2d)

    # Variables
    for svar in vars_his:
        his_plot(svar)



# ===== plot ================================================== CHJ =====
def his_plot(svar):
# ============================================================= CHJ =====

    global out_title,out_fname,cs_min,cs_max,sfld2d
    global cs_cmap,lb_ext,tick_ln,tick_wd,tlb_sz,n_rnd
    global nts,nys,nxs

    print(' ===== '+svar+' ===== history2d ===============================')
    # Extract data array
    sfld=np.ma.masked_invalid(his2d[svar].data)

    (nts,nys,nxs)=sfld.shape
    print(' time+2D: nts=',nts,' nys=',nys,' nxs=',nxs)

    if prt_tlvl>=nts:
        sys.exit('ERROR: prt_tlvl >= max. time level !!!')


    sfld2d=sfld[prt_tlvl,:,:]

    out_title=out_title_base+svar+'::time level='+str(prt_tlvl)
    out_fname=out_fname_base+svar

    cs_cmap='jet'
    lb_ext='neither'
    tick_ln=1.5
    tick_wd=0.45
    tlb_sz=5
    n_rnd=5
    cmap_range='round'
 

    print(' Plotting field=',svar)

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


    print(' Plotting the entire domain ... ')
    plot_domain(svar)

    print(' Plotting boundaries ... ')
    plot_bndry(svar)



# Plot: Boundary ========================================== CHJ =====
def plot_bndry(svar):
# ========================================================= CHJ =====

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


    len_bndr_m1=len_bndr-1

    # tick: bottom
    tick_b_x=np.arange(0,nxs,dtick_l)
    tick_lbl_b_x=np.arange(1,nxs+1,dtick_l)
    tick_b_y=np.arange(0,len_bndr,dtick_s)
    tick_lbl_b_y=np.arange(1,len_bndr+1,dtick_s)

    # tick: top
    tick_t_x=np.arange(0,nxs,dtick_l)
    tick_lbl_t_x=np.arange(1,nxs+1,dtick_l)
    tick_t_y=np.arange(0,len_bndr,dtick_s)
    tick_lbl_t_y=np.arange(nys-len_bndr_m1,nys+1,dtick_s)
   
    # tick: left
    tick_l_x=np.arange(0,len_bndr,dtick_s)
    tick_lbl_l_x=np.arange(1,len_bndr+1,dtick_s)
    tick_l_y=np.arange(0,nys,dtick_l)
    tick_lbl_l_y=np.arange(1,nys+1,dtick_l)
 
    # tick: right
    tick_r_x=np.arange(0,len_bndr,dtick_s)
    tick_lbl_r_x=np.arange(nxs-len_bndr_m1,nxs+1,dtick_s)
    tick_r_y=np.arange(0,nys,dtick_l)
    tick_lbl_r_y=np.arange(1,nys+1,dtick_l)
 

    # bottom (North)
    ax_b.pcolormesh(sfld2d[:len_bndr-1,:],cmap=cs_cmap,rasterized=True,vmin=cs_min,vmax=cs_max)
    ax_b.tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz)
    ax_b.set_xticks(tick_b_x)
    ax_b.set_xticklabels(tick_lbl_b_x)
    ax_b.set_yticks(tick_b_y)
    ax_b.set_yticklabels(tick_lbl_b_y)

    # top (South)
    ax_t.pcolormesh(sfld2d[-len_bndr_m1:,:],cmap=cs_cmap,rasterized=True,vmin=cs_min,vmax=cs_max)
    ax_t.tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz,
                     bottom=False,labelbottom=False,top=True,labeltop=False)
    ax_t.set_xticks(tick_t_x)
    ax_t.set_xticklabels(tick_lbl_t_x)
    ax_t.set_yticks(tick_t_y)
    ax_t.set_yticklabels(tick_lbl_t_y)

    # left (East)
    ax_l.pcolormesh(sfld2d[len_bndr:-len_bndr,:len_bndr-1],cmap=cs_cmap,rasterized=True,vmin=cs_min,vmax=cs_max)
    ax_l.tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz,
                     labelbottom=False)

    ax_l.set_xticks(tick_l_x)
    ax_l.set_xticklabels(tick_lbl_l_x)
    ax_l.set_yticks(tick_l_y)
    ax_l.set_yticklabels(tick_lbl_l_y)

    # right (West)
    ax_r.pcolormesh(sfld2d[len_bndr:-len_bndr,-len_bndr_m1:],cmap=cs_cmap,rasterized=True,vmin=cs_min,vmax=cs_max)
    ax_r.tick_params(direction='out',length=tick_ln,width=tick_wd,
                     left=False,labelleft=False,labelbottom=False,right=True)
    ax_r.set_xticks(tick_r_x)
    ax_r.set_xticklabels(tick_lbl_r_x)
    ax_r.set_yticks(tick_r_y)
    ax_r.set_yticklabels(tick_lbl_r_y)


    # center
    ax_c.axis([0,10,0,10])
    ax_c.axis('off')

    ax_c.text(5,0.05,'North (bottom)',fontsize=tlb_sz,horizontalalignment='center',
              verticalalignment='center')
    ax_c.text(5,9.95,'South (top)',fontsize=tlb_sz,horizontalalignment='center',
              verticalalignment='center')
    ax_c.text(0.05,5,'East (left)',fontsize=tlb_sz,horizontalalignment='center',
              verticalalignment='center',rotation='vertical')
    ax_c.text(9.95,5,'West (right)',fontsize=tlb_sz,horizontalalignment='center',
              verticalalignment='center',rotation='vertical')

    ax_c.text(5,7,out_title_base,fontsize=tlb_sz+1.5,color='black',
              horizontalalignment='center',verticalalignment='center')
    ax_c.text(5,6,'Variable='+svar,fontsize=tlb_sz+1.5,color='black',
              horizontalalignment='center',verticalalignment='center')
    ax_c.text(5,5,'Time level='+str(prt_tlvl),fontsize=tlb_sz+1.5,color='black',
              horizontalalignment='center',verticalalignment='center')
    ax_c.text(5,4,'Boundary cells='+str(len_bndr),fontsize=tlb_sz+1.5,color='black',
              horizontalalignment='center',verticalalignment='center')


    out_fname_bndr=out_fname+'_bndr'

    # Output figure
    ndpi=300
    out_file(out_fname_bndr,ndpi)




# Plot: Entire domain ===================================== CHJ =====
def plot_domain(svar):
# ========================================================= CHJ =====

    tick_x=np.arange(0,nxs,dtick_l)
    tick_lbl_x=np.arange(1,nxs+1,dtick_l)
    tick_y=np.arange(0,nys,dtick_l)
    tick_lbl_y=np.arange(1,nys+1,dtick_l)

    fig,ax=plt.subplots(1,1,figsize=(3.5,2.85))
    ax.set_title(out_title, fontsize=tlb_sz+1)

    cs=ax.pcolormesh(sfld2d,cmap=cs_cmap,rasterized=True,vmin=cs_min,vmax=cs_max)
    ax.tick_params(direction='out',length=tick_ln,width=tick_wd,labelsize=tlb_sz)
    ax.set_xticks(tick_x)
    ax.set_xticklabels(tick_lbl_x)
    ax.set_yticks(tick_y)
    ax.set_yticklabels(tick_lbl_y)

    divider=make_axes_locatable(ax)
    ax_cb=divider.new_horizontal(size="3%",pad=0.1,axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar=plt.colorbar(cs,cax=ax_cb,extend=lb_ext)
    cbar.ax.tick_params(labelsize=tlb_sz)
    cbar.set_label(svar,fontsize=tlb_sz)

    out_fname_domain=out_fname+'_domain'

    # Output figure
    ndpi=300
    out_file(out_fname_domain,ndpi)

  

# Output file ============================================= CHJ =====
def out_file(out_fname,ndpi):
# ========================================================= CHJ =====
    # Output figure
    plt.savefig(out_fig_dir+out_fname+'.png',dpi=ndpi,bbox_inches='tight')
    plt.close('all')



# Main call ================================================ CHJ =====
if __name__=='__main__':
    main()

