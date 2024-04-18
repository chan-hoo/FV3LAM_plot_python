###################################################################### CHJ #####
## Name		: plot_fv3lam_grb2his.py
## Language	: Python 3.7
## Usage	: Plot time-history of max. precipitation/reflectivity from 'BGDAWP(grib2)'
## Input files  : BGDAWPXX.grib2
## NOAA/NWS/NCEP/EMC
## History ===============================
## V000: 2021/03/08: Chan-Hoo Jeon : Preliminary version
## V001: 2021/03/10: Chan-Hoo Jeon : Change annotation format
## V002: 2021/03/10: Chan-Hoo Jeon : Change the fields
###################################################################### CHJ #####

import os, sys
import pygrib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np


# HPC machine ('hera','wcoss_dell','wcoss_cray','orion')
machine='hera'

print(' You are on', machine)

#### Machine-specific input data ==================================== CHJ =====
# out_fig_dir: directory where the output files are created

if machine=='hera':
    out_fig_dir="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/tools/fv3sar_pre_plot/Fig/"
elif machine=='wcoss_dell':
    out_fig_dir="/gpfs/dell2/emc/modeling/noscrub/Chan-Hoo.Jeon/tools/fv3_lam_plot/"
elif machine=='wcoss_cray':
    out_fig_dir="/gpfs/hps3/emc/meso/noscrub/Chan-Hoo.Jeon/tools/fv3_lam_plot/"
elif machine=='orion':
    out_fig_dir="/work/noaa/fv3-cam/chjeon/tools/Fig/"
else:
    sys.exit('ERROR: the machine name does not exsit !!!')

plt.switch_backend('agg')

# Case-dependent input =============================================== CHJ =====
# Path to the directory where the input NetCDF files are located.
dnm_data="/scratch2/NCEPDEV/fv3-cam/Chan-hoo.Jeon/srw_dev/expt_dirs/test_GST/2019061500/postprd/"

# Input file hour range (start, end, interval)
fnm_hr_s=0
fnm_hr_e=28
fnm_hr_i=1

# basic forms of title and file name
out_title='FV3LAM :: Time-history of Max. Temperature :: '+machine.upper()
out_fname='fv3lam_out_grb2his_'+machine



# Main part (will be called at the end) ==================== CHJ =====
def main():
# ========================================================== CHJ =====

    his_tpr_max=[]
    his_rfl_max=[]
    his_xlabel=[]
    for fhr in range(fnm_hr_s,fnm_hr_e+1,fnm_hr_i):
        fnm_hr = str(fhr).zfill(2)
        his_xlabel.append(fnm_hr)
        # Input file
        fnm_in='rrfs.t00z.bgdawpf0'+fnm_hr+'.tm00.grib2'
        print(' === '+fnm_hr+' === in: '+fnm_in+' ====================')

        # open the data file
        fname=os.path.join(dnm_data,fnm_in)
        try: grbs=pygrib.open(fname)
        except: raise Exception('Could NOT find the file',fname)

#        for grb in grbs:
#            print(grb.name)
#            print(grb.typeOfLevel)
#            print(grb.level)
#            print(grb.Nx)
#            print(grb.Ny)
#        exit()

#        print(' *** 2 meter Temperature *** ')
        grb_tpr=grbs.select(name="2 metre temperature")[0].values
        grb_tpr[grb_tpr<=0]=np.nan
        grb_tpr=(grb_tpr-273.15)*1.8+32.0
        tpr_max=np.nanmax(grb_tpr)
        print(' var1 max = ',tpr_max)
        his_tpr_max.append(tpr_max)

#        print(' *** Temperature at 100m *** ')
        grb_rfl=grbs.select(name="Temperature",typeOfLevel="heightAboveGround")[3].values
        grb_rfl[grb_rfl<=0]=np.nan
        grb_rfl=(grb_rfl-273.15)*1.8+32.0
        rfl_max=np.nanmax(grb_rfl)
        print(' var2 max = ',rfl_max)
        his_rfl_max.append(rfl_max)

    print(his_tpr_max)
    print(his_rfl_max)

    x1=np.arange(len(his_xlabel))

    ax_tpr_max = np.nanmax(his_tpr_max)*1.1
    ax_tpr_min = 75.0
    ax_rfl_max = np.nanmax(his_rfl_max)*1.1
    ax_rfl_min = 75.0

    txt_fnt=6
    bar_wdth=0.07
    ann_fnt=3

    fig,(ax1,ax2)=plt.subplots(2,1,figsize=(6,3))
    fig.suptitle(out_title,fontsize=txt_fnt+1,y=1.01)
    rects1=ax1.bar(x1,his_tpr_max,bar_wdth,alpha=0.7)
    ax1.set_ylabel('Max. 2m temperature (F)',fontsize=txt_fnt-1)
    ax1.tick_params(axis="y",labelsize=txt_fnt-2)
    ax1.set_xticks(x1)
    ax1.set_xticklabels(his_xlabel,fontsize=txt_fnt-2)
    ax1.set_ylim([ax_tpr_min,ax_tpr_max])

    rects2=ax2.bar(x1,his_rfl_max,bar_wdth,alpha=0.7)
    ax2.set_ylabel('Max. temperature at 100m (F)',fontsize=txt_fnt-1)
    ax2.set_xlabel('Output file hour (fhr)',fontsize=txt_fnt-1)
    ax2.tick_params(axis="y",labelsize=txt_fnt-2)
    ax2.set_xticks(x1)
    ax2.set_xticklabels(his_xlabel,fontsize=txt_fnt-2)
    ax2.set_ylim([ax_rfl_min,ax_rfl_max])

    def autolabel(ax,rects):
        for rect in rects:
            height=rect.get_height()
            ax.annotate('{:.2f}'.format(height),
            xy=(rect.get_x()+rect.get_width()/2,height),
            xytext=(0,0.5), textcoords="offset points",
            ha='center', va='bottom',color='red',fontsize=ann_fnt)

    autolabel(ax1,rects1)
    autolabel(ax2,rects2)
    fig.tight_layout()


# Output figure
    ndpi=300
    out_file(out_fname,ndpi) 


      
# Output file ============================================= CHJ =====
def out_file(out_file,ndpi):
# ========================================================= CHJ =====
    # Output figure
    plt.savefig(out_fig_dir+out_file+'.png',dpi=ndpi,bbox_inches='tight')
    plt.close('all')



# Main call ================================================ CHJ =====
if __name__=='__main__':
    main()

