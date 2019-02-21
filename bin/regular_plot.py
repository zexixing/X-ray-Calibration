import _mypath
import os
import numpy as np
import pylab as pl
from lib.file_oper import read_screened_column
from lib.list_data import bin_list, stat_list
from lib.plot.relation import relation_err

thisdir = os.path.dirname(__file__)
def count_time():
    """
    have a look for the event file
    """
    fig_save_path = '../docs/tycho/count_time/'
    data_to_thisdir_path = '../data_output/tycho/'
    list_name = 'list'
    table_list = os.path.join(thisdir, data_to_thisdir_path, list_name)
    f = open(table_list, 'r')
    L = f.readlines()
    f.close()
    for obs in L:
        print obs[:-1]
        table_path = os.path.join(thisdir, data_to_thisdir_path, obs[:-1])
        # get time and energy
        time_obs = read_screened_column(table_path, 0, 1, 0, 10000)
        energy_obs = read_screened_column(table_path, 1, 1, 0, 10000)
        energy_obs = np.array(energy_obs)/1000.0
        relation_err(time_obs, energy_obs,
                     ylim_low=0, ylim_high=10,
                     xlabel='time(s)', ylabel='energy(keV)',
                     title='energy-time for '+str(obs[:-4]),
                     savepath=fig_save_path, savename='energy-time_'+str(obs[:-4]),
                     )

def clip():
    """
    return flagged bg countrate txts and pngs
    """
    obs_low = 0.0
    obs_high = 10000.0
    bg_low = 8000.0
    bg_high = 10000.0
    time_step = 20.0
    countrate_ave_critical = 0.2
    average = []
    savetxt_to_thisdir_path = '../docs/tycho/clip_txt/'
    savefig_to_thisdir_path = '../docs/tycho/clip_fig/'
    data_to_thisdir_path = '../data_output/tycho/'
    list_name = 'list'
    savetxt_path = os.path.join(thisdir, savetxt_to_thisdir_path)
    savefig_path = os.path.join(thisdir, savefig_to_thisdir_path)
    table_list = os.path.join(thisdir, data_to_thisdir_path, list_name)
    f = open(table_list, 'r')
    L = f.readlines()
    f.close()
    for obs in L:
        print obs[:-1]
        table_path = os.path.join(thisdir, data_to_thisdir_path, obs[:-1])
        # get and bin bg time
        time_obs_photon = read_screened_column(table_path, 0, 1, obs_low, obs_high)
        time_bg_photon = read_screened_column(table_path, 0, 1, bg_low, bg_high)
        count_bg, time_bg = bin_list(time_bg_photon,
                                     time_obs_photon[0],
                                     time_obs_photon[-1],
                                     time_step,
                                     )
        countrate_bg = np.array(count_bg)/float(time_step)
        # define initial and final iteration by clipping once
        d = len(time_bg)
        F = 1
        flag = np.zeros(d)
        countrate_fin = []
        countrate_stat_ini = stat_list(countrate_bg, time_step)
        countrate_ave_ini = countrate_stat_ini['ave']
        countrate_err_ini = countrate_stat_ini['err']
        countrate_err_ave_ini = countrate_stat_ini['err_ave']
        average.append(countrate_ave_ini)
        for i in range(0,d):
            if countrate_bg[i] - countrate_ave_ini > 3*countrate_err_ini[i]:
                flag[i] = F
            else:
                countrate_fin.append(countrate_bg[i])
        d_fin = len(countrate_fin)
        if d_fin > 0:
            countrate_stat_fin = stat_list(countrate_fin, time_step)
            countrate_ave_fin = countrate_stat_fin['ave']
            countrate_err_fin = countrate_stat_fin['err']
            countrate_err_ave_fin = countrate_stat_fin['err_ave']
            average.append(countrate_ave_fin)
        elif d_fin == 0:
            savetxt_obs = str(obs[:-4])+'_clip_flag.txt'
            savetxt_obs_path = os.path.join(savetxt_path, savetxt_obs)
            f = open(savetxt_obs_path, 'w')
            for i in range(0, d):
                f.write(str(time_bg[i])+' '
                        +str(countrate_bg[i])+' '
                        +str(flag[i])+'\n')
                f.close()
        # loop to convergence, which means
        # the difference of initial ave and final ave is small, and
        # they are both less than critical average value
        while (abs(countrate_ave_ini - countrate_ave_fin) > countrate_err_ave_ini
               or countrate_ave_ini > countrate_ave_critical
               or countrate_ave_fin > countrate_ave_critical) \
               and countrate_ave_ini != countrate_ave_fin:
            countrate_ave_ini = countrate_ave_fin
            countrate_err_ave_ini = countrate_err_ave_fin
            countrate_fin = []
            F += 1
            for i in range(0, d):
                if flag[i] == 0:
                    if countrate_bg[i]-countrate_ave_ini>3*countrate_err_ini[i]:
                        flag[i] = F
                    else:
                        countrate_fin.append(countrate_bg[i])
            d_fin = len(countrate_fin)
            if d_fin > 0:
                countrate_stat_fin = stat_list(countrate_fin, time_step)
                countrate_ave_fin = countrate_stat_fin['ave']
                countrate_err_ave_fin = countrate_stat_fin['err_ave']
                average.append(countrate_ave_fin)
        # Or else, we can define a critical value.
        if countrate_ave_ini == countrate_ave_fin \
           and countrate_ave_ini > countrate_ave_critical:
            F += 1
            for i in range(0, d):
                if flag[i] == 0:
                    if countrate_bg[i]-countrate_ave_critical>3*countrate_err_ini[i]:
                        flag[i] = F
            average.append(countrate_ave_critical)
        # write in file
        savetxt_obs = str(obs[:-4])+'_clip_flag.txt'
        savetxt_obs_path = os.path.join(savetxt_path, savetxt_obs)
        f = open(savetxt_obs_path, 'w')
        for i in range(0, d):
            f.write(str(time_bg[i])+' '
                    +str(countrate_bg[i])+' '
                    +str(flag[i])+'\n')
        f.close()
        # plot
        t_0=[]
        cr_0=[]
        t_1=[]
        cr_1=[]
        t_2=[]
        cr_2=[]
        t_3=[]
        cr_3=[]
        for j in range(0,d):
            if flag[j]==1:
                t_1.append(time_bg[j])
                cr_1.append(countrate_bg[j])           
            if flag[j]==2:
                t_2.append(time_bg[j])
                cr_2.append(countrate_bg[j])
            if flag[j]>2:
                t_3.append(time_bg[j])
                cr_3.append(countrate_bg[j])
            if flag[j]==0:
                t_0.append(time_bg[j])
                cr_0.append(countrate_bg[j])
        pl.plot(t_0,cr_0,'ko',MarkerSize=2,label='unmasked')
        if len(t_1)>0:
            pl.plot(t_1,cr_1,'ro',MarkerSize=2,label='masked in the 1st round')
            pl.plot([time_bg[0],time_bg[-1]],[average[1],average[1]], color ='red', linewidth=1.0, linestyle="--")
        if len(t_2)>0:
            pl.plot(t_2,cr_2,'bo',MarkerSize=2,label='masked in the 2st round')
            pl.plot([time_bg[0],time_bg[-1]],[average[2],average[2]], color ='blue', linewidth=1.0, linestyle="--")
        if len(t_3)>0:
            pl.plot(t_3,cr_3,'yo',MarkerSize=2,label='masked in other rounds')
            pl.plot([time_bg[0],time_bg[-1]],[average[3],average[3]], color ='y', linewidth=1.0, linestyle="--")
        if countrate_ave_critical in average:
            pl.plot([time_bg[0],time_bg[-1]],[countrate_ave_critical,countrate_ave_critical], color ='k', linewidth=1.0, linestyle="--")
        pl.legend(loc='upper left',prop={'size':6})
        pl.xlabel('time (s)')
        pl.ylabel('counts rate per 20s (counts/20s)')
        pl.title('counts rate - time of '+str(obs[:-4]))
        savefig_obs = str(obs[:-4])+'_clip_flag.png'
        savefig_obs_path = os.path.join(savefig_path, savefig_obs)
        pl.savefig(savefig_obs_path)
        pl.close('all')

def lightcurve():
    """
    return energies and the lightcurve in good time
    """
    energy_low = 500
    energy_high = 5000
    CO=[]
    TI=[]
    TM=[]
    CR=[]
    ER=[]
    good = 0
    bad = 0
    fig_to_thisdir_path = '../docs/casa/'
    energy_to_thisdir_path = '../docs/casa/energy'
    data_to_thisdir_path = '../data_output/casa/'
    clipped_data_to_thisdir_path = '../docs/casa/clip_txt/'
    list_name = 'list'
    fig_save_path = os.path.join(thisdir, fig_to_thisdir_path)
    energy_save_path = os.path.join(thisdir, energy_to_thisdir_path)
    table_list = os.path.join(thisdir, data_to_thisdir_path, list_name)
    f = open(table_list, 'r')
    L = f.readlines()
    f.close()
    f = open(energy_save_path,'w')
    f.close()
    for i in L:
        obs = i[:-1]
        clip = i[:-4]+'_clip_flag.txt'
        clip = os.path.join(thisdir, clipped_data_to_thisdir_path, clip)
        obs = os.path.join(thisdir, data_to_thisdir_path, obs)
        t_clip = np.loadtxt(clip).T[0]
        flag = np.loadtxt(clip).T[2]
        t_obs = np.loadtxt(obs).T[0]
        e_obs = np.loadtxt(obs).T[1]
        step=t_clip[1]-t_clip[0]
        d_clip = len(t_clip)
        d_obs = len(t_obs)
        
        t_good_mark = []
        t_bad_mark = []
        t_good = []
        t_bad = []
        e_good = []
        e_bad = []

        for j in range(0,d_clip):
            if flag[j]==0:
                t_good_mark.append(t_clip[j])
                good += 1
            else:
                t_bad_mark.append(t_clip[j])
                bad += 1
        for j in range(0,d_obs):
            if len(t_bad_mark)!=0:
                cri = min(abs(np.array(t_bad_mark)-t_obs[j]))
                if cri<(step/2.0):
                    t_bad.append(t_obs[j])
                    e_bad.append(e_obs[j])
                elif (t_obs[j]>t_clip[0]-step/2.0
                      or t_obs[j]==t_clip[0]-step/2.0) \
                     and \
                     (t_obs[j]<t_clip[-1]+step/2.0
                      or t_obs[j]==t_clip[-1]+step/2.0
                      ):
                    t_good.append(t_obs[j])
                    e_good.append(e_obs[j])
            elif (t_obs[j]>t_clip[0]-step/2.0
                  or t_obs[j]==t_clip[0]-step/2.0) \
                 and \
                 (t_obs[j]<t_clip[-1]+step/2.0
                  or t_obs[j]==t_clip[-1]+step/2.0
                  ):
                t_good.append(t_obs[j])
                e_good.append(e_obs[j])
        d_good = len(t_good)
        count = 0
        for j in range(0,d_good):
            if e_good[j]>energy_low and e_good[j]<energy_high:
                count+=1
        tinterval = len(t_good_mark)*20.0
        CO.append(count)
        TI.append(tinterval)
        TM.append((t_obs[0]+t_obs[-1])/2.0)
        CR.append(count/tinterval)
        ER.append(np.power(count,0.5)/tinterval)       
        f = open(energy_save_path,'a')
        for j in range(0,d_good):
            f.write(str(e_good[j])+'\n')
        f.close()
    RSD = np.std(CR)/np.mean(CR)
    relation_err(TM, CR, err=ER, rsd=RSD,
                 xlabel='time(s)', ylabel='clipped count rate per observation (counts/obs)',
                 title='count rate-time',
                 savepath=fig_save_path, savename='lightcurve',
                 text_xloc=0.7, text_yloc=0.1,
                 )
    print '%.4f'%(float(good)/(float(good+bad)))

def sed_total():
    """
    return SED in good time
    """
    clipped_data_to_thisdir_path = '../docs/casa/clip_txt/'
    clipped_data_path = os.path.join(thisdir, clipped_data_to_thisdir_path, 'list')
    save_to_thisdir_path = '../docs/casa/'
    energy_path = os.path.join(thisdir, save_to_thisdir_path, 'energy')
    savefig_path = os.path.join(thisdir, save_to_thisdir_path)
    # calculate good time
    f = open(clipped_data_path, 'r')
    cliplist = f.readlines()
    f.close()
    total_good_time = 0
    for i in cliplist:
        clip = os.path.join(thisdir, clipped_data_to_thisdir_path, i[:-1])
        t = np.loadtxt(clip).T[0]
        flag = np.loadtxt(clip).T[2]
        d = len(flag)
        step = t[1] - t[0]
        for j in range(0, d):
            if flag[j]==0:
                total_good_time += step
    # calculate sed using histogram
    bin_low = 0.0
    bin_high = 10000.0
    bin_step = 20.0
    delta_e = bin_step/1000.0
    energy_in_good_time = np.loadtxt(energy_path)
    count, energy_bin = bin_list(energy_in_good_time, bin_low, bin_high, bin_step)
    countrate = count/(total_good_time*delta_e)
    error = np.sqrt(np.array(count))/(total_good_time*delta_e)
    energy_bin = energy_bin/1000.0
    relation_err(energy_bin, countrate, err=error,
                 xlim_low=0.4, xlim_high=10,
                 ylim_high=1.5,
                 xscale='log', yscale='log',
                 xlabel='kev', ylabel='couts s-1 kev-1',
                 title='Spectrum of casa',
                 savepath=savefig_path, savename='sed_casa',
                 )
    # output in text file
    savetxt_path = os.path.join(thisdir, save_to_thisdir_path, 'sed_casa_txt')
    f = open(savetxt_path, 'w')
    for i in range(0, len(error)):
        if error[i] != 0.0:
            f.write(str(energy_bin[i])+' '+str(countrate[i])+' '+str(error[i])+'\n')
    f.close()
    
sed_total()
