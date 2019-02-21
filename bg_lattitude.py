import _mypath
import os
import numpy as np
from lib.file_oper import read_screened_column, read_screened_index
from lib.list_data import bin_list, splited_list_in_box
from lib.plot.relation import relation_err

thisdir = os.path.dirname(__file__)
def bg_latitude():
    """
    plot relationship between background and latitude of the satellite
    """
    data_to_thisdir_path = '../data_output/crab/'
    table_list = os.path.join(thisdir, data_to_thisdir_path, 'list')
    f = open(table_list, 'r')
    L = f.readlines()
    f.close()
    la_bg = []
    latitude_low = -90
    latitude_high = 90
    latitude_step = 1
    bg_low = 8000
    bg_high = 10000
    for obs in L:
        print obs[:-1]
        table_path = os.path.join(thisdir, data_to_thisdir_path, obs[:-1])
        # get all latitudes filtered by energy
        latitude_screened = read_screened_column(table_path, 2, 1, bg_low, bg_high)
        la_bg = la_bg + latitude_screened
        # get unscreened observed table in every latitude box
        time_obs = read_screened_column(table_path, 0, 1, 0, 10000)
        latitude_obs = read_screened_column(table_path, 2, 1, 0, 10000)
        index_obs = read_screened_index(table_path, 1, 0, 10000)
        obs_table = np.column_stack((latitude_obs, time_obs, index_obs))
        boxes = np.arange(latitude_low,latitude_high + latitude_step,latitude_step)
        binned_obs = splited_list_in_box(obs_table, 0, boxes)
        # get obs time in every latitude
        bin_num = len(binned_obs)
        while obs == L[0]:
            binned_total_time = np.zeros((bin_num))
            break
        for bin_i in range(0, bin_num):
            binned_obs_i = binned_obs[bin_i]
            count_num = len(binned_obs_i)
            time_bin_i = 0
            if count_num > 1:
                for cou_j in range(1, count_num):
                    if binned_obs_i[cou_j][2]-binned_obs_i[cou_j-1][2] == 1:
                        delta_time_j = binned_obs_i[cou_j][1]-binned_obs_i[cou_j-1][1]
                        time_bin_i = time_bin_i + delta_time_j
            binned_total_time[bin_i] = binned_total_time[bin_i] + time_bin_i
    # bin the count list according to latitude
    bg_count, latitude = bin_list(la_bg, latitude_low, latitude_high, latitude_step)
    # get the count rate list
    bg_count = np.array(bg_count)
    bg_count_rate = []
    bg_err = []
    for bin_i in range(0, bin_num):
        if binned_total_time[bin_i] != 0:
            bg_count_rate.append(bg_count[bin_i]/binned_total_time[bin_i])
            bg_err.append(np.power(bg_count[bin_i],0.5)/binned_total_time[bin_i])
        else:
            bg_count_rate.append(0)
            bg_err.append(0)
    bg_rsd = np.std(bg_count_rate)/np.mean(bg_count_rate)
    # plot
    fig_save_path = os.path.join(thisdir, '../docs/')
    relation_err(latitude, bg_count_rate, err=bg_err, rsd=bg_rsd,
                 xlim_low=latitude_low, xlim_high=latitude_high,
                 ylim_low=-0.1, ylim_high=1.8,
                 xlabel='latitude (degree)', ylabel='background count rate',
                 title='background count rate to latitude (2017/01-2018/03)',
                 savepath=fig_save_path, savename='bg-cr_latitude.png',
                 )
    # savefile
    txt_save_name = os.path.join(fig_save_path, 'bg-cr_latitude.txt')
    np.savetxt(txt_save_name,
               np.column_stack((latitude,
                                bg_count,
                                binned_total_time,
                                bg_count_rate,
                                bg_err,
                                )),
               fmt = '%.4f %.4f %.4f %.4f %.4f',
               )


bg_latitude()
