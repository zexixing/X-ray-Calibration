"""
operate data in structure of list
"""
import numpy as np

def bin_list(data_list,
             bin_low,
             bin_high,
             bin_step
             ):
    """
    return a binned list and related bin space list

    this function use the principle of histogram and
    can be used to plot SED, count-latitude, etc..
    """
    bin_ends = np.arange(bin_low,bin_high+bin_step,bin_step)
    bin_height, bin_ends = np.histogram(data_list, bin_ends)
    bin_list = (bin_ends[:-1]+bin_ends[1:])/2.0
    return bin_height, bin_list
    
def splited_list_in_box(data_list, box_index, boxes):
    """
    return a 3-d array in box order
    e.g., [[box1:[elem1],[elem2]..], [box2:[e1],[e2]..],..]

    data_list should be a 2-d array,
    e.g., [[e1: box_1, a_1, b_1], [e2: box_2, a_2, b_2],..]
    where value of box should be data_list[:][box_index]

    if data_list which we have is not a 2-d array but n 1-d lists
    we can use the following codes to merge them
    data_list = np.column_stack((list1,list2,list3,...))
    """
    elems_in_every_box = []
    box_num = len(boxes) -1 # eg. the fmt of input box: [a, b, c], 3 elements
                            # thus the boxes are [a, b] and [b, c], 2 boxes
    for i in range(0, box_num):
        box_i_begin = boxes[i]
        box_i_end = boxes[i+1]
        if_in_box_low = data_list[:, box_index] > box_i_begin
        if_in_box_high = data_list[:, box_index] < box_i_end
        if_in_box = np.logical_and(if_in_box_low, if_in_box_high)
        elems_in_box_i = data_list[if_in_box]
        elems_in_every_box = elems_in_every_box + [elems_in_box_i]
    return elems_in_every_box

def stat_list(data_list, step):
    """
    return useful statistical parameters in astro,
    such as average, error, mean error, ...
    updating
    """
    d = len(data_list)
    stat_para = {}
    stat_para['ave'] = np.mean(data_list)
    stat_para['err'] = np.sqrt(data_list)/np.sqrt(step)
    stat_para['err_ave'] = np.sqrt(np.sum(np.power(data_list,2)))/float(d)
    return stat_para
