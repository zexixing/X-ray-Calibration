"""
read data from files or write data into files
"""

import numpy as np

def read_screened_column(table_path,
                         column_number,
                         screen_number,
                         screen_low,
                         screen_high,
                         ):
    """
    return a (screened) list from a table
    """
    data_screened = []
    column_output = np.loadtxt(table_path).T[column_number]
    column_screen = np.loadtxt(table_path).T[screen_number]
    dim = len(column_output) # read the number of rows of the table
                             # which can be used as number of cycles
    for i in range(0, dim):
        if (column_screen[i] > screen_low
            and column_screen[i] < screen_high):
            data_screened.append(column_output[i])
    return data_screened

def read_screened_index(table_path,
                        screen_number,
                        screen_low,
                        screen_high,
                        ):
    """
    return a list of indexes of elements screened according to some para
    """
    index_screened = []
    column_screen = np.loadtxt(table_path).T[screen_number]
    dim = len(column_screen) # read the number of rows of the table
                             # which can be used as number of cycles
    for i in range(0, dim):
        if (column_screen[i] > screen_low
            and column_screen[i] < screen_high):
            index_screened.append(i)
    return index_screened
