# -*- coding: utf-8 -*-
"""
simple relation plots with 2-dimension information
"""
import matplotlib.pylab as pl

def relation_err(x, y,
                 xlabel, ylabel, title,
                 savepath, savename,
                 style = 'ko', size = 1, trans = 0.3,
                 xscale = 'linear', yscale = 'linear',
                 text_xloc = 0.7, text_yloc = 0.9,
                 **kwargs
                 ):
    """
    plot the most common figures and save
    
    errors can be plotted if needed
    xlim and ylim can be setted if needed
    """
    fig = pl.figure()
    ax = fig.add_subplot(1,1,1)
    pl.plot(x, y, style, MarkerSize=size)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)
    if len(kwargs) != 0:
        if 'xlim_low' in kwargs:
            pl.xlim(left=kwargs['xlim_low'])
        if 'xlim_high' in kwargs:
            pl.xlim(right=kwargs['xlim_high'])
        if 'ylim_low' in kwargs:
            pl.ylim(bottom=kwargs['ylim_low'])
        if 'ylim_high' in kwargs:
            pl.ylim(top=kwargs['ylim_high'])
        if 'err' in kwargs:
            pl.errorbar(x, y, yerr=[kwargs['err'],kwargs['err']],
                        fmt=style, MarkerSize=size, alpha = trans)
        if 'rsd' in kwargs:
            pl.text(text_xloc, text_yloc,
                    r'$\sigma_{relative}$='+str('%.4f'%kwargs['rsd']),
                    transform=ax.transAxes)
    pl.title(title)
#    pl.savefig(savepath+savename)
    pl.show()

import numpy as np
path = '/Users/zexixing/research_xpnav-01_spag/pulsar/B0531+21/positive-0.2-clip/result/EA/ea.txt'
ener = np.loadtxt(path).T[0]
area = np.loadtxt(path).T[1]
relation_err(ener, area,
             'eV', 'cm2', 'eff area',
             ' ',' ', style = 'k-', 
             xscale='log', yscale='log',
             )
