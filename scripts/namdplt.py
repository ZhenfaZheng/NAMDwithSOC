#!/usr/bin/env python

import math
import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
from glob import glob


def main():

    Eref = 9.3494

    inp = read_inp('inp')
    soctype = int(inp['SOCTYPE'])
    potim = float(inp['POTIM'])

    filcoup = 'NATXT'
    figname = 'COUPLE_NA.png'
    coup = read_couple(filcoup, ctype=1)
    coup_av = np.average(np.abs(coup), axis=0)
    hbar = 0.6582119281559802 # ev.fs
    coup_av = coup_av * hbar / ( 2 * potim )
    n = coup_av.shape[0]
    for ii in range(n):
        coup_av[ii,ii] = 0.0
    plot_couple(coup_av, figname)

    if (soctype==2):
        filcoup = 'SOTXT'
        figname = 'COUPLE_SO.png'
        coup = read_couple(filcoup, ctype=1)
        coup_av = np.average(np.abs(coup), axis=0)
        n = coup_av.shape[0]
        for ii in range(n):
            coup_av[ii,ii] = 0.0
        plot_couple(coup_av, figname)

    figname = 'NAMD.png'
    filshps = glob('SHPROP.*')
    shp = readshp(filshps)
    ntsteps = shp.shape[0]
    ksen = np.loadtxt('EIGTXT')[0:ntsteps, :]
    plot_tdprop(shp, Eref, lplot=2, ksen=ksen, figname=figname)


def read_inp(infile='inp'):
    '''
    Read input parameters.

    Parameters:
    infile: string, coupling file.

    Returns: dictionary, the keys and values are all strings.
    '''

    text = [line for line in open(infile) if line.strip()]

    inp = {}
    for line in text:
        if (line[0]=='&' or line[0]=='/' or line[0]=='!'):
            continue
        temp = line.split('=')
        key = temp[0].strip()
        value = temp[1].strip().strip('\'').strip('\"')
        inp[key] = value

    return inp


def read_couple(filcoup='NATXT', ctype=0):
    '''
    This function loads data from NATXT file.

    Parameters:
    filcoup  : string, coupling file.
    ctype    : integer, different forms of NATXT files.
               0: origin type, NAC data are real number;
               1: NAC are complex, and restore in two real numbers;

    Returns: ndarray, coupling data in forms of coup[nsw-1, nb, nb]
    '''

    if ctype==0:
        coup = np.loadtxt(filcoup)
        try:
            nt = int(coup.shape[0])
            nb = int( np.sqrt(coup.shape[1]) )
        except IndexError:
            nt = 1
            nb = int( np.sqrt(coup.shape[0]) )

        coup.resize(nt,nb,nb)

    elif ctype==1:
        data = np.loadtxt(filcoup)
        try:
            nt = int(data.shape[0])
            nb = int( np.sqrt(data.shape[1]/2) )
            coup = data[:,0::2] + data[:,1::2]*(1.0j)
        except IndexError:
            nt = 1
            nb = int( np.sqrt(data.shape[0]/2) )
            coup = data[0::2] + data[1::2]*(1.0j)

        coup.resize(nt,nb,nb)

    return coup


def read_inicon(filinicon='INICON'):
    '''
    This function loads data from INICON file.
    '''
    inicon = np.loadtxt(filinicon)
    return inicon


def readshp(filshps):
    '''
    This function loads data from SHPROP.xxx files.

    Parameters:
    filshps: a list of strings, file names if SHPROP.xxx files. such as
             ['SHPROP.1', 'SHPROP.5']

    Returns: ndarray, average data of SHPROP.xxx files, in forms of
             shp[ntsteps, nb+2].
    '''
    shps = np.array( [ np.loadtxt(filshp) for filshp in filshps ] )
    shp = np.average(shps, axis=0)
    return shp


def plot_couple(coup, figname='COUPLE.png'):
    '''
    This function plots average couplings.

    Parameters:
    coup: ndarray, average coupling data in forms of coup[nb, nb]
    figname: string, output figure file name.
    '''

    fig = plt.figure()
    figsize_x = 4.8
    figsize_y = 3.6 # in inches
    fig.set_size_inches(figsize_x, figsize_y)

    cmap = 'bwr'
    n = coup.shape[0]
    coup *= 1000.0 # change unit to meV
    Bmin = 0.5; Bmax = n + 0.5
    cmin = 0.0; cmax = np.max(coup)
    norm = mpl.colors.Normalize(cmin,cmax)
    plt.imshow(coup, cmap=cmap, origin='lower', norm=norm,
        extent=(Bmin,Bmax,Bmin,Bmax), interpolation='none')

    cbar = plt.colorbar()
    # cbar.ax.set_title('   meV')
    cbar.set_label('Coupling (meV)')
    plt.tight_layout()
    plt.savefig(figname, dpi=400)


def plot_tdprop(shp, Eref=0.0, lplot=1, ksen=None, figname='tdshp.png'):
    '''
    This function loads data from SHPROP.xxx files,
    and plot average evolution of energy & surface hopping proportion of
    electronic states.

    Parameters:
    shp    : ndarray, average data of SHPROP.xxx files, in forms of
             shp[ntsteps, nb+2].
    Eref   : float, energy reference. Make sure shp & ksen have same Eref!!!
    lplot  : integer, fig type to plot. 1: plot average proportions; 2: plot
             average energy evolution with proportions.
    ksen   : ndarray, KS energies in forms of ksen[ntsteps, nbands]. Here we
             suppose ksen do not change by time.
    figname: string, file name of output figure.
    '''

    figsize_x = 4.8
    figsize_y = 3.2 # in inches
    namdtime = shp[-1,0]

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    if lplot==1:
        ylabel = 'SHPROP'
        ax.plot(shp[:,0], shp[:,2:])
        ax.set_ylim(0.0,1.0)
        ax.set_ylabel(ylabel)
    else:
        cmap = 'hot_r'
        dotsize = 20
        ylabel = 'Energy (eV)'
        ntsteps = shp.shape[0]
        nbands = shp.shape[1] -2
        cmin = 0.0; cmax = np.max(shp[:,2:])
        cmax = math.ceil(cmax*10)/10
        norm = mpl.colors.Normalize(cmin,cmax)

        if (ksen.shape[1]!=nbands):
            print('\nNumber of ksen states doesn\'t match with SHPROP data!\n')
        T = np.tile(shp[:,0], nbands).reshape(nbands,ntsteps).T
        sc = ax.scatter(T, ksen-Eref, s=dotsize, c=shp[:,2:], lw=0,
                        norm=norm, cmap=cmap)
        ax.plot(shp[:,1]-Eref, 'r', lw=1, label='Average Energy')
        plt.colorbar(sc)

        ax.set_ylabel(ylabel)

    ax.set_xlim(0,namdtime)
    ax.set_xlabel('Time (fs)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)


if __name__=='__main__':
    main()
