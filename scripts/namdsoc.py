#/usr/bin/env python

import math, os, re
import numpy as np
import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
from glob import glob


def main():

    #                       Basic parameters setting                      #
    #######################################################################

    which_plt = [1, 2, 3, 4]
    '''
    Select which figures to plot.
    1: COUPLE_NA.png (COUPLE_SO.png); 2: TDEN.png; 3: TDPOP.png.
    4: TDKSEN.png
    '''

    ldata = True
    pathD = 'Data/'
    '''
    if ldata=True, will extract band information from run directory,
    then create folder pathD, and save the data in pathD directory.
    if ldata=False, will read data from pathD.
    '''

    #######################################################################

    inp = read_inp('inp')

    if ldata and (not which_plt==[1]):
        loadData(inp, pathD)

    if (1 in which_plt):

        coup = read_couple('NATXT', ctype=1)
        coup_av = np.average(np.abs(coup), axis=0)
        hbar = 0.6582119281559802 # ev.fs
        potim = float(inp['POTIM'])
        coup_av = coup_av * hbar / ( 2 * potim )
        n = coup_av.shape[0]
        for ii in range(n):
            coup_av[ii,ii] = 0.0
        plot_couple(coup_av, figname='COUPLE_NA.png')

        stype = int(inp['SOCTYPE'])
        if (stype==2):
            coup = read_couple('SOTXT', ctype=1)
            coup_av = np.average(np.abs(coup), axis=0)
            n = coup_av.shape[0]
            for ii in range(n):
                coup_av[ii,ii] = 0.0
            plot_couple(coup_av, figname='COUPLE_SO.png')


    if (2 in which_plt) or (3 in which_plt):
        filshps = glob('SHPROP.*')
        shp, ksen, cw = data_proc(inp, pathD, filshps)

    if (2 in which_plt):
        plot_tden(shp, ksen, cw=cw, figname='TDEN.png')

    if (3 in which_plt):
        plot_tdpop(shp, cw, figname='TDPOP.png')

    if (4 in which_plt):
        plot_tdksen(pathD, emin=-2.0, emax=3.0, figname='TDKSEN.png')

    print("\nDone!\n")


###############################################################################


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


def loadData(inp, pathD='Data', atomsA=None, atomsB=None):

    whichK = 1
    nsw = int(inp['NSW'])
    ndigit = len(inp['NSW'])
    pathrun = inp['RUNDIR']

    fil = '{:0>{width}}'.format(1, width=ndigit)
    Ns = ReadNs(os.path.join(pathrun, fil))

    print('\nReading band energies from OUTCAR...\n')

    E, Ef = [], []
    for ii in range(nsw):

        fil = '{:0>{width}}'.format(ii+1, width=ndigit)
        path = os.path.join(pathrun, fil)
        energy, fermi = ReadEnergy(path, whichK, Ns)
        E.append(energy); Ef.append(fermi)

        if (ii+1) % (nsw//10) == 0:
            print('%6d OUTCARs have been read.'%(ii+1))

    E = np.array(E); Ef = np.array(Ef)

    stype = int(inp['SOCTYPE'])
    if stype==2:
        W = np.ones_like(E, dtype=float)
        n = int( E.shape[1] / 2 )
        W[:,n:] = -1.0

    else:
        print('\nReading weight from PROCAR...\n')

        wtype = 33
        # 31, 32, 33: spin weight for x, y and z directions,
        # respectively. (For VASP_ncl)

        W = []
        for ii in range(nsw):
            fil = '{:0>{width}}'.format(ii+1, width=ndigit)
            path = os.path.join(pathrun, fil)
            weight = ReadWeight(path, atomsA, atomsB, whichK, Ns, wtype)
            W.append(weight)

            if (ii+1) % (nsw//10) == 0:
                print('%6d PROCARs have been read.'%(ii+1))

        W = np.array(W)

    if not os.path.exists(pathD):
        os.mkdir(pathD)

    averW = np.average(W, axis=0)
    averE = np.average(E, axis=0) - np.average(Ef)
    average = np.vstack((averE,averW)).T

    path = os.path.join(pathD, 'energy.dat')
    np.savetxt(path, E, fmt='%10.6f')

    path = os.path.join(pathD, 'weight.dat')
    np.savetxt(path, W, fmt='%10.6f')

    path = os.path.join(pathD, 'fermi.dat')
    np.savetxt(path, Ef, fmt='%10.6f')

    path = os.path.join(pathD, 'average.dat')
    np.savetxt(path,average,fmt="%12.6f")

    print('\nBand energy information have been saved!')


def ReadNs(path):
    '''
    Read number of k-points, bands and ions from  PROCAR.
    '''

    infile = os.path.join(path, 'PROCAR')

    with open(infile) as f:
        for line in f:
            if 'k-points' in line:
                Ns = re.findall(r'\d+',line)
                #Ns = [nkpts, nbands, nions]
                break

    Ns = [int(x) for x in Ns]

    return Ns
    

def ReadEnergy(path, whichK, Ns):
    '''
    Read energies of states and fermi level from OUTCAR.
    '''
    
    infile = os.path.join(path, 'OUTCAR')
    out = [line for line in open(infile) if line.strip()]

    energy = []; fermi = []
    nkpts, nbands, nions = Ns
    jj = -nbands - 1
    for ii, line in enumerate(out):

        if ii in range(jj, jj+nbands):
            energy.append(float(line.split()[1]))
        elif 'band No.' in line:
            jj = ii + 1
        elif 'E-fermi' in line:
            fermi.append(float(line.split()[2]))

    nspin = int( len(energy)/nkpts/nbands )
    energy = np.array(energy).reshape(nspin, nkpts, nbands)
    energy = energy[:,whichK-1,:].flatten()
    fermi = np.array(fermi)

    return energy, fermi


def ReadWeight(path, atomsA, atomsB, whichK, Ns, wtype):
    '''
    Read weight of atomsA in atomsA & atomsB of each orbital from PROCAR.
    wtype = 1 or 31, 32, 33
    1: spatial weight, need to set parameters atomsA/B bellow.
    31, 32, 33: spin weight for x, y and z directions, respectively. (For VASP_ncl)
    '''

    infile = os.path.join(path, 'PROCAR')
    procar = [line for line in open(infile) if line.strip()]

    weight = []
    nkpts, nbands, nions = Ns

    if wtype==1:

        for line in procar:
            if not re.search('[a-zA-Z]',line):
                weight.append(float(line.split()[-1]))

        nspin = int( len(weight)/(nkpts*nbands*nions) )
        weight = np.array(weight).reshape(nspin,nkpts,nbands,nions)
        weight = weight[:,whichK-1,:,:]

        weightA = np.sum(weight[:,:,atomsA], axis=-1)
        weightB = np.sum(weight[:,:,atomsB], axis=-1)
        weight = weightA/(weightA+weightB)
        weight = weight.flatten()

    else:
        for line in procar:
            if 'tot ' in line:
                weight.append(float(line.split()[-1]))
        weight = np.array(weight).reshape(nkpts, nbands, 4)
        weight = weight[whichK-1,:,:]

        ispin = int( wtype%10 )
        totspin = np.linalg.norm(weight[:,1:], axis=-1) + 1.0e-5
        weight = weight[:, ispin] / totspin
        weight = weight.flatten()

    return weight


def data_proc(inp, pathD, filshps):

    if not filshps:
        print("\nERROR: SHPROP files are not found!\n")
        os._exit(0)

    path = os.path.join(pathD, 'energy.dat')
    energy = np.loadtxt(path)

    path = os.path.join(pathD, 'weight.dat')
    weight = np.loadtxt(path)

    path = os.path.join(pathD, 'fermi.dat')
    fermi  = np.loadtxt(path)

    Eref = np.average(fermi)

    stype = int(inp['SOCTYPE'])
    if (stype==1):
        bmin = int(inp['BMIN'])
        bmax = int(inp['BMAX'])
        bands = np.arange(bmin-1, bmax)
    elif (stype==2):
        bminU = int(inp['BMINU'])
        bmaxU = int(inp['BMAXU'])
        bandsU = np.arange(bminU-1, bmaxU)
        bminD = int(inp['BMIND'])
        bmaxD = int(inp['BMAXD'])
        bandsD = np.arange(bminD-1, bmaxD)
        nbands = int(inp['NBANDS'])
        bands = np.hstack((bandsU, bandsD + nbands))
        weight[:,:nbands] =  1.0
        weight[:,nbands:] = -1.0

    energy = energy[:, bands]
    weight = weight[:, bands]

    # if (weight.min() > 0): weight = weight * 2.0 - 1.0

    nsample = len(filshps)

    for ii in range(nsample):

        fil = filshps[ii]
        shp = np.loadtxt(fil)
        c = shp[:, 2:]

        if ii==0:
            nts = shp.shape[0]
            namdtime = shp[-1,0]
            potim = namdtime / nts

        it1 = int( float(fil.split('.')[-1]) / potim ) - 1
        it2 = it1 + nts

        if ii==0:
            shp_avg  = shp * 1.0
            ksen     = energy[it1:it2, :] * 1.0
            cw       = c * weight[it1:it2, :]
        else:
            shp_avg += shp
            ksen    += energy[it1:it2, :]
            cw      += c * weight[it1:it2, :]

    shp_avg /= nsample
    ksen    /= nsample
    cw      /= nsample

    ksen -= Eref
    shp_avg[:,1] = shp_avg[:,1] - Eref

    return shp_avg, ksen, cw


def plot_couple(coup, figname='COUPLE.png'):

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
    print("\n%s has been saved."%figname)


def plot_tden(shp, ksen, cw=None, figname='TDEN.png'):

    figsize_x = 4.8
    figsize_y = 3.2 # in inches
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    if (cw is None):
        cmap = 'hot_r'
        cpop = shp[:, 2:]
        cmin = 0.0; cmax = np.max(cpop)
        cmax = math.ceil(cmax*10)/10
    else:
        cmap = 'bwr'
        cpop = cw
        cmin = -1.0; cmax = 1.0
        aa = np.abs(cw.min())
        bb = np.abs(cw.max())
        if aa > bb:
            cmin = -aa; cmax = aa
        else:
            cmin = -bb; cmax = bb

    dotsize = 5
    namdtime = shp[-1,0]
    nts = shp.shape[0]
    nbands = shp.shape[1] -2
    norm = mpl.colors.Normalize(cmin,cmax)

    if (ksen.shape[1]!=nbands):
        print('\nNumber of ksen states doesn\'t match with SHPROP data!\n')
    T = np.tile(shp[:,0], nbands).reshape(nbands,nts).T
    sc = ax.scatter(T, ksen, s=dotsize, c=cpop, lw=0,
                    norm=norm, cmap=cmap)
    ax.plot(shp[:,1], 'r', lw=1, label='Average Energy')
    plt.colorbar(sc)

    ax.set_xlim(0,namdtime)
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Energy (eV)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)


def plot_tdpop(shp, cw, figname='TDPOP.png'):

    figsize_x = 4.8
    figsize_y = 3.2 # in inches
    namdtime = shp[-1,0]

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    time = shp[:,0]
    pop = np.sum(cw, axis=1) / 2.0 + 0.5

    lbs = ['Component 1', 'Component 2']
    ax.plot(time,   pop, 'r', lw=1.0, label=lbs[0])
    ax.plot(time, 1-pop, 'b', lw=1.0, label=lbs[1])

    ax.legend(loc=1)
    ax.set_xlim(0, namdtime)
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Population')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)


def plot_tdksen(pathD, emin, emax, potim=1.0, figname='TDKSEN.png'):

    path = os.path.join(pathD, 'energy.dat')
    energy = np.loadtxt(path)

    path = os.path.join(pathD, 'weight.dat')
    weight = np.loadtxt(path)

    path = os.path.join(pathD, 'fermi.dat')
    fermi  = np.loadtxt(path)

    Eref = np.average(fermi)
    energy -= Eref;

    avgE = np.average(energy, axis=0)

    figsize_x = 4.8
    figsize_y = 3.2 # in inches
    fig, ax = plt.subplots()
    fig.set_size_inches(figsize_x, figsize_y)
    mpl.rcParams['axes.unicode_minus'] = False

    E    = energy.T[(avgE>emin) & (avgE<emax)].T
    cpop = weight.T[(avgE>emin) & (avgE<emax)].T

    nsw = E.shape[0]
    nbs = E.shape[1]

    T = np.mgrid[0:nsw, 0:nbs][0] * potim

    cmap = 'bwr'
    cmin = -1; cmax = 1
    norm = mpl.colors.Normalize(cmin,cmax)

    ax.scatter(T, E, c=cpop, s=1, lw=0.0, alpha=1.0,
               norm=norm, cmap=cmap)

    ax.set_ylim(emin, emax)
    ax.set_xlim(0, nsw*potim)
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Energy (eV)')

    plt.tight_layout()
    plt.savefig(figname, dpi=400)
    print("\n%s has been saved."%figname)


if __name__=='__main__':
    main()
