import numpy as np
import pygimli as pg
from pygimli.physics.ert.importData import (importAsciiColumns,
                                            importRes2dInv,
                                            importData)


def importTDIPdata(filename, verbose=False):
    """Read in TDIP data.

    Supported formats
    -----------------
    TXT - ABEM LS or Syscal Pro Ascii output
    BIN - Syscal binary format
    GDD - GDD format
    2DM - Ares II Ascii format
    AMP - ABEM SAS4000 format
    TX2 - Aarhus Workbench data
    DIP - AarhusInv (processed) data
    DAT - Res2dInv format

    Returns
    -------
    data : DataContainerERT
        data container
    MA : numpy.array (ngates, ndata)
        spectral chargeability matrix
    t : iterable
        time vector
    header : dict
        dictionary with supporting information
    """
    header = {}
    if isinstance(filename, str):
        ext = filename[filename.rfind('.')+1:]
        if ext.lower() in ['txt', 'tx2', 'gdd']:
            data, header = importAsciiColumns(
                filename, verbose=verbose, return_header=True)
        elif ext.lower() in ['dat', 'res2dinv']:
            data, header = importRes2dInv(filename, verbose=verbose,
                                          return_header=True)
        # elif ext.lower() == 'amp': # outdated?
        #     return importABEM(filename, return_all=True)
        elif ext.lower() == 'dip':
            return importDIP(filename, return_all=True)
        elif ext.lower() == '2dm':
            return importAres2(filename, return_all=True)
        else:
            data = importData(filename)
    elif isinstance(filename, pg.DataContainer):
        data = filename
    else:
        raise TypeError("Cannot use this type:"+type(filename))

    ipkey = ''
    testkeys = ['IP_#{}(mV/V)', 'M{}', 'ip{}']
    for key in testkeys:
        if data.exists(key.format(1)):
            ipkey = key
    MA = []
    i = 1
    while data.exists(ipkey.format(i)):
        ma = data(ipkey.format(i))
        if max(ma) <= 0:
            break

        MA.append(ma)
        i += 1

    MA = np.array(MA)
    t = np.arange(MA.shape[0]) + 1  # default: t is gate number
    if 'ipGateT' in header:  # Syscal
        t = header['ipGateT'][:-1] + np.diff(header['ipGateT'])/2
    elif data.exists('Gate1') and data.exists('Ngates'):
        ngates = int(data('Ngates')[0])
        dt = np.array([data('Gate'+str(i+1))[0] for i in range(ngates)]) * 1e-3
        delay = 0.001  # hard-coded as lost in import
        t = np.cumsum(dt) - dt/2 + delay
        header['ipGateT'] = np.cumsum(np.hstack((0, dt))) + delay
    elif 'IP_WindowSecList' in header:  # new Terrameter LS 2
        dt = np.array([float(sdt) for sdt in
                       header["IP_WindowSecList"].split()])
        gatet = np.cumsum(dt)
        t = (gatet[:-1] + gatet[1:]) / 2
        header['ipGateT'] = gatet

    testkeys = ['TM{}']
    for key in testkeys:
        if data.exists(key.format(1)):
            tt = [data(key.format(i+1))[0] for i in range(
                MA.shape[0])]
            dt = np.array(tt)
            if sum(dt) > 100:
                dt *= 0.001

            delay = 0.0
            if data.exists(key.format(0)):
                delay = data(key.format(0))[0]
                if delay > 1:
                    delay *= 0.001

            t = np.cumsum(dt) - dt/2 + delay
            header['ipGateT'] = np.cumsum(np.hstack((delay, dt)))

    return data, MA, t, header

def importDIP(filename, verbose=True, return_header=False, return_all=False):
    """Import ERT data from Ares II (GF Brno) multielectrode instrument."""
    header = {}
    data = pg.DataContainerERT()
    with open(filename) as fid:
        lines = fid.readlines()
        tokens = None
        for i, line in enumerate(lines):
            if line.find('NumElect') > 0 or line.find('rho') > 0:
                tokens = line[1:].replace('IP ', 'IP_').split()
            if tokens is not None and line[0] != '/':
                break

        A = np.genfromtxt(lines[i:-1], names=tokens)

        data = pg.DataContainerERT()
        ndata = len(A)
        # ux = np.unique(np.hstack((A['XA'], A['XB'], A['XM'], A['XN'])))
        XX = np.concatenate([A["X"+tok] for tok in ["A", "B", "M", "N"]])
        ux, ix = np.unique(XX, return_index=True)
        if "ZA" in A.dtype.names:
            ZZ = np.concatenate([A["Z"+tok] for tok in ["A", "B", "M", "N"]])
            uz = ZZ[ix]
        else:
            uz = np.zeros_like(ux)
        for xi, zi in zip(ux, uz):
            data.createSensor([xi, 0, zi])

        for i in range(ndata):
            nn = [int(np.nonzero(ux == A['X'+s][i])[0]) for s in
                  ['A', 'B', 'M', 'N']]
            data.createFourPointData(i, *nn)

        data.set('valid', A['InUse'])
        # FW: Avoid f-strings to keep Py 3.5 compatibility
        nst = ['IP_{}'.format(i) for i in range(1, 10)]
        i = 10
        while 'IP{}_data'.format(i) in tokens:
            nst.append('IP{}'.format(i))
            i += 1

        MAall = np.column_stack([A[nsi+'_data'] for nsi in nst]).T
        VA = np.column_stack([A[nsi+'_inUse'] for nsi in nst]).T
        MA = np.ma.MaskedArray(MAall, np.isclose(VA, 0))

        t = np.array([A[nsi+'_center'][0] for nsi in nst])
        dt = np.array([A[nsi+'_width'][0] for nsi in nst])
        header['ipDT'] = dt
        header['delay'] = t[0] - dt[0]/2
        header['ipGateT'] = np.hstack((t-dt/2, t[-1]+dt[-1]/2))
        data.set('rhoa', A['Rho'])

        if return_all:  # the variant needed for the TDIP class
            return data, MA, t, header
        elif return_header:  # might become the default for all
            return data, header
        else:  # just the data container as currently used by importData
            return data  # unless all importers support returning a header


def importAres2(filename, verbose=True, return_header=False, return_all=False):
    """Import ERT data from Ares II (GF Brno) multielectrode instrument."""
    header = {}
    data = pg.DataContainerERT()
    with open(filename) as fid:
        lines = fid.readlines()
        for i, line in enumerate(lines):
            lines[i] = re.sub("\*[0-9]", "", line.rstrip())

        for i, line in enumerate(lines):
            if line.startswith('C1'):
                break
            if ':' in line and len(line) <= 80:
                sp = line.split(':')
                header[sp[0]] = sp[1]

        # first remove *1 in elecrode spacings & find out what it means
        cols = np.genfromtxt(lines[i:], names=True, delimiter='\t',
                             autostrip=True)

        if verbose:
            print(header)
        if 'Electrode distance' in header:
            dx = float(header['Electrode distance'].split('m')[0])
        elif 'Distance' in header:
            dx = float(header['Distance'].split('m')[0])
        else:
            raise Exception('Electrode distance cannot be determined.')

        if 'Profile length' in header:
            plen = float(header['Profile length'].split('m')[0])
        elif 'Length' in header:
            plen = float(header['Length'].split('m')[0])

        nel = int(plen/dx)+1
        for i in range(nel):
            data.createSensor([i*dx, 0])

        data.resize(len(cols))
        tokIn = ['C1', 'C2', 'P1', 'P2']
        tokOut = ['a', 'b', 'm', 'n']
        names = cols.dtype.names
        for i in range(4):
            tok = tokIn[i]+'el'
            if tok not in names:
                tok = tokIn[i]+'_el'
            col = cols[tok]
            col[np.isinf(col)] = -1
            data.set(tokOut[i], col)

        if 'ImA' in names:
            data.set('i', cols['ImA'] / 1000)
        elif 'I_mA' in names:
            data.set('i', cols['I_mA'] / 1000)
        if 'UmV' in names:
            data.set('u', cols['UmV'] / 1000)
        elif 'V_mV' in names:
            data.set('u', cols['V_mV'] / 1000)
        if 'Stdev' in names:
            data.set('err', cols['Stdev'] / 100)
        elif 'Stdev_' in names:
            data.set('err', cols['Stdev_'] / 100)
        if 'AppResOhmm' in names:
            data.set('rhoa', cols['AppResOhmm'])
        elif 'AppRes_Ohmm' in names:
            data.set('rhoa', cols['AppRes_Ohmm'])

        MA = []
        for i in range(33):
            sti = "IP" + str(i+1)
            if sti in cols.dtype.names:
               MA.append(cols[sti] * 10)  # in % instead of mV/V

        t = []
        delay = 0.005  # fixed in instrument
        if "IP windows" in header:
            dt = np.fromstring(header["IP windows"][:-2], dtype=float,
                               sep="\t") * 0.001  # ms
            tGate = delay + np.hstack((0, np.cumsum(dt)))
            header["ipGateT"] = tGate
            t = tGate[:-1] + dt/2

        data.markValid(data('rhoa') > 0)
        if return_all:  # the variant needed for the TDIP class
            return data, np.array(MA), t, header
        elif return_header:
            return data, header
        else:
            return data

