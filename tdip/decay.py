import numpy as np
import matplotlib.pyplot as plt
import pygimli as pg
from pygimli.frameworks import MarquardtInversion
from pygimli.physics.SIP import SIPSpectrum
from .modelling import ColeColeTD, MultiDebyeTDModelling, CCTDModelling


class Decay(object):
    """Decay (time-domain decaying function)."""

    def __init__(self, t=None, v=None, tau=None):
        """Initialize decay (and load).

        Parameters
        ----------
        t : iterable
            time vector
        v : iterable
            value (voltage) vector
        tau : iterable
            vector of characteristic relaxation times
        """
        if isinstance(t, str):  # a file
            t, v = np.genfromtxt(t)

        self.t = t
        self.v = v
        self.tau = tau

    def __repr__(self):  # for print function
        """Readable representation of the class."""
        return "Decay with {:d} points ({:.3f}..{:.3f}s)".format(
            len(self.t), min(self.t), max(self.t))

    def removeMaskedData(self):
        """Remove masked data from decay."""
        if isinstance(self.v, np.ma.masked_array):
            self.t = self.t[~self.v.mask]
            self.v = self.v.data[~self.v.mask]

    def noMask(self):
        pg.deprecated("noMask", "removeMaskedData")
        return self.removeMaskedData()

    def filter(self, tmin=0, tmax=100, vmin=0, vmax=10000):
        """Filter data (remove data outside criteria).

        Parameters
        ----------
        tmin : float
            minimum time
        tmax : float
            maximum time
        vmin : float
            minimum value
        vmax : float
            maximum value
        """
        self.noMask()
        ind = ((self.t >= tmin) & (self.t <= tmax) &
               (self.v > vmin) & (self.v < vmax))
        self.t = self.t[ind]
        self.v = self.v[ind]

    @property
    def tau(self):
        """Time constants vector."""
        return self._tau

    @tau.setter
    def tau(self, tau):
        """Set time constants vector."""
        self._tau = tau
        self.initFopInv()

    def initFopInv(self):
        """Initialize forward operator and inversion."""
        self.fDD = MultiDebyeTDModelling(self.t, self._tau, verbose=False)
        self.invDD = pg.Inversion(fop=self.fDD, verbose=False)
        self.tLog = pg.trans.TransLog()
        self.invDD.transModel = self.tLog
        self.invDD.transData = self.tLog

    def decompose(self, v=None, tau=None, error=0.03, show=False, **kwargs):
        """Debye decomposition.

        Parameters
        ----------
        v : iterable [optional]
            decay values to decompose (otherwise self.v is used)
        tau : iterable [optional]
            time constant vector to use (otherwise try to use self.tau)
        error : float | array
            error model for weighting the inversion
        show : bool
            plot data with forward response
        startModel : float | array
            starting model (otherwise taken from decay energy)
        **kwargs : dict
            inversion/regularization options to be passed to inversion instance
        """
        if v is None:
            if isinstance(self.v, np.ma.masked_array):
                v = np.array(self.v.data)
                v[self.v.mask] = -0.001
            else:
                v = self.v

        if tau is not None:
            self.tau = tau

        if self.tau is None:
            self.tau = self.t

        if kwargs.pop("nnls", True):
            from scipy.optimize import nnls
            G = np.zeros([len(self.t), len(self.tau)])
            for i, tau in enumerate(self.tau):
                G[:, i] = np.exp(-self.t/tau)

            self.modelDD, *_ = nnls(G, v)
            self.invDD.response = G.dot(self.modelDD)
        else:
            kwargs.setdefault("startModel", sum(v) / len(self.tau))
            if isinstance(error, float):
                error = np.ones_like(v) * error

            error[v <= 1e-5] = 10000.
            v[v <= 1e-5] = 0.0001
            v = np.abs(v)

            # self.INV.setRegularization(limits=[])
            self.modelDD = self.invDD.run(v, error, **kwargs)
            if show:
                self.invDD.echoStatus()

        if show:
            return self.showModel()

    def logMeanTau(self):
        """Return log-mean tau value."""
        return np.exp(np.sum(np.log(self.fDD.tau)*self.modelDD) /
                      np.sum(self.modelDD))

    def invert(self, **kwargs):
        """Invert for Cole-Cole parameters."""
        self.fCC = CCTDModelling(self.t)
        self.invCC = MarquardtInversion(fop=self.fCC)
        startModel = kwargs.pop("startModel", np.array([self.v[0], 0.3, 0.25]))
        self.modelCC = self.invCC.run(self.v, np.ones_like(self.v)*0.01,
                                      startModel=startModel, **kwargs)

    def showModel(self, ax=None, **kwargs):
        """Show Debye decomposition model."""
        if ax is None:
            fig, ax = plt.subplots()

        ax.semilogx(self.tau, self.modelDD, **kwargs)
        ax.grid(True)
        ax.set_xlabel(r"$\tau$ (s)")
        ax.set_ylabel(r"$m$ (-)")
        return ax

    def show(self, v=None, ax=None, xScale=None, yScale=None, **kwargs):
        """Show decay.

        Parameters
        ----------
        ax : matplotlib axes
            axes object
        xScale/yScale : str
            x and y scale: ['linear', 'log', 'symlog', 'logit']
        kwargs : dict
            keyword arguments to be passed to the plot
        """
        if ax is None:
            fig, ax = plt.subplots()

        v = v or self.v
        if kwargs.pop('dp', False):
            v = self.dp()

        ax.plot(self.t, v*1000, **kwargs)
        if xScale:
            ax.set_xscale(xScale)
        if yScale:
            ax.set_yscale(yScale)

        ax.grid(True)
        ax.set_xlabel("t (s)")
        ax.set_ylabel("v (mV/V)")
        return ax

    def showAll(self, details=False, **kwargs):
        """Show data with model response.

        Parameters
        ----------
        details : bool [False]
            plot individual Debye decays
        *kwargs : dict
            keywords passed to plot function
        """
        kwargs.setdefault("marker", "*")
        kwargs.setdefault("color", "blue")
        kwargs.setdefault("linestyle", " ")
        kwargs.setdefault("label", "data")
        ax = self.show(**kwargs)
        if hasattr(self, "invDD"):
            response = kwargs.pop("response", self.invDD.response)
            ax.plot(self.t, response*1000, "r-", label="Debye")
            if details:
                for i in range(len(self.fDD.tau)):
                    ax.plot(self.t, self.fDD.T.col(i)*self.modelDD[i]*1000,
                            "r--")
        if hasattr(self, "invCC"):
            ax.plot(self.t, self.invCC.response*1000, "g-", label="Cole-Cole")

        ax.legend()
        # ax.set_ylim(0.001, max(self.v)*1.1)
        return ax

    def showOld(self, ax=None, **kwargs):
        """Show decay."""
        if ax is None:
            fig, ax = plt.subplots()

        ax.loglog(self.t, self.v, **kwargs)
        ax.grid(True)
        return ax

    def simulate(self, m=1.0, tau=None, c=1, t=None, **kwargs):
        """Generate Debye (c=1) or Cole-Cole (c<1) model.

        Parameters
        ----------
        tau : float | iterable
            Time constant (relaxation time)
        m : float | iterable
            chargeability, needs to match time constant
        c : float
            relaxation exponent (for Cole-Cole model)
        t : iterable
            time discretization (gate midpoint) vector
        """
        tau = tau or self.tau
        if t is None:
            t = self.t
        if c == 1:
            if hasattr(tau, "__iter__") and hasattr(m, "__iter__"):
                u = np.zeros_like(t)
                for taui, mi in zip(tau, m):
                    u += np.exp(-t / taui) * mi
            else:
                u = np.exp(-t / tau) * m
        else:
            u = ColeColeTD(t=t, m=m, tau=tau, c=c)

        return u

    def dp(self):
        """Compute differential polarisability."""
        logT = np.log(self.t)
        dlogT = np.diff(logT)
        dV = np.diff(self.v)
        dVdlogT = - dV / dlogT
        dp = np.hstack((dVdlogT[0], (dVdlogT[:-1]+dVdlogT[1:])/2, dVdlogT[-1]))
        return dp


    def convertDDToSpectrum(self, f=None, rho=1.):
        """Convert Debye decomposition into frequency-domain spectrum.

        Parameters
        ----------
        f : iterable [logspace(-3, 3, 41)]
            frequency vector
        rho : float [1]
            DC resistivity/impedance
        """
        if self.modelDD is None:
            self.decompose()

        if f is None:
            f = np.logspace(-3, 3, 41)

        T, W = np.meshgrid(self.tau, f * 2. * np.pi)
        A = 1 - 1. / (W*T * 1j + 1)
        Z = (1. - A.dot(self.modelDD)) * rho

        return SIPSpectrum(f=f, amp=np.abs(Z), phi=-np.angle(Z))


if __name__ == "__main__":
    t = np.logspace(-2, np.log10(4), 12)
    print(t)
    mydecay = Decay(t=t)
    # %%
    mydecay.v = mydecay.simulate(tau=0.2, c=0.5)
    mydecay.show(dp=True, xScale='log', yScale='log')
