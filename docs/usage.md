The tdip package holds some functions for importing data, and for modelling IP effects, all based upon [pyGIMLi](https://pygimli.org).

Main entry, however, is the class TDIP, to be imported by

```python
from tdip import TDIP
```

It can be initialized empty or directly with a data file

```python
pro = TDIP("datafile.ohm")
```

Valid formats are:

* Syscal Pro export file (*.txt)
* ABEM TXT export file (*.txt or raw time series)
* Syscal Pro binary file (*.bin)
* GDD format (*.gdd)
* Ares II format (*.2dm)
* Aarhus Workbench processed data (*.tx2 and *.dip)
* res2dinv data

You can also use `pro.loadData(filename)`.

As a result, the instance holds an ERT data file `pro.data`, which can be used and processed with the pyGIMLi ERT module, and a matrix of apparent chargeabilities `pro.MA`.

For most formats, it should also load the gate times `pro.t`, if not, they can be set by `pro.setGates()`.

## Preprocessing

`pro.showRhoa()` shows the apparent resistivity. If it is not there yet, it can be created by computing the geometric factor, either analytically or numerically, by `ensureRhoa()`. This should be done after changing the electrode positions after loading.

`.fitDataDecays()` fits every decay by a simple Debye (exponential decay) or Cole-Cole (possibly also CPA) model to have a rough guess of the quality of the decay.

The function `.showData()` can show every property like `'a'` or `'rhoa'`, including he fitted `'m0'`, `'tau'` or the fit quality `'fit'`.

You can show the apparent chargeability for individual gates by `.showApparentChargeability()`.
The function `pro.generateDataPDF()` generates a multipage PDF with all included data fields.


Decays can be shown by `pro.showDecay()` passing an AB or MN combination (or both).
The function `pro.generateDecayPDF()` generates a multipage PDF with all decays grouped for the same AB.
You can also work with individual decays after
`decay=pro.getDecay(abmn=...)`

For further manipulation, the functino `getDataIndex(abmn=[...])` can help you identifying specific numbers.

Most important preprocessing tool is `pro.filter()`.
It works either along the data sequence axis (like `data.remove()`) or along the time axis and removes complete data (including all gates) or gates (for all data).
Important keywords are

```python
rmin=0, rmax=1e9, emax=1e9, # apparent resistivity and error
kmax=1e6, electrode=None, forward=False,
umin=0, umax=9e9, nr=[],
a=None, b=None, m=None, n=None, ab=None, mn=None, corrSID=1,
tmin=0, tmax=1e9, # gate times
fitmax=1e9, m0max=None, m0min=None, taumin=-9e99, taumax=9e99, # fitted decays
```

If you don't want to remove the whole data point, but only disregard some erroneous points (like outliers), you can mask these using `.mask()` with the keywords `mamin` and `mamax`.
You can also mask the whole decay by `filter(..., mask=True)` in order to not to loose the apparent resistivity information.

## Inversion

First, you start with inverting the apparent resistivity through
`pro.invertRhoa()`

This triggers an `ERTManager()` and you can pass all options.

You can invert a single apparent chargeability by `pro.invertMa()` specifying its number (0 for the fitted m0).

To invert all gates sequentially, use
`pro.individualInversion()`
or simultaneously, use
`pro.simultaneousInversion()`

There is a special inversion into Debye models for every cell.
`pro.invertDebye()`

Class instances can be accessed by `pro.ERT`, `pro.invIP`

## Postprocessing

As a result of apparent chargeability inversion, you have a matrix of chageabilities `pro.M` for every model cell, matching the inversion mesh `pro.pd`.

`.getModelDecay()` returns a decay class instance for a number or position (using `getCellID()`).


`.fitModelDecays()`

`fitModelDecays()` can be used to fit a Cole-Cole model to every model decay, and display `m0`, `tau` and `c` as sections.

`showColeColeResults()`

`saveFigures()`

`saveResults()`

`loadResults()`

so that you can do postprocessing independent of inversion by

```python
pro = TDIP(filename)
pro.loadResults()
...
```

## Working in 3D

While visualization is designed for 2D profiles, numerical computations should work in 3D as well, provided a 3D mesh is created for the ERT manager.

## Conversion to FD

`convertToFD()` converts the whole data set to a frequency-domain instance `FD` (package FDIP required).