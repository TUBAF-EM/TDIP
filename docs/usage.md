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

