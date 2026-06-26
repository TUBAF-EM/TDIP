# Time-Domain Induced-Polarization (TDIP) data processing

Formerly part of BERT repository, but since 2026 independent.
The code is based upon pyGIMLi (Rücker et al. 2017) and its ERT module based on Günther et al. (2006).
Whereas the ERT module only includes only single-frequency or chargeability inversion, TDIP is designed for analysing the spectral induced polarization (SIP) in the time domain. For frequency domain SIP, please use [FDIP](https://github.com/TUBAF-EM/FDIP)
For details and examples, we refer to Martin et al. (2020) where it was first used, and the accompagnying Zenodo data sets with codes.
Further application was already published by Rossi et al. (2018) or Bazin et al. (2018).

To use the package and keep updated, just

1. `git clone` the repository and go to its location with bash, PowerShell etc.
2. install the code editable by `pip install -e .`
3. Update at any later time by `git pull` (no further step needed)

If you use UV, the simplest way is to create a virtual environment by typing `uv run` in the project folder.
The venv is created in the main folder and if you open the folder in VSCode it is chosen as default environment.

## References

* Martin, T., Günther, T., Orozco, A.F. & Dahlin, T. (2020): Evaluation of spectral induced polarization field measurements in time and frequency domain, J. Appl. Geophys. 180, 104141, [doi:10.1016/j.jappgeo.2020.104141](https://doi.org/10.1016/j.jappgeo.2020.104141).
* Rücker, C., Günther, T., Wagner, F.M. (2017): pyGIMLi: An open-source library for modelling and inversion in geophysics, Computers & Geosciences 109, 106-123, [doi:10.1016/j.cageo.2017.07.011](http://doi.org/10.1016/j.cageo.2017.07.011).
* Rossi, M., Dahlin, T., Olsson, P.-I. & Günther, T. (2018): Data acquisition, processing and filtering for reliable 3D resistivity and time-domain induced polarization tomography in an urban area: field example of Vinsta, Stockholm, Near Surface Geophysics 16(3), 220-229, [doi:10.3997/1873-0604.2018014](https://doi.org/10.3997/1873-0604.2018014)
* Bazin, S., Lysdahl, A.K., Viezzoli, A., Günther, T., Anschütz, H., Scheibz, J., Pfaffhuber, A.A., Radic, T. & Fjermestad, H. (2018): Resistivity and chargeability survey for tunnel investigation: a case study on toxic black shale in Norway. Near Surface Geophysics 16(1), 1-11, [doi:10.3997/1873-0604.2017036](https://doi.org/10.3997/1873-0604.2017036)
* Günther, T. & Martin, T. (2016): Spectral two-dimensional inversion of frequency-domain induced polarisation data from a mining slag heap. Journal of Applied Geophysics 135, 436-448, [doi:10.1016/j.jappgeo.2016.01.008](https://doi.org/10.1016/j.jappgeo.2016.01.008).