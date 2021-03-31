# StratoVari
Here we conduct a series of large ensemble simulations (300 members each) to explore the robustness of the stratospheric response to Arctic sea-ice loss. We follow the new PAMIP protocol using the SC-WACCM4 and E3SM models.

###### Under construction... ```[Python 3.7]```

## Contact
Zachary Labe - [Research Website](http://sites.uci.edu/zlabe/) - [@ZLabe](https://twitter.com/ZLabe)

## Description

+ ```Scripts/```: Main [Python](https://www.python.org/) scripts/functions used in data analysis and plotting
+ ```requirements.txt```: List of environments and modules associated with the most recent version of this project. A Python [Anaconda3 Distribution](https://docs.continuum.io/anaconda/) was used for our analysis. All AGCM experiments were processed through resources on CISL's [Cheyenne](https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne) supercomputer. Tools including [NCL](https://www.ncl.ucar.edu/), [CDO](https://code.mpimet.mpg.de/projects/cdo), and [NCO](http://nco.sourceforge.net/) were also used for initial data manipulation.

## Data
+ Polar Amplification Model Intercomparison Project (PAMIP) : [[DATA]](https://pcmdi.llnl.gov/CMIP6/ArchiveStatistics/esgf_data_holdings/PAMIP/index.html)
    + Smith, D. M., Screen, J. A., Deser, C., Cohen, J., Fyfe, J. C., García-Serrano, J., ... & Peings, Y. (2019). The Polar Amplification Model Intercomparison Project (PAMIP) contribution to CMIP6: investigating the causes and consequences of polar amplification. Geoscientific Model Development, 12(3), 1139-1164. doi:10.5194/gmd-12-1139-2019 [[Publications]](https://www.geosci-model-dev.net/12/1139/2019/)
+ Energy Exascale Earth System Model (E3SMv1) : [[DATA]](https://e3sm.org/data/)
    + Golaz, J. C., Caldwell, P. M., Van Roekel, L. P., Petersen, M. R., Tang, Q., Wolfe, J. D., ... & Baldwin, S. A. (2019). The DOE E3SM coupled model version 1: Overview and evaluation at standard resolution. Journal of Advances in Modeling Earth Systems, doi:10.1029/2018MS001603 [[Publication]](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018MS001603)
+ Whole Atmosphere Community Climate Model (WACCM4) : [[DATA]](http://www.cesm.ucar.edu/working_groups/Whole-Atmosphere/code-release.html)
    + Marsh, D. R., M. J. Mills, D. E. Kinnison, J.-F. Lamarque, N. Calvo, and L. M. Polvani, 2013: Climate change from 1850 to 2005 simulated in CESM1 (WACCM). J. Climate, 26, 7372–7391, doi:10.1175/JCLI-D-12-00558.1 [[Publication]](http://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-13-00255.1)
    + Smith, K. L., R. R. Neely, D. R. Marsh, and L. M. Polvani (2014): The Specified Chemistry Whole Atmosphere Community Climate Model (SC-WACCM), Journal of Advances in Modeling Earth Systems, 6(3), 883–901, doi:10.1002/2014MS000346 [[Publication]](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014MS000346)


## Publications
[1] Peings, Y., **Z.M. Labe**, and G. Magnusdottir (2021), Are 100 ensemble members enough to capture the remote atmospheric response to +2°C Arctic sea ice loss? <em>Journal of Climate</em>, DOI:10.1175/JCLI-D-20-0613.1 [[HTML]](https://journals.ametsoc.org/view/journals/clim/aop/JCLI-D-20-0613.1/JCLI-D-20-0613.1.xml)[[BibTeX]](https://sites.uci.edu/zlabe/files/2021/02/PeingsLabeMagnusdottir_ENS2021_JCLI_BibTeX.pdf)

## Conferences
+ [1] Magnusdottir, G., Y. Peings, and **Z.M. Labe**. Response to sea-ice loss under the Polar Amplification MIP protocol in extended ensembles of simulations, *2019 American Geophysical Union Annual Meeting*, San Francisco, CA (Dec 2019). [[Abstract]](https://agu.confex.com/agu/fm19/meetingapp.cgi/Paper/553470)
+ [2] **Labe, Z.M.**, Y. Peings, and G. Magnusdottir. Detection of Signal in the Large-Scale Circulation Response to Arctic Sea-Ice Decline, *33rd Conference on Climate Variability and Change*, Boston, MA (Jan 2020). [[Abstract]](https://ams.confex.com/ams/2020Annual/meetingapp.cgi/Paper/367289)
+ [3] Peings, Y., **Z.M. Labe**, and G. Magnusdottir. Are 100 ensemble members enough to capture the remote atmospheric response to +2°C Arctic sea ice loss?, *2020 American Geophysical Union Annual Meeting*, Virtual Conference (Dec 2020). [[Abstract]](https://agu.confex.com/agu/fm20/meetingapp.cgi/Paper/685885)
+ [4] Peings, Y., G. Magnusdottir., and **Z.M. Labe**. Are 100 ensemble members enough to capture the remote atmospheric response to +2°C Arctic sea ice loss?, *34th Conference on Climate Variability and Change*, Virtual Conference (Jan 2021). [[Abstract]](https://ams.confex.com/ams/101ANNUAL/meetingapp.cgi/Paper/382352)
+ [5] Peings, Y., **Z.M. Labe**, and G. Magnusdottir. Influence of internal variability: how to ensure results are robust? *Polar Amplification Model Intercomparison (PAMIP) Virtual Workshop*, Virtual Conference (Mar 2021).
