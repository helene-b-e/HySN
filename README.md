# HySN
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1435010.svg)](https://doi.org/10.5281/zenodo.1435010)

A simple python script to create a high resolution HYbrid SeNorge dataset of daily near surface humidity and incident shortwave and longwave radiation created by merging reanalysis data (Era-interim) with a national 1-by-1 km gridded temperature dataset (SeNorge). The assumptions and methods used to downscale humidity and longwave radiation are similar to those used in WATCH/PGMFD/NLDAS-1. The script will run on a desktop computer, and as it is currently configured might take about half a day to produce 97 billion data points compressed to about 45 GB of data consisting of daily estimates of specific humidity, surface pressure, incident shortwave radiation, and incident longwave radiation between 1979 and 2016. 

The current version of the files from 1979-2017 is stored at Zenodo:[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1970170.svg)](https://doi.org/10.5281/zenodo.1970170). Monthly means are also avaliable from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1993870.svg)](https://doi.org/10.5281/zenodo.1993870)

### Monthly means of humidty (q<sub>2</sub>), and incident shortwave (SW) and longwave (LW) radiation (1979-2017): 
![monthlyhysn_q2](https://user-images.githubusercontent.com/23070665/49574679-3f00af00-f941-11e8-9f85-a539a530cc51.png)
![monthlyhysnsw](https://user-images.githubusercontent.com/23070665/49574711-4de76180-f941-11e8-88af-eb0ee96eef8f.png)
![monthlyhysnlw](https://user-images.githubusercontent.com/23070665/49574715-50e25200-f941-11e8-8c9b-554a090f0fac.png)
