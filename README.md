# HySN
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1435010.svg)](https://doi.org/10.5281/zenodo.1435010)

A simple python script to create a high resolution HYbrid SeNorge dataset of daily near surface humidity and incident shortwave and longwave radiation created by merging reanalysis data (Era-interim) with a national 1-by-1 km gridded temperature dataset (SeNorge). The assumptions and methods used to downscale humidity and longwave radiation are similar to those used in WATCH/PGMFD/NLDAS-1. The script will run on a desktop computer, and as it is currently configured might take about half a day to produce 97 billion data points compressed to about 45 GB of data consisting of daily estimates of specific humidity, surface pressure, incident shortwave radiation, and incident longwave radiation between 1979 and 2016. 

The current version of the files from 1979-2016 is stored at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1434802.svg)](https://doi.org/10.5281/zenodo.1434802)

### Monthly means of humidty (q<sub>2</sub>), and incident shortwave (SW) and longwave (LW) radiation in 1984: 
![monthlyhysn_q2](https://user-images.githubusercontent.com/23070665/46022268-8ddccd80-c0e2-11e8-9ffa-08a782015ce0.png)
![monthlyhysnsw](https://user-images.githubusercontent.com/23070665/46022565-2a06d480-c0e3-11e8-98cc-db29e7a0b2f0.png)
![monthlyhysnlw](https://user-images.githubusercontent.com/23070665/46022934-d8ab1500-c0e3-11e8-86b3-26da299ffe0c.png)
