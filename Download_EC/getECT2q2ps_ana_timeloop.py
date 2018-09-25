#!/usr/bin/env python
import calendar
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

def retrieve_interim():
    """       
    A function to demonstrate how to iterate efficiently over several years and months etc for a particular interim_request.      
       Change the variables below to adapt the iteration to your needs. 
       You can use the variable 'target' to organise the requested data in files as you wish.
       In the example below the data are organised in files per year. (eg "t2m_d2m_ps_anaEI_1984.nc")
    http://apps.ecmwf.int/codes/grib/param-db
    134 is mean surface pressure
    """
    yearStart = 1979
    yearEnd = 2018

    for year in list(range(yearStart, yearEnd + 1)):
        startDate = '%04d%02d%02d' % (year, 1, 1)
        #numberOfDays = calendar.monthrange(year, month)[1]
        lastDate = '%04d%02d%02d' % (year, 12, 31)
        target = "t2m_d2m_ps_anaEI_%04d.nc" % (year)
        requestDates = (startDate + "/TO/" + lastDate) 
        print('extracting',requestDates)
        print('going to', target)
        interim_request(requestDates, target)

def interim_request(requestDates, target):
    """       
        An ERA interim request for analysis surface Td, T, ps data.
        Change the keywords below to adapt it to your needs.
        (eg to add or to remove  levels, parameters, times etc)
        Request cost per day is 112 fields, 14.2326 Mbytes 
    """

    server.retrieve({
        "class": "ei",
        "dataset": "interim",
        "date": requestDates, 
        "expver": "1",
        "grid": "0.75/0.75",
        "levtype": "sfc",
        "param": "134.128/167.128/168.128",
        "step": "0",
        "stream": "oper",
        "time": "00:00:00/06:00:00/12:00:00/18:00:00",
        "type": "an",
        "area": "72.75/-1.5/57/33", 
        "format": "netcdf",
        "target": target
    })

if __name__ == '__main__':

    retrieve_interim()
