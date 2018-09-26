import matplotlib.pyplot as plt
import xarray as xr
import calendar
from matplotlib import rcParams
datadir="data/"
EIpath=datadir+'ERAI/'
savedir=datadir+"HySN/"
yr='1984'
#map_proj = ccrs.UTM('33N')
Q2SN=xr.open_dataset(savedir+'HySN_Near_Surface_Specific_Humidity_'+str(yr)+'.nc', chunks={'time': 50})
Q2month=Q2SN['huss'].groupby('time.month').mean('time')*1000.
Q2month.name='$q_2$ [g/kg]'
p = Q2month.plot(col='month',cmap=plt.cm.gist_heat, col_wrap=4)
i=1
for ax in p.axes.flat:
    mont=calendar.month_abbr[i]
    ax.set_title('$q_2$  '+mont)
    ax.axis('off')
    i=i+1

plt.savefig('MonthlyHySN_Q2.png')
plt.show()


SWSN=xr.open_dataset(savedir+'HySN_Surface_Downwelling_Shortwave_Radiation_'+str(yr)+'.nc', chunks={'time': 30})
SWmonth=SWSN['rsds'].groupby('time.month').mean('time')
SWmonth.name='$SW\downarrow$ [$Wm^{-2}$]'
p = SWmonth.plot(col='month',cmap=plt.cm.gist_heat, col_wrap=4)

i=1
for ax in p.axes.flat:
    mont=calendar.month_abbr[i]
    ax.set_title('$SW\downarrow$  '+mont)
    ax.axis('off')
    i=i+1
plt.savefig('MonthlyHySNSW.png')
plt.show()

LWSN=xr.open_dataset(savedir+'HySN_Surface_Downwelling_Longwave_Radiation_'+str(yr)+'.nc', chunks={'time': 30})
LWmonth=LWSN['rlds'].groupby('time.month').mean('time')
LWmonth.name='$LW\downarrow$  [$Wm^{-2}$]'
p = LWmonth.plot(col='month',cmap=plt.cm.gist_heat, col_wrap=4)

i=1
for ax in p.axes.flat:
    mont=calendar.month_abbr[i]
    ax.set_title('$LW\downarrow$  '+mont)
    ax.axis('off')
    i=i+1
plt.savefig('MonthlyHySNLW.png')
plt.show()
