from lsreduce import io, astrometry, reduction
import datetime
import ephem
import numpy as np

siteinfo = io.read_siteinfo('/home/talens-irex/Research/MASCARA/conf/LS/siteinfo.dat', 'LaSilla')
wcspars, polpars, astromask = io.read_astromaster('/home/talens-irex/Research/MASCARA/conf/LS/astromasterLSEnight.hdf5')
astro = astrometry.Astrometry(wcspars, polpars, astromask)

twilight=5
t = datetime.datetime.utcnow() + datetime.timedelta(minutes=15)  # Add some grace for safety.
print t
t = astrometry.last_sunset(siteinfo, t=t, horizon=twilight)
print t
t = astrometry.last_sunset(siteinfo, t=t, horizon=-10)
print t

station = np.recarray(50, dtype=[('jd', 'float'), ('lst', 'float')])
station['jd'] = 2459323.18926 + np.arange(50)*6.4/(60*60*24)
station['lst'] = np.arange(50)*6.4/(60*60)
print station

station = reduction.sunmoon2station(station, siteinfo, astro)

print station

print astrometry.sun_position(siteinfo)