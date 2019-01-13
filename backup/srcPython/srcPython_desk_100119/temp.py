isc = 6
nb_tle_for_sc = len(ecc[isc])
nb_sec_since_date_ref = np.zeros([nb_tle_for_sc])
for itle in range(nb_tle_for_sc):
    nb_sec_since_date_ref[itle] = ( datetime.strptime( date[isc][itle].split('.')[0], "%y%j" ) - date_ref ).total_seconds() + np.float( '0.' + date[isc][itle].split('.')[1] ) * 24 * 3600
save_last_tle_date[isc] = nb_sec_since_date_ref[-1]
x = nb_sec_since_date_ref 
y = ( np.array(sma[isc]) ) / ( np.array(sma[isc][0]) )

nslope = len(x) - 1
yslope = np.zeros([nslope])
xslope = x[:-1]
for itle in range(nslope):
    yslope[itle] = ( y[itle+1] - y[itle] ) / ( x[itle+1] - x[itle] )
