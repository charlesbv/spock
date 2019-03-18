link_in = "https://iswa.ccmc.gsfc.nasa.gov/iswa_data_tree/model/heliosphere/wsa-enlil-cone/velocity-density-timeline-DATA/"
#link_in = "https://iswa.ccmc.gsfc.nasa.gov/iswa_data_tree/model/heliosphere/wsa-enlil-cone/velocity-density-timeline-DATA/2018/07/"
from bs4 import BeautifulSoup
import ipdb
import requests
import os
if link_in[-1] != '/':
    link_in = link_in + '/'
r  = requests.get(link_in)
data = r.text
soup = BeautifulSoup(data)
for link in soup.find_all('a'):
    print 'YYYYYYYYYYYYYYY'
    yy = link.get('href')
    if (( yy[0] != '/') & (yy[0] != '?')): # yy could be wierd stuff - we only want if it's a year
        print yy
        linkyear = link_in + yy
        ryear  = requests.get(linkyear)
        datayear = ryear.text
        soupyear = BeautifulSoup(datayear)
        for linkmm in soupyear.find_all('a'):
            mm = linkmm.get('href')
            if (( mm[0] != '/') & (mm[0] != '?')): # mm could be wierd stuff - we only want if it's a month
                print mm
                linkmonth = linkyear + mm
                rmonth  = requests.get(linkmonth)
                datamonth = rmonth.text
                soupmonth = BeautifulSoup(datamonth)
                for linkdd in soupmonth.find_all('a'):
                    dd = linkdd.get('href')
                    if (( dd[0] != '/') & (dd[0] != '?')): # dd could be wierd stuff - we only want if it's a day
                        linkday = linkmonth + dd
                        os.system("wget " + linkday)

