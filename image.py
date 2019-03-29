# program to get an image of an area
# importing required packages
import numpy as np
from selenium import webdriver
import time
import os
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord


catname = str(sys.argv[1])
hitnum = int(sys.argv[2])
cat = np.genfromtxt(catname)
#cat = np.genfromtxt('Mars_\(499\).txt')
print cat

ra = cat[hitnum][4]
dec = cat[hitnum][5]
ra = ra*u.radian
dec = dec*u.radian
print ra
print dec
coords = SkyCoord(ra=ra, dec=dec)
print coords.ra.hms.h
print coords.dec.dms

coordstring = str(coords.ra.hms.h)+' '
if coords.ra.hms.h == 0:
	coordstring += str(coords.ra.hms.m)+' '
else:
	coordstring += str(abs(coords.ra.hms.m))+' '
coordstring += str(np.round(abs(coords.ra.hms.s), decimals=2))+' '
coordstring += str(coords.dec.dms.d)+' '
if coords.dec.dms.d ==0:
	coordstring += str(coords.dec.dms.m)+' '
else:
	coordstring += str(abs(coords.dec.dms.m))+' '
coordstring += str(np.round(abs(coords.dec.dms.s), decimals=2))
print coordstring
browser = webdriver.Chrome()
browser.get('http://www-wfau.roe.ac.uk/sss/pixel.html')
coordElem = browser.find_element_by_name("coords")
coordElem.send_keys(coordstring)
sizeElem = browser.find_element_by_name("size")
sizeElem.send_keys('15')
sizeElem.submit()
time.sleep(3)
dload = browser.find_element_by_link_text('here')
dload.click()
time.sleep(2)
browser.quit()
os.system('./get.sh')