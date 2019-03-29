import numpy as np
import math
import datetime
import time
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.jplhorizons import Horizons

np.set_printoptions(suppress=True)
# function reads catalogue and returns a list where each element is a plate
# replaces spaces with zeros to avoid errors
def getData():
	cat = open('platelist.txt', 'r')
	cat = cat.readlines()
	for i in range(0,len(cat)):
		cat[i] = cat[i].replace(" ","0")
	return cat

def getDataprint():
	cat = open('catprint.txt', 'r')
	cat = cat.readlines()
	for i in range(0,len(cat)):
		cat[i] = cat[i].replace(" ","0")
	return cat

# converts the catalogue from string list to data array
def mkArray(cat, lng):
	data = np.zeros((len(cat), 8))
	for i in range(0, len(cat)):
		data[i][0] = cat[i][2:7]
		data[i][1] = getJD(cat[i][30:40],lng)
		data[i][2] , data[i][3] = RaDecRad(cat[i][20:30], cat[i][2:7])
		data[i][4] = cat[i][52:56]
	invalid = []
	for i in range(0,len(cat)):
		if cat[i][7]!='0':
			invalid.append(i)
	for i in sorted(invalid, reverse=True):
		data = np.delete(data,i,0)
	return data

# determines which plates are to be discarded
def magLims(plist, cat):
	for i in range(0,len(cat)):
		if plist[i][1] == 'U':
			maglim = 100
		elif plist[i][1] == 'B':
			maglim = 21+2.5*np.log10(cat[i][4]/600.)
		elif plist[i][1] == 'J':
			maglim = 22.5+2.5*np.log10(cat[i][4]/600.)
		elif plist[i][1] == 'V':
			maglim = 21+2.5*np.log10(cat[i][4]/600.)
		elif plist[i][1] == 'R':
			maglim = 22+2.5*np.log10(cat[i][4]/900.)
		elif plist[i][1] == 'I':
			maglim = 19.5+2.5*np.log10(cat[i][4]/900.)
		else:
			maglim = 100
		cat[i][7] = maglim
	return cat


# calculates JDate (currently only sensible years)
def getJDate(y,m,d):
	if m <=2:
		y_ = y-1
		m_ = m+12
	else:
		y_ = y
		m_ = m
	A = np.trunc(y_/100.)
	B = 2 - A + np.trunc(A/4.)
	C = math.trunc(365.25*y_)
	D = math.trunc(30.6001*(m_+1))
	JDate = B+C+D+d+1720994.5
	return JDate

#converts JD at 0h and GST to UT
def GSTUT(JD, GST):
	S =JD-2451545.0
	T = S/36525.0
	T0 = 6.697374558+(2400.051336*T)+(0.000025862*T**2)
	if T0<0:
		while T0<0: T0=T0+24
	elif T0>=24:
		while T0>24: T0=T0-24
	else: T0=T0
	A = GST-T0
	if A<0:
		while A<0: A=A+24
	elif A>=24:
		while A>=24: A=A-24
	else: A=A
	UT = A*0.9972675663
	return UT

# takes date and time string and longitude in degrees and gives JD
def getJD(Raw, lng):
	if int(Raw[0:2])>=50: year = 1900+float(Raw[0:2])
	else: year = 2000+float(Raw[0:2])
	month = float(Raw[2:4])
	day = float(Raw[4:6])
	hour = Raw[6:8]
	hour = float(hour.replace(" ","0"))
	minute = float(Raw[8:10])
	# converting LST to fractional hours
	LST = hour+(minute/60.)
	lng = lng/15.
	GST = LST-lng
	#GST in hours
	if GST > 24: GST=GST-24
	elif GST <= 0: GST=GST+24
	else: GST=GST
	#finding the Julian Date
	JDate = getJDate(year,month,day)
	#if LST<lng: JDate=JDate-1
	#else: JDate=JDate
	UT = GSTUT(JDate, GST)
	JD = JDate + UT/24.
	return JD

# takes radec as it appears in the cat and returns them seperately and in radians
def RaDecRad(radec, count):
	dmin = int(radec[8:10])/60.
	dec = int(radec[6:8])+dmin
	dec = math.pi*dec/180.
	if radec[5]=='+':
		Dec=dec
	else: Dec = -dec
	ramin = int(radec[2:5])/600.
	ra = int(radec[0:2])+ ramin
	ra = ra*15
	gc = SkyCoord(ra=ra*u.degree, dec=Dec*u.radian, frame='fk4')
	gc = gc.transform_to('icrs')
	Ra = gc.ra.radian
	Dec = gc.dec.radian
	print count
	return Ra, Dec

# converts from spherical to tangent plane coordinates wrt a reference point z
# returns a value for xi eta in radians
"""
Adapted from Starlink SLALIB Positional Astronomy Library
"""
def ds2tp(raz,decz,ra,dec):
	radif = ra-raz
	denom = math.sin(dec)*math.sin(decz)+math.cos(dec)*math.cos(decz)*math.cos(radif)
	xi = math.cos(dec)*math.sin(radif)/denom
	eta = ((math.sin(dec)*math.cos(decz)-math.cos(dec)*math.sin(decz)*math.cos(radif))/denom)
	TINY = 10**-6
	if denom > TINY:
		j = 0
	elif denom >= 0:
		j = 1
	elif denom > -TINY:
		j = 2
	else:
		j = 3
	return xi , eta, j

# checks if an object will appear on a plate
def match(xi,eta,thresh, j):
	if abs(xi)<=thresh and abs(eta)<=thresh and j==0:
		return 1
	else:
		return 0

# returns x-y position of object on plate in actual size
def onPlate(hitpos, psize, pscale, cat):
	for i in range(0,len(hitpos)):
		xi = np.degrees(hitpos[i][2])*60*60
		x = (xi/pscale)+psize/2.
		eta = np.degrees(hitpos[i][3])*60*60
		y = (eta/pscale)+psize/2.
		hitpos[i][2] = x
		hitpos[i][3] = y
	return hitpos

# changes the JD to be half way through the exposure time
def midExposure(cat):
	for i in range(0, len(cat)):
		print 'length in min ' + str(cat[i][4]/10.)
		extime = (cat[i][4]/10.)*60
		print 'length in seconds ' +str(extime)
		extimedays = extime/86400.
		print extimedays
		cat[i][1] += (extimedays/2.)
	return cat

# submits horizons queries and returns radec of object for every epoch
# places the values inside the data array
def horizons(location, target, cat, id_type):
	extra = np.zeros((len(cat)))
	extra2 = np.zeros((len(cat)))
	extra3 = np.zeros((len(cat)))
	cat = np.column_stack((cat, extra, extra2, extra3))
	Rarray = np.empty(1)
	Decarray = np.empty(1)
	galbarray = np.empty(1)
	gallarray = np.empty(1)
	magarray = np.empty(1)
	for i in range(0,len(cat),350):
		progress = round((i*100/float(len(cat))), 1)
		print str(progress) + '%'
		epoch = cat[i:i+350]
		epochs = []
		for j in range(0,len(epoch)):
			epochs.append(epoch[j][1])
		epochs = np.array(epochs)
		obj = Horizons(id=target, location=location, id_type=id_type, epochs=epochs)
		eph = obj.ephemerides()
		if i==0:
			targetname = eph['targetname'][i]
		newRA = np.array(eph['RA'])
		newDec = np.array(eph['DEC'])
		newgalb = np.array(eph['GlxLat'])
		newgall = np.array(eph['GlxLon'])
		newMag = np.array(eph['V'])
		Rarray = np.concatenate((Rarray, newRA))
		Decarray = np.concatenate((Decarray, newDec))
		galbarray = np.concatenate((galbarray, newgalb))
		gallarray = np.concatenate((gallarray, newgall))
		magarray = np.concatenate((magarray, newMag))
	for i in range(0, len(cat)):
		cat[i][5] = Rarray[i+1]
		cat[i][6] = Decarray[i+1]
		cat[i][8] = gallarray[i+1]
		cat[i][9] = galbarray[i+1]
		cat[i][10] = magarray[i+1]
	cat[:,5] = np.radians(cat[:,5])
	cat[:,6] = np.radians(cat[:,6])
	targetname = targetname.replace(' ', '_')
	return cat , targetname

# gives a list of plate numbers with potential hits on them
def hitList(cat, thresh):
	hitnumlist = []
	hitnum = 0
	hitList = []
	hitposx = []
	hitposy = []
	hitra = []
	hitdec = []
	hitgall = []
	hitgalb = []
	hitmag = []
	hithormag = []
	for i in range(0, len(cat)):
		raz = cat[i][2]
		decz = cat[i][3]
		ra = cat[i][5]
		dec = cat[i][6]
		xi , eta , j = ds2tp(raz, decz, ra, dec)
		hit = match(xi, eta, thresh, j)
		if hit == 1:
			hitnumlist.append(hitnum)
			hitnum += 1
			hitList.append(cat[i][0])
			hitposx.append(xi)
			hitposy.append(eta)
			hitra.append(cat[i][5])
			hitdec.append(cat[i][6])
			hitgall.append(cat[i][7])
			hitgalb.append(cat[i][8])
			# limit - obs
			magabove = cat[i][7]-cat[i][10]
			if cat[i][7] == 100:
				magabove = 100+cat[i][10]
			hitmag.append(magabove)
			hithormag.append(cat[i][10])
	hitpos = np.column_stack((hitnumlist,hitList,hitposx,hitposy,hitra,hitdec,hitgall,hitgalb,hitmag,hithormag))
	return hitpos

def offPlane(hitpos):
	for i in xrange(len(hitpos)-1, -1, -1):
		if abs(hitpos[i][7]) < 0.1:
			hitpos = np.delete(hitpos, i, axis=0)
			print 'Deleted'
		else:
			print hitpos[i][7]
	for i in range(0, len(hitpos)):
		hitpos[i][0] = i
	return hitpos

def catradec(hitpos,cat):
	extra = np.zeros((1,len(hitpos)))
	hitpos = np.column_stack((hitpos,extra))
	extra2 = np.copy(extra)
	extra3 = np.copy(extra)
	hitpos = np.column_stack((hitpos,extra2))
	for i in range(0,len(hitpos)):
		hitpos[i][3] = cat[hitpos[i][0]][5] 
		hitpos[i][4] = cat[hitpos[i][0]][6]
		hitpos[i][5] = i
	return hitpos

def export(hitpos, targetname):
	header , filename = getHeader(targetname)
	fmt = '%.0f %.0f %.2f %.2f %.9f %.9f %.2f %.2f %.3f %.3f'
	np.savetxt(filename, hitpos, fmt=fmt, header=header)
	return

def getHeader(targetname):
	header = str(targetname)+'\n'
	filename = './objects/'+targetname+'.txt'
	now = datetime.datetime.now()
	header += str(now)+'\n'
	header += 'index p_num xpos(mm) ypos(mm) RA(rad) Dec(rad) gLon(deg) gLat(deg) mag\n'
	return header , filename

#method to write cat to a file
def catout(cat):
	fname = 'cat.txt'
	header = 'cat.txt\n'
	header += 'catalogue in array form\n'
	header += 'plate epoch plateRA plateDec extime horizonsRA horizonsDec'
	np.savetxt(fname, cat, header=header)

def catin():
	cat = np.genfromtxt('cat.txt')
	return cat

def catinprint():
	cat = np.genfromtxt('catprint.txt')
	return cat

# method for setting up the run
def init():
	lng = 149.0661
	thresh = np.radians(3.31)
	id_type = int(input('Body type? smallbody = 1 majorbody = 2 : '))
	if id_type ==1:
		id_type= 'smallbody'
	else: id_type= 'majorbody'
	body = str(input('Body Number : '))
	pscale = 67.12
	psize = 355.
	return lng , thresh , id_type , body , pscale , psize 

"""
# main for extracting a catalogue
def main():
	begin = time.time()
	lng, thresh, id_type, body, pscale, psize = init()
	plist = getData()
	cat = mkArray(plist, lng)
	cat = magLims(plist, cat)
	catout(cat)
	print "success"
	end = time.time()
	duration = end-begin
	print duration
"""

# main for running with a working catalogue
def main():
	begin = time.time()
	lng, thresh, id_type, body, pscale, psize = init()
	cat = catin()
	cat , targetname = horizons(260,body,cat,id_type)
	print cat[500]
	hitpos = hitList(cat, thresh)
	print hitpos
	hitposDown = onPlate(hitpos, psize, pscale, cat)
	hitposDown = offPlane(hitposDown)
	print 'emulsion down'
	print hitposDown
	export(hitposDown, targetname)
	print 'success'
	end = time.time()
	duration = end-begin
	print duration

main()