import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt
import argparse
import sys

def periastron_advance_constant(p_orb,ecc,omdot,domdot):

	p_orb=p_orb*24*3600
	cnt=3*(6.67430e-11/299792458**3)**(2/3)*(p_orb/(2*np.pi))**(-5/3)/(1-ecc**2)
	omdot=omdot*(np.pi/180)*(1/(31557600))
	domdot=domdot*(np.pi/180)*(1/(31557600))

	mtot=((omdot/cnt)**(3/2))
	dmtot=(domdot*3/2)*((omdot**(1/2))/(cnt**(3/2)))

	mtot=mtot/1.9891e30    #In solar massaes
	dmtot=dmtot/1.9891e30

	return mtot,dmtot

def einstein_delay_constant(p_orb,ecc,gamma,dgamma):

	p_orb=p_orb*24*3600
	cnt=ecc*(6.67430e-11/299792458**3)**(2/3)*(p_orb/(2*np.pi))**(1/3)

	Mgamma=gamma/cnt  #In kg
	dMgamma=dgamma/cnt

	return Mgamma/(1.9891e30**(2/3)),dMgamma/(1.9891e30**(2/3)) #In solar masses

def orbital_decay_constant(p_orb,ecc,pbdot,dpbdot):

	p_orb=p_orb*24*3600
	cnt=(192*np.pi/5)*(6.67430e-11/299792458**3)**(5/3)*(p_orb/(2*np.pi))**(-5/3)*(1+(73/23)*ecc**2+(37/96)*ecc**4)/(1-ecc**2)**(7/2)

	Mchirp=pbdot/cnt  #In kg
	dMchirp=dpbdot/cnt

	return Mchirp/(1.9891e30**(5/3)),dMchirp/(1.9891e30**(5/3)) #In solar masses

def constraint_from_h3_stig(p_orb,x,stig,dstig,h3,dh3):

	s=2*stig/(1+stig**2)
	ds=2*dstig*(1/(1+stig**2)-2*(stig**2)/((1+stig**2)**2))

	mcomp=(299792458**3)*(h3/(stig**3))/6.67430e-11
	dmcomp=np.sqrt(np.square((299792458**3)*(dh3/(stig**3))/6.67430e-11)+np.square(dstig*(299792458**3)*(h3/(stig**4))/6.67430e-11))

	# In solar masses.
	return mcomp/1.9891e30,dmcomp/1.9891e30,s,ds

def read_ephemeris(file):
	ephemeris=open(file,"r")
	ecc=0
	omega=0
	for line in ephemeris:
		line=line.replace(" ","")
		line=line.replace("	","")
		if line[:2]=="F0":
			f0=float(line[2:])
		elif line[:2]=="PB":
			p_orb=float(line[2:])
		elif line[:2]=="A1":
			x=float(line[2:])
		elif line[:3]=="ECC":
			ecc=float(line[3:])
		elif line[:2]=="OM":
			omega=float(line[2:])
		elif line[:2]=="T0":
			periastron=float(line[2:])
	if f0 and p_orb and x and periastron:
		print("Pulsar parameters loaded from {}:".format(file))
		print("- Pulsar frequency: {} Hz".format(f0))
		print("- Orbital period: {} days".format(p_orb))
		print("- Projected axis: {} ls".format(x))
		print("- Eccentricity: {}".format(ecc))
		print("- Periastron angle: {} degrees".format(omega))
		print("- Periastron passage: {} (MJD)".format(periastron))
		print("")
	else:
		sys.exit("The ephemeris file doesn't include all the necessary parameters.")
	return f0,p_orb,x,ecc,omega,periastron

parser=argparse.ArgumentParser(description="Take in orbital parameters and post-Keplerian effects and constrain the mass and inclination angle.")
parser.add_argument("--ephemeris",help="Fitorbit-format ephemeris. If given all other inputs will be ignored.")
parser.add_argument("-p","--period",type=float,help="Orbital period in days")
parser.add_argument("-x","--axis",type=float,help="Projected semimajor axis in ls.")
parser.add_argument("-e","--eccentricity",type=float,help="Excentricity. If not given, assumed 0. Needed for higher order estimations.")
parser.add_argument("--omdot",help="Periastron advance (value+/-uncertainty) in º/yr")
parser.add_argument("--gamma",help="Einstein delay (value+/-uncertainty) in s.")
parser.add_argument("--pbdot",help="Spin-down (value+/-uncertainty) in s/s.")
parser.add_argument("--h3",help="Orthometric amplitude amplitude of saphiro delay (value+/-uncertainty) in s.")
parser.add_argument("--stig",help="Orthometric amplitude amplitude of saphiro delay (value+/-uncertainty).")
parser.add_argument("--demodulate",type=float,help="If set at any value, it will print out all mass ranges from 90 to 30º.")
parser.add_argument("-i","--inclination",type=float,help="Force custom inclination angle into the mass function (degrees).")
parser.add_argument("-v","--verbose",action="store_true")
args = parser.parse_args()

print(" ")

if args.ephemeris:
	(f0,p_orb,x,ecc,omega,periastron)=read_ephemeris(args.ephemeris)
else:
	if args.period:
		p_orb=args.period
		if args.verbose==True:
			print("Orbital period: {} days".format(p_orb))
			print(" ")
	else:
		sys.exit("Please specify orbital period in days with -p")
	if args.axis:
		x=args.axis
		if args.verbose==True:
			print("Projected axis: {} ls".format(x))
			print(" ")
	else:
		sys.exit("Please specify projected semimajor axis in days with -x")
	if args.eccentricity:
		ecc=args.eccentricity
		if args.verbose==True:
			print("Eccentricity: {} ls".format(ecc))
			print(" ")
	else:
		ecc=0
		if args.verbose==True:
			print("Eccentricity assumed to be 0")
			print(" ")

if args.omdot:
	omdot=float(args.omdot.split("+/-")[0])
	domdot=float(args.omdot.split("+/-")[1])
	if args.verbose==True:
		print("Constraining from omdot {} º/yr.".format(args.omdot))
	(mtot,dmtot)=periastron_advance_constant(p_orb,ecc,omdot,domdot)
	if args.verbose==True:
		print("Total mass: {}+/-{} solar masses.".format(mtot,dmtot))
		print("")

if args.gamma:
	gamma=float(args.gamma.split("+/-")[0])
	dgamma=float(args.gamma.split("+/-")[1])
	if args.verbose==True:
		print("Constraining from gamma {} s.".format(args.gamma))
	(mgamma,dmgamma)=einstein_delay_constant(p_orb,ecc,gamma,dgamma)
	if args.verbose==True:
		print("Einstein delay mass function: {}+/-{} (solar masses)^(2/3).".format(mgamma,dmgamma))
		print("")

if args.pbdot:
	pbdot=float(args.pbdot.split("+/-")[0])
	dpbdot=float(args.pbdot.split("+/-")[1])
	if args.verbose==True:
		print("Constraining from pbdot.")
	(mchirp,dmchirp)=orbital_decay_constant(p_orb,ecc,pbdot,dpbdot)
	if args.verbose==True:
		print("Orbital decay mass function: {}+/-{} (solar masses)^(5/3)".format(mchirp,dmchirp))
		print("Chirp mass: {}+/-{} solar masses".format(mchirp**(3/5),(3/5)*dmchirp*mchirp**(-2/5)))
		print("")

if args.h3 and args.stig:
	h3=float(args.h3.split("+/-")[0])
	dh3=float(args.h3.split("+/-")[1])
	stig=float(args.stig.split("+/-")[0])
	dstig=float(args.stig.split("+/-")[1])
	if args.verbose==True:
		print("Constraining from Shapiro delay.")
	(mcomp,dmcomp,s,ds)=constraint_from_h3_stig(p_orb,x,stig,dstig,h3,dh3)
	if args.verbose==True:
		print("Companion mass: {}+/-{}".format(mcomp,dmcomp))
		print("Inclination angle: {}+/-{}".format(np.arcsin(s)*(180/np.pi),(np.arcsin(s+ds)-np.arcsin(s-ds))*(180/np.pi)/2))
		print("")

# Start computing constraints like a madman:

if args.omdot:

	print("Masses from periastron advance an mass function:")
	print(" ")
	print("Mtot= {} +/- {} ({} - {}) solar masses".format(mtot,dmtot,mtot-dmtot,mtot+dmtot))
	print(" ")

	if args.demodulate:

		print("i(º), sin(i), Mcomp(Msun), Mpulsar(Msun)")

		file=open("inclinations.txt","w")

		for s in np.arange(1,0.51,-args.demodulate):

			mcomp=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)*(mtot**2)/1.9891e30)**(1/3)/s
			mpulsar=mtot-mcomp

			print("{}, {}, {}, {}".format(180*np.arcsin(s)/np.pi,s,mcomp,mpulsar))

			file.write("{},{}\n".format(s,mcomp))
	
		file.close()
		print(" ")

	elif args.inclination:

		s=np.sin(args.inclination*np.pi/180)
		mcomp=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)*(mtot**2)/1.9891e30)**(1/3)/s
		dmcomp=(2*dmtot*mtot**(-1/3)/3)*(299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)/1.9891e30)**(1/3))/s
		mpulsar=mtot-mcomp
		dmpulsar=np.sqrt(dmtot**2+dmcomp**2)

		print("At i= {}º".format(args.inclination))
		print("Mcomp= {} +/ {} ({} - {}) solar masses".format(mcomp,dmcomp,mcomp-dmcomp,mcomp+dmcomp))
		print("Mpulsar= {} +/- {} ({} - {}) solar masses".format(mpulsar,dmpulsar,mpulsar-dmpulsar,mpulsar+dmpulsar))
		print(" ")

	else:

		mcomp=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)*(mtot**2)/1.9891e30)**(1/3)
		dmcomp=(2*dmtot*mtot**(-1/3)/3)*(299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)/1.9891e30)**(1/3))
		mpulsar=mtot-mcomp
		dmpulsar=np.sqrt(dmtot**2+dmcomp**2)

		print("At i= 90º")
		print("Mcomp= {} +/ {} ({} - {}) solar masses".format(mcomp,dmcomp,mcomp-dmcomp,mcomp+dmcomp))
		print("Mpulsar= {} +/- {} ({} - {}) solar masses".format(mpulsar,dmpulsar,mpulsar-dmpulsar,mpulsar+dmpulsar))
		print(" ")

		mcomp=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)*(mtot**2)/1.9891e30)**(1/3)/np.sin(60*np.pi/180)
		dmcomp=(2*dmtot*mtot**(-1/3)/3)*(299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)/1.9891e30)**(1/3))/np.sin(60*np.pi/180)
		mpulsar=mtot-mcomp
		dmpulsar=np.sqrt(dmtot**2+dmcomp**2)

		print("At i= 60º")
		print("Mcomp= {} +/ {} ({} - {}) solar masses".format(mcomp,dmcomp,mcomp-dmcomp,mcomp+dmcomp))
		print("Mpulsar= {} +/- {} ({} - {}) solar masses".format(mpulsar,dmpulsar,mpulsar-dmpulsar,mpulsar+dmpulsar))
		print(" ")

		mcomp=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)*(mtot**2)/1.9891e30)**(1/3)/np.sin(45*np.pi/180)
		dmcomp=(2*dmtot*mtot**(-1/3)/3)*(299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)/1.9891e30)**(1/3))/np.sin(45*np.pi/180)
		mpulsar=mtot-mcomp
		dmpulsar=np.sqrt(dmtot**2+dmcomp**2)

		print("At i= 45º")
		print("Mcomp= {} +/ {} ({} - {}) solar masses".format(mcomp,dmcomp,mcomp-dmcomp,mcomp+dmcomp))
		print("Mpulsar= {} +/- {} ({} - {}) solar masses".format(mpulsar,dmpulsar,mpulsar-dmpulsar,mpulsar+dmpulsar))
		print(" ")
		print(" ")

if args.omdot and args.gamma:
	
	mcomp=(-mtot+np.sqrt(mtot**2+4*mtot**(4/3)*mgamma))/2
	mpulsar=mtot-mcomp
	s=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)*((mtot**2)/(1.9891e30*mcomp**3)))**(1/3)

	print("Masses from periastron advance and Einstein delay:")
	print(" ")
	print("Mcomp= {} solar masses".format(mcomp))
	print("Mtot= {} solar masses".format(mtot))
	print("Mpulsar= {} solar masses".format(mpulsar))
	print("Inclination angle from mass function: s={} (i={} º)".format(s,np.arcsin(s)*(180/np.pi)))
	print(" ")
	print(" ")

if args.omdot and args.pbdot:

	# Neither pdot or pbdot are able to give a proper constraint.
	
	m1=(mtot-np.sqrt(mtot**2-4*mtot**(1/3)*mchirp))/2
	m2=(mtot+np.sqrt(mtot**2-4*mtot**(1/3)*mchirp))/2
	s1=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)*((mtot**2)/(1.9891e30*m1**3)))**(1/3)
	s2=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67430e-11)*((mtot**2)/(1.9891e30*m2**3)))**(1/3)

	print("Masses from periastron advance and orbital decay:")
	print(" ")
	print("M1= {} solar masses".format(m1))
	print("M2= {} solar masses".format(m2))
	print("Mtot= {} solar masses".format(mtot))
	print("Inclination angle if 1 is companion: s={} (i={} º)".format(s1,np.arcsin(s1)*(180/np.pi)))
	print("Inclination angle if 2 is companion: s={} (i={} º)".format(s2,np.arcsin(s2)*(180/np.pi)))
	print(" ")
	print(" ")
