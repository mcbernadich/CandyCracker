import numpy as np
import matplotlib.pyplot as plt
#from numba import jit
import argparse
import sys

def normal_distribution(y,nominal_y,dy):
	return np.exp(-((y-nominal_y)**2)/(2*(dy**2)))

def periastron_advance_constraint(p_orb,ecc,x,omdot,domdot):

	p_orb=p_orb*24*3600
	cnt=(1.9891e30**(2/3))*(180/np.pi)*(24*3600*365)*3*(6.67408e-11/299792458**3)**(2/3)*(p_orb/(2*np.pi))**(-5/3)/(1-ecc**2)
	nominal_omdot=omdot
	domdot=domdot
	one_sigma=normal_distribution(nominal_omdot+domdot,nominal_omdot,domdot)

	#Initialize the m2-m1 axes.
	m1_ax=np.arange(1.0,2.4,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(m1,m2)=np.meshgrid(m1_ax,m2_ax)

	#Evaluate omdot at each point of the m1-m2 grid and compute the distribution.
	omdot=cnt*((m1+m2)**(2/3))
	m1_m2_distribution=normal_distribution(omdot,nominal_omdot,domdot)

	#Initialize the m2-cosi axes.
	cosi_ax=np.arange(0.005,1,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(cosi,m2)=np.meshgrid(cosi_ax,m2_ax)

	#Recompute the m1 grid out of m2 and cosi, but with the same dimensionality.
	m1=(m2**3*(6.67408e-11*1.9891e30/(4*np.pi**2))*((np.sqrt(1-cosi**2)**3)*(p_orb**2)/(x*299792458)**3))**(1/2)-m2

	#Evaluate omdot at each point of the cosi-m2 grid and compute the distribution.
	omdot=cnt*((m1+m2)**(2/3))
	cosi_m2_distribution=normal_distribution(omdot,nominal_omdot,domdot)

#	sini=299792458*x*((1/p_orb**2)*(4*np.pi**2/6.67408e-11)*(((m1+m2)**2)/(1.9891e30*m2**3)))**(1/3)
#	sini[sini>1]=1.
#	cosi=np.sqrt(1-sini**2)

	return m1_m2_distribution,cosi_m2_distribution,one_sigma

def constraint_from_h3(p_orb,x,h3,dh3):

	p_orb=p_orb*24*3600
	nominal_h3=h3
	dh3=dh3
	one_sigma=normal_distribution(nominal_h3+dh3,nominal_h3,dh3)

	#Initialize the cosi-m2 axes.
	cosi_ax=np.arange(0.005,1,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(cosi,m2)=np.meshgrid(cosi_ax,m2_ax)

	#Evaluate h3 at each point of the m2-cosi grid and compute the distribution.
	h3=(6.67408e-11*1.9891e30*m2/(299792458**3))*(((1-cosi)/(1+cosi))**(3/2))
	h3=h3.astype("float")
	cosi_m2_distribution=normal_distribution(h3,nominal_h3,dh3)

	#Initialize the m2-m1 axes.
	m1_ax=np.arange(1.0,2.4,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(m1,m2)=np.meshgrid(m1_ax,m2_ax)

	#Recompute the cosi grid out of m1 and m2
	sini=299792458*x*((1/p_orb**2)*(4*np.pi**2/6.67408e-11)*(((m1+m2)**2)/(1.9891e30*m2**3)))**(1/3)
	sini[sini>1]=1.
	cosi=np.sqrt(1-sini**2)

	#Evaluate omdot at each point of the m1-m2 grid and compute the distribution.
	h3=(6.67408e-11*1.9891e30*m2/(299792458**3))*(((1-cosi)/(1+cosi))**(3/2))
	h3=h3.astype("float")
	m1_m2_distribution=normal_distribution(h3,nominal_h3,dh3)

	# In solar masses.
	return m1_m2_distribution,cosi_m2_distribution,one_sigma

def constraint_from_stig(p_orb,x,stig,dstig):

	p_orb=p_orb*24*3600
	nominal_stig=stig
	dstig=dstig
	one_sigma=normal_distribution(nominal_stig+dstig,nominal_stig,dstig)

	#Initialize the cosi-m2 axes.
	cosi_ax=np.arange(0.005,1,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(cosi,m2)=np.meshgrid(cosi_ax,m2_ax)

	#Evaluate stig at each point of the m2-cosi grid and compute the distribution.
	stig=((1-cosi)/(1+cosi))**(1/2)
	cosi_m2_distribution=normal_distribution(stig,nominal_stig,dstig)

	#Initialize the m2-m1 axes.
	m1_ax=np.arange(1.0,2.4,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(m1,m2)=np.meshgrid(m1_ax,m2_ax)

	#Recompute the cosi grid out of m1 and m2
	sini=299792458*x*((1/p_orb**2)*(4*np.pi**2/6.67408e-11)*(((m1+m2)**2)/(1.9891e30*m2**3)))**(1/3)
	sini[sini>1]=1.
	cosi=np.sqrt(1-sini**2)

	#Evaluate stig at each point of the m1-m2 grid and compute the distribution.
	stig=((1-cosi)/(1+cosi))**(1/2)
	m1_m2_distribution=normal_distribution(stig,nominal_stig,dstig)

	# In solar masses.
	return m1_m2_distribution,cosi_m2_distribution,one_sigma

def constraint_from_gamma(p_orb,ecc,x,gamma,dgamma):

	p_orb=p_orb*24*3600
	cnt=ecc*(1.9891e30*6.67408e-11/299792458**3)**(2/3)*(p_orb/(2*np.pi))**(1/3)
	nominal_gamma=gamma
	dgamma=dgamma
	one_sigma=normal_distribution(nominal_gamma+dgamma,nominal_gamma,dgamma)

	#Initialize the m2-m1 axes.
	m1_ax=np.arange(1.0,2.4,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(m1,m2)=np.meshgrid(m1_ax,m2_ax)

	#Evaluate gamma at each point of the m1-m2 grid and compute the distribution.
	gamma=cnt*m2*(m1+2*m2)/((m1+m2)**(4/3))
	m1_m2_distribution=normal_distribution(gamma,nominal_gamma,dgamma)

	#Initialize the m2-cosi axes.
	cosi_ax=np.arange(0.005,1,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(cosi,m2)=np.meshgrid(cosi_ax,m2_ax)

	#Recompute the m1 grid out of m2 and cosi, but with the same dimensionality.
	m1=(m2**3*(6.67408e-11*1.9891e30/(4*np.pi**2))*((np.sqrt(1-cosi**2)**3)*(p_orb**2)/(x*299792458)**3))**(1/2)-m2

	#Evaluate gamma at each point of the cosi-m2 grid and compute the distribution.
	gamma=cnt*m2*(m1+2*m2)/((m1+m2)**(4/3))
	cosi_m2_distribution=normal_distribution(gamma,nominal_gamma,dgamma)

	# In solar masses.
	return m1_m2_distribution,cosi_m2_distribution,one_sigma

def constraint_from_DDGR(p_orb,x,m2,mT,dm2,dmT):

	p_orb=p_orb*24*3600
	nominal_m2=m2
	nominal_mT=mT
	one_sigma=normal_distribution(nominal_m2+dm2,nominal_m2,dm2)*normal_distribution(nominal_mT+dmT,nominal_mT,dmT)
	print(one_sigma)

	#Initialize the m2-m1 axes.
	m1_ax=np.arange(1.0,2.4,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(m1,m2)=np.meshgrid(m1_ax,m2_ax)

	#Evaluate m2 and mtot at each point of the m1-m2 grid and compute the distribution
	M2=m2
	MTOT=m1+m2
	m1_m2_distribution=normal_distribution(M2,nominal_m2,dm2)*normal_distribution(MTOT,nominal_mT,dmT)

	#Initialize the m2-cosi axes.
	cosi_ax=np.arange(0.005,1,0.0005)
	m2_ax=np.arange(0.005,3,0.0005)
	(cosi,m2)=np.meshgrid(cosi_ax,m2_ax)

	#Recompute the m1 grid out of m2 and cosi, but with the same dimensionality.
	m1=(m2**3*(6.67408e-11*1.9891e30/(4*np.pi**2))*((np.sqrt(1-cosi**2)**3)*(p_orb**2)/(x*299792458)**3))**(1/2)-m2

	#Evaluate m2 and mtot at each point of the m1-m2 grid and compute the distribution.
	M2=m2
	MTOT=m1+m2
	cosi_m2_distribution=normal_distribution(M2,nominal_m2,dm2)*normal_distribution(MTOT,nominal_mT,dmT)

	return m1_m2_distribution,cosi_m2_distribution,one_sigma

#@jit(nopython=True)
def find_sigma_contours(array):

	#Compute the total contents of the array, and the maximum value.
	total=np.sum(array)
	max_val=np.max(array)

	#Find the 1-sigma, 2-sigma and 3-sigma contour values.
	one_sigma=total*(1-0.6827)
	two_sigma=total*(1-0.9545)
	three_sigma=total*(1-0.9973)
	sigmas=[three_sigma,two_sigma,one_sigma]

	#Loop over accumulaed values.
	fractions=10**np.arange(-0.5,3.003,0.002)
	i=1
	accumulated_value=np.zeros(shape=len(fractions))
	contour_values=[]
	for fraction in fractions:
		accumulated_value[i-1] = np.sum(array[ array < (max_val/1000)*fraction ])
		if i/50==int(i/50):
			print(i,(max_val/1000)*fraction)
		i=i+1
	accumulated_value=np.array(accumulated_value)

	#Find the contour values.
	for sigma in sigmas:
		argument_of_closest=np.argmin(abs(accumulated_value-sigma))
		contour_values.append((max_val/1000)*fractions[argument_of_closest])

	print(contour_values)
	return contour_values

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
parser.add_argument("-r","--range",help="Range of Shapiro delay in s. It supersedes h3.")
parser.add_argument("-s","--shape",help="Shape of Shapiro delay. It supersedes stig.")
parser.add_argument("--h3",help="Orthometric amplitude amplitude of saphiro delay (value+/-uncertainty) in s.")
parser.add_argument("--stig",help="Orthometric amplitude amplitude of saphiro delay (value+/-uncertainty).")
parser.add_argument("--m2",help="Companion mass from the DDGR model (value+/-uncertainty) in solar masses. It will be plotted in the diagrams.")
parser.add_argument("--mtot",help="Total mass from the DDGR model (value+/-uncertainty) in solar masses. It will be used for plotting cosi and M1.")
parser.add_argument("--plot_product",type=bool,help="If set, compute and plot a 3-sigma contour from the product of all post-Keplerian distributions.",default=False)
parser.add_argument("-v","--verbose",action="store_true")
args = parser.parse_args()

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
	omdot_significance=int(100*omdot/domdot)
	if args.verbose==True:
		print("Constraining from omdot {} º/yr.".format(args.omdot))
	(omdot_m1_m2,omdot_cosi_m2,omdot_one_sigma)=periastron_advance_constraint(p_orb,ecc,x,omdot,domdot)

if args.gamma:
	gamma=float(args.gamma.split("+/-")[0])
	dgamma=float(args.gamma.split("+/-")[1])
	gamma_significance=int(100*gamma/dgamma)
	if args.verbose==True:
		print("Constraining from gamma {} s.".format(args.gamma))
	(gamma_m1_m2,gamma_cosi_m2,gamma_one_sigma)=constraint_from_gamma(p_orb,ecc,x,gamma,dgamma)

if args.h3:
	h3=float(args.h3.split("+/-")[0])
	dh3=float(args.h3.split("+/-")[1])
	h3_significance=int(100*h3/dh3)
	if args.verbose==True:
		print("Constraining from Shapiro delay (h3).")
	(h3_m1_m2,h3_cosi_m2,h3_one_sigma)=constraint_from_h3(p_orb,x,h3,dh3)

if args.stig:
	stig=float(args.stig.split("+/-")[0])
	dstig=float(args.stig.split("+/-")[1])
	stig_significance=int(100*stig/dstig)
	if args.verbose==True:
		print("Constraining from Shapiro delay (stig).")
	(stig_m1_m2,stig_cosi_m2,stig_one_sigma)=constraint_from_stig(p_orb,x,stig,dstig)

if args.m2 and args.mtot:
	m2_value=float(args.m2.split("+/-")[0])
	dm2=float(args.m2.split("+/-")[1])
	mtot_value=float(args.mtot.split("+/-")[0])
	dmtot=float(args.mtot.split("+/-")[1])
	if args.verbose==True:
		print("Plotting the DDGR model masses.")
	(DDGR_m1_m2,DDGR_cosi_m2,DDGR_one_sigma)=constraint_from_DDGR(p_orb,x,m2_value,mtot_value,dm2,dmtot)


# Start plotting constraints like a madman:
m1_ax=np.arange(1.0,2.4,0.0005)
m2_ax=np.arange(0.005,3,0.0005)
(m1,m2)=np.meshgrid(m1_ax,m2_ax)
total_m1_m2=(m1/m1+m2/m2)/2

if args.omdot:

	plt.contourf(m1,m2,omdot_m1_m2,[omdot_one_sigma,10*omdot_one_sigma],colors=['khaki'],zorder=omdot_significance,alpha=0.4)
	plt.contour(m1,m2,omdot_m1_m2,[omdot_one_sigma],colors=['khaki'],zorder=omdot_significance)
	plt.plot([],[],"yo",label="$\dot\omega$")
	total_m1_m2=total_m1_m2*omdot_m1_m2

if args.gamma:

	plt.contourf(m1,m2,gamma_m1_m2,[gamma_one_sigma,10*gamma_one_sigma],colors=['lightcoral'],zorder=gamma_significance,alpha=0.4)
	plt.contour(m1,m2,gamma_m1_m2,[gamma_one_sigma],colors=['lightcoral'],zorder=gamma_significance)
	plt.plot([],[],"ro",label="$\gamma$")
	total_m1_m2=total_m1_m2*gamma_m1_m2

if args.h3:

	plt.contourf(m1,m2,h3_m1_m2,[h3_one_sigma,10*h3_one_sigma],colors=['lightskyblue'],zorder=h3_significance,alpha=0.4)
	plt.contour(m1,m2,h3_m1_m2,[h3_one_sigma],colors=['lightskyblue'],zorder=h3_significance)
	plt.plot([],[],"bo",label="$h_3$")	
	total_m1_m2=total_m1_m2*h3_m1_m2

if args.stig:
	plt.contourf(m1,m2,stig_m1_m2,[stig_one_sigma,10*stig_one_sigma],colors=['lightgreen'],zorder=stig_significance,alpha=0.4)
	plt.contour(m1,m2,stig_m1_m2,[stig_one_sigma],colors=['lightgreen'],zorder=stig_significance)
	plt.plot([],[],"go",label="$\\varsigma$")
	total_m1_m2=total_m1_m2*stig_m1_m2

if args.m2 and args.mtot:
#	plt.contour(m1,m2,DDGR_m1_m2,[DDGR_one_sigma],colors=['k'],zorder=99999999999999999999999999999)
	plt.plot([mtot_value-m2_value],[m2_value],"ko",label="DDGR",zorder=99999999999999999999999999999)

#Plot the constraints on the mass up to two sigma
if args.plot_product==True:
	plt.contour(m1,m2,total_m1_m2,find_sigma_contours(total_m1_m2),colors=['lightgrey','darkgrey','grey'],zorder=9999999999999999999999999999)

# Plot the mass function limit.
m2=np.arange(0.005,3,0.0005)
limit=(m2**3*(6.67408e-11*1.9891e30/(4*np.pi**2))*(((p_orb*24*3600)**2)/(x*299792458)**3))**(1/2)-m2
m1=limit
plt.fill_between(m1,m2,color="darkgray",zorder=99999999)

# Plot the total mass distribution

plt.xlabel("Pulsar mass (M$_\odot$)")
plt.ylabel("Companion mass (M$_\odot$)")
plt.xlim((1.0,2.4))
plt.ylim((0.0,3.0))
plt.grid()
plt.legend()
plt.show()

# Start plotting constraints like a madman:
cosi_ax=np.arange(0.005,1,0.0005)
m2_ax=np.arange(0.005,3,0.0005)
(cosi,m2)=np.meshgrid(cosi_ax,m2_ax)
total_cosi_m2=(cosi/cosi+m2/m2)/2

if args.omdot:

	plt.contourf(cosi,m2,omdot_cosi_m2,[omdot_one_sigma,10*omdot_one_sigma],colors=['khaki'],zorder=omdot_significance,alpha=0.4)
	plt.contour(cosi,m2,omdot_cosi_m2,[omdot_one_sigma],colors=['khaki'],zorder=omdot_significance)
	plt.plot([],[],"yo",label="$\dot\omega$")
	total_cosi_m2=total_cosi_m2*omdot_cosi_m2

if args.gamma:

	plt.contourf(cosi,m2,gamma_cosi_m2,[gamma_one_sigma,10*gamma_one_sigma],colors=['lightcoral'],zorder=gamma_significance,alpha=0.4)
	plt.contour(cosi,m2,gamma_cosi_m2,[gamma_one_sigma],colors=['lightcoral'],zorder=gamma_significance)
	plt.plot([],[],"ro",label="$\gamma$")
	total_cosi_m2=total_cosi_m2*gamma_cosi_m2

if args.h3:

	plt.contourf(cosi,m2,h3_cosi_m2,[h3_one_sigma,10*h3_one_sigma],colors=['lightskyblue'],zorder=h3_significance,alpha=0.4)
	plt.contour(cosi,m2,h3_cosi_m2,[h3_one_sigma],colors=['lightskyblue'],zorder=h3_significance)
	plt.plot([],[],"bo",label="$h_3$")	
	total_cosi_m2=total_cosi_m2*h3_cosi_m2

if args.stig:

	plt.contourf(cosi,m2,stig_cosi_m2,[stig_one_sigma,10*stig_one_sigma],colors=['lightgreen'],zorder=stig_significance,alpha=0.4)
	plt.contour(cosi,m2,stig_cosi_m2,[stig_one_sigma],colors=['lightgreen'],zorder=stig_significance)
	plt.plot([],[],"go",label="$\\varsigma$")
	total_cosi_m2=total_cosi_m2*stig_cosi_m2

if args.m2 and args.mtot:
#	plt.contour(cosi,m2,DDGR_cosi_m2,[DDGR_one_sigma],colors=['k'],zorder=99999999999999999999999999999)
	sini_value=299792458*x*((1/(p_orb*24*3600)**2)*(4*np.pi**2/6.67408e-11)*(((mtot_value)**2)/(1.9891e30*m2_value**3)))**(1/3)
	cosi_value=np.sqrt(1-sini_value**2)
	plt.plot([cosi_value],[m2_value],"ko",label="DDGR",zorder=99999999999999999999999999999)

if args.plot_product==True:
	plt.contour(cosi,m2,total_cosi_m2,find_sigma_contours(total_cosi_m2),colors=['lightgrey','darkgrey','grey'],zorder=9999999999999999999999999999)

plt.xlabel("cos(i)")
plt.ylabel("Companion mass (M$_\odot$)")
plt.xlim((0,1))
plt.ylim((0.0,3.0))
plt.grid()
plt.legend()
plt.show()
	