import numpy as np
import subprocess
from os.path import exists
import argparse
import glob
import multiprocessing as multi
from functools import partial
from contextlib import closing


# libstempo is unable to read CHI2R from the PAR, so I did it myself.
def read_chi2r_from_par(parFile):
	par_read=open(parFile,"r")
	for line in par_read:
		if line == "" or line == " " or line == "	": #Fixing minor bugs due to empty lines.
			line="1 1"
		if line == "\n": #Fixing minor bugs due to empty lines.
			line="1 1\n"
		chunks = line.strip().split()
		if chunks[0]=="CHI2R":
			chi2r=float(chunks[1])
	par_read.close()
	return chi2r

# Read a higher resolution version of the chi2r.
def read_chi2r_from_logs(tempo2logs):
	for line in tempo2logs.split("\n"):
		if line[0:11]=="Fit Chisq =":
			chi2r=float(line.split("=")[3].split("	")[0])
			break
	return chi2r

# libstempo is unable to read CHI2R from the PAR, so I did it myself.
def compatible(parFile,ref_pos,tol):

	import astropy.units as u
	from astropy.coordinates import SkyCoord

	compatible=True

	par_read=open(parFile,"r")
	for line in par_read:
		if line == "" or line == " " or line == "	": #Fixing minor bugs due to empty lines.
			line="1 1"
		if line == "\n": #Fixing minor bugs due to empty lines.
			line="1 1\n"
		chunks = line.strip().split()
		if chunks[0]=="RAJ":
			ra=chunks[1]
			dra=chunks[3]
		if chunks[0]=="DECJ":
			dec=chunks[1]
			ddec=chunks[3]
	par_read.close()

	# Measure distance between sources.

#	print(ra,dec)

	solution = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
	reference = SkyCoord(ref_pos.split(",")[0], ref_pos.split(",")[1], unit=(u.hourangle, u.deg), frame='icrs')

	print(solution.arcsec)
	print(reference.arcsec)

#	distance = solution.separation(reference)

	if ( abs(solution.arcsec[0]-reference.arcsec[0])-dra > tolerance/np.cos(reference.rad[0]) ) or ( abs(solution.arcsec[1]-reference.arcsec[1])-ddec > tolerance ):

		compatible=False

	return compatible

# Read a PAR file, add JUMPS MJD statements for each osbervation except the last one,
# and fit all of the parameters (with 1 in the PAR) and the JUMP MJD statements.
# This marks the initialization of the script as well.
def add_jumps_and_fit(parFile,timFile,skipJumps,nFits):

	# Open files and create a new PAR file that will contain the jumps.
	par_read=open(parFile,"r")
	parFile_jumps=parFile.split(".")[0]+"_jumps.par"
	par_write=open(parFile_jumps,"w")
	for line in par_read:
		par_write.write(line)
	tim_read=open(timFile,"r")
	old_obs=""
	i=0
	jumps=0
	start_end_pairs=[]
	time_intervals=[]
	phase_jumps_times=[]

	# Compute the jumps.
	print("Reading ToAs and computing jumps.")
	par_write.write("\n")
	for line in tim_read:
		if line=="FORMAT 1\n" or line=="MODE 1\n" or line[0]=="C":
			dummy=1
		else:
			line=line.strip().split()
			obs=line[0]
			time=float(line[2])
			if i==0:
				time_start=time
			if obs!=old_obs and i>1:
				time_end=time_old
				start_end_pairs.append([time_start-0.001,time_end+0.001])
				time_start=time
				time_intervals.append(time_start-time_end)
				phase_jumps_times.append((time_end+time_start)/2)
				jumps=jumps+1
			old_obs=obs
			time_old=float(time)
			i=i+1

	nToAs=i

	par_read.close()
	tim_read.close()

	# Add the time jumps.
	# Unlike in tempo, tempo2 can not jump at a particular MJD, but it can jump the data within an interval [MJD1,MJD2].
	# We put each observation except the last one within its own jump interval.
	# The last observation, which has no MJD interval jump is therefore the reference one.
	# This is why it is important to not jump it. However, it can have a jump if it is not a MJD one and stays fiexd (e.g., an instrumental jump). 
	print("")
	print("Adding time jumps to "+parFile_jumps)

	start_end_pairs=np.array(start_end_pairs)
	time_intervals=np.array(time_intervals)
	phase_jumps_times=np.array(phase_jumps_times)

	i=0
	for pair in start_end_pairs:
		print("JUMP"+str(i)+" MJD "+str(pair[0]-0.001)+" "+str(pair[1]+0.001)+" 0.0 1")
		if skipJumps!=True:
			par_write.write("JUMP MJD "+str(pair[0]-0.001)+" "+str(pair[1]+0.001)+" 0.0 1\n")
		i=i+1

	par_write.close()

	# Load the PAR with jumps and fit for everything.
	print("")
	print("Fitting "+parFile.split(".")[0]+"_jumps.par with "+timFile)

	for i in range(1,nFits+1):
		subprocess.run(["tempo2","-f",parFile_jumps,timFile,"-outpar",parFile_jumps],stdout=subprocess.DEVNULL)

	print(" ")
	print("Number of jumps:",jumps)
	print("Number of ToAs:",nToAs)
	chi2r=read_chi2r_from_par(parFile_jumps)
	print("CHI2R:",chi2r)

	return jumps,chi2r,time_intervals,phase_jumps_times

# This removes one time jump and adds a phase jump in the PAR.
def remove_jump_add_phase(parFile_jumps,phase_jump_times,jump_index,phase,max_chi2r):

	par_read=open(parFile_jumps,"r")
	if phase>=0:
		parFile_phases=parFile_jumps.split(".")[0]+"_"+str(jump_index)+"+"+str(phase)+".par"
	if phase<0:
		parFile_phases=parFile_jumps.split(".")[0]+"_"+str(jump_index)+str(phase)+".par"
	par_write=open(parFile_phases,"w")

	# Read the array of time jumps.
	print("")
	print("Removing JUMP"+str(jump_index)+" (MJD= "+str(phase_jumps_times[jump_index])+")")
	jumps=[]
	for line in par_read:
		if line == "" or line == " " or line == "	": #Fixing minor bugs due to empty lines.
			line="1 1"
		if line == "\n": #Fixing minor bugs due to empty lines.
			line="1 1\n"
		chunks = line.strip().split()
		if chunks[0]+" "+chunks[1]=="JUMP MJD":
			jumps.append([float(chunks[2]),float(chunks[3]),float(chunks[4])])
		else:
			par_write.write(line)
	jumps=np.array(jumps)
	# Create a new array with removed time jumps according to phase jumps position.
	# Unlike in tempo, we can't remove jumps between observations in tempo2.
	# Instead, we join consecutive [MJD1,MJD2] jump intervals to join 2 observations together in the same jump with respect the last observation.
	i=1
	skip_step=False
	for element in jumps:
		if skip_step==True:
			skip_step=False
			i=i+1
			continue
		# The second-to-last observation (last [MJD1,MJD2] jump intervals) can only be joined with the last obs, which has no jump.
		# Therefore, in such case, the [MJD1,MJD2] jump interval is just removed. This is repeated for every jump that comes lust in time ordering.
		if (element==jumps[-1]).all():
			if phase_jump_times[jump_index]>element[1]:
				dummy=1 #Do not write the time jump in the new file.
			else:
				par_write.write("JUMP MJD "+str(element[0])+" "+str(element[1])+" "+str(element[2])+" 1\n")
		# If we don't have the first or last [MJD1,MJD2] jump intervals, then we can start joining them.
		else:
			# Check if studied phase jump lies in between to time jumps. If so, join them. Also, skip the next step as it is part of the joined interval.
			if phase_jump_times[jump_index]>element[1] and phase_jump_times[jump_index]<jumps[i,0]:
				par_write.write("JUMP MJD "+str(element[0])+" "+str(jumps[i,1])+" 0.0 1\n")
				element=par_read.readline()
				skip_step=True
			# Otherwise, leave them untouched.
			else:
				par_write.write("JUMP MJD "+str(element[0])+" "+str(element[1])+" "+str(element[2])+" 1\n")
				element_old=element
		i=i+1

	par_read.close()

	# Add the desired phase jump.
	if phase>=0:
		print("Adding phase turn: +"+str(phase))
		par_write.write("PHASE +"+str(phase)+" "+str(phase_jump_times[jump_index])+"\n")
		print("PHASE +"+str(phase)+" "+str(phase_jump_times[jump_index]))
	if phase<0:
		print("Adding phase turn:",str(phase))
		par_write.write("PHASE "+str(phase)+" "+str(phase_jump_times[jump_index])+"\n")
		print("PHASE "+str(phase)+" "+str(phase_jump_times[jump_index]))

	par_write.close()

	# Fit! Also, store the output to read the chi2r, whih has a higher resolution thatn the one in the parameter file.
	subrun=subprocess.run(["tempo2","-f",parFile_phases,timFile,"-outpar",parFile_phases.split(".")[0]+"_solution_candidate.par"],stdout=subprocess.PIPE)

	# Check if the fit has not crashed (solution_candidate exists.)
	solution_exists=exists(parFile_phases.split(".")[0]+"_solution_candidate.par")

	# Read the resulting chi2. Preserve the par file only if chi2r<2.
	subprocess.run(["mv",parFile_phases.split(".")[0]+"_solution_candidate.par",parFile_phases],stdout=subprocess.DEVNULL)
	chi2r=read_chi2r_from_par(parFile_phases) #Preliminarily, read the chi2r from the parameter file
	if chi2r<max_chi2r and solution_exists==True:
		chi2r=read_chi2r_from_logs(subrun.stdout.decode()) #Get the higher res version of the chi2r
		# We remove the phase jump for the future use of the par file.
		par_read=open(parFile_phases,"r")
		par_write=open(parFile_phases+"a","w")
		for line in par_read:
			chunks = line.strip().split()
			if chunks[0]!="PHASE":
				par_write.write(line)
		par_read.close()
		par_write.close()
		subprocess.run(["mv",parFile_phases+"a",parFile_phases],stdout=subprocess.DEVNULL)
#		subprocess.run(["tempo2","-f",parFile_phases,timFile,"-outpar",parFile_phases],stdout=subprocess.DEVNULL)
		print("CHI2R=",chi2r)
	elif solution_exists==True:
		chi2r=read_chi2r_from_logs(subrun.stdout.decode()) #Get the higher res version of the chi2r
		subprocess.run(["rm",parFile_phases],stdout=subprocess.DEVNULL)
		print("CHI2R=",chi2r)
	else:
		subprocess.run(["rm",parFile_phases],stdout=subprocess.DEVNULL)

	return chi2r,solution_exists

#Identify from 3 chi2r values whether we are in a minimum or a slope, and then run until the valley within chi2r<2 is identified.
def find_chi2r_interval(parFile,phase_jump_times,jump_index,max_chi2r,max_solutions):

	print("")
	print("Computing ramifications from solution "+parFile+".")

	if max_solutions/2 == int(max_solutions/2): #Make sure that the maximum number of solutions is even (symmetry in parabola, including the middle point).
		max_solutions=max_solutions+1
	chi2r=[]
	phase=[-1,0,1]

	(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,-1,max_chi2r)
	if exists==True:
		chi2r.append(instant_chi2r)
		left_chi2r=instant_chi2r
	else:
		print("Tempo2 can't fit a solution for this phase jump. Using dummy CHI2R value of 999999999.")
		chi2r.append(999999999.0)
		left_chi2r=999999999.0
	(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,0,max_chi2r)
	if exists==True:
		chi2r.append(instant_chi2r)
		middle_chi2r=instant_chi2r
	else:
		print("Tempo2 can't fit a solution for this phase jump. Using dummy CHI2R value of 999999999.")
		chi2r.append(999999999.0)
		middle_chi2r=999999999.0
	(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,+1,max_chi2r)
	if exists==True:
		chi2r.append(instant_chi2r)
		right_chi2r=instant_chi2r
	else:
		print("Tempo2 can't fit a solution for this phase jump. Using dummy CHI2R value of 999999999.")
		chi2r.append(999999999.0)
		right_chi2r=999999999.0

	# Check which component presents the lowest chi2r.
	if chi2r[0]!=chi2r[1] and chi2r[1]!=chi2r[2]:
		direction=chi2r.index(min(chi2r))
		starting_phase=2
		print("")
		print("The slope is clear from the three chi2r values.")

	# Make sure that we don't have a disaster here (erase ambiguity in case of similar numbers).
	if (chi2r[0]==chi2r[1] and chi2r[1]<chi2r[2]) or (chi2r[0]<chi2r[1] and chi2r[1]==chi2r[2]):
		direction=0
		starting_phase=2
		print("")
		print("Some ambiguity in the chi2r values, but a prefered direction is found.")
		max_solutions=max_solutions+2

	if (chi2r[1]==chi2r[2] and chi2r[1]<chi2r[0]) or (chi2r[1]>chi2r[2] and chi2r[1]==chi2r[0]):
		direction=2
		starting_phase=2
		print("")
		print("Some ambiguity in the chi2r values, but a prefered direction is found.")
		max_solutions=max_solutions+2

	if chi2r[1]==chi2r[2] and chi2r[0]==chi2r[1]:   #Special case. Perhaps we are on a very smooth slope. It requires some thought.

		print("")
		print("All chi2r values are the same. Doing some further tests.")

		(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,-2,max_chi2r)
		if exists==True:
			phase.insert(0,-2)
			chi2r.insert(0,instant_chi2r)
			left_chi2r=instant_chi2r
		else:
			print("Tempo2 can't fit a solution for this phase jump. Using dummy CHI2R value of 999999999.")
			phase.insert(0,-2)
			chi2r.insert(0,999999999.0)
			left_chi2r=999999999.0

		(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,2,max_chi2r)
		if exists==True:
			phase.append(2)
			chi2r.append(instant_chi2r)
			right_chi2r=instant_chi2r
		else:
			print("Tempo2 can't fit a solution for this phase jump. Using dummy CHI2R value of 999999999.")
			phase.append(2)
			chi2r.append(999999999.0)
			left_chi2r=999999999.0

		if (chi2r[0]>chi2r[1] and chi2r[4]>chi2r[3]) or (chi2r[0]==chi2r[1] and chi2r[4]==chi2r[3]): #If it's all the same for 5 in a row, we must be in hell of a deep minimum.

			direction=1
			print("")
			print("Minimum found or ambiguity persists. Assuming that we are at a minimum.")
			max_solutions=max_solutions+2

		else:

			dummy_array=[left_chi2r,999999999.0,right_chi2r]
			direction=dummy_array.index(min(dummy_array))
			print("")
			print("A prefered direction has been found.")
			max_solutions=max_solutions+4

		starting_phase=3

	# Choose a direction based on this.

	if direction==0: #We look for the minimum of the parabola on the left.

		print("Looking for the minimum at negative phase turns.")

		# Walk backwards until we hit the chi2r or the max solutions wall.
		i=-starting_phase
		instant_chi2r=left_chi2r
		previous_chi2r=middle_chi2r
		min_phase=-9999999999 #Just a dummy value that makes sure we don't have min_phase-i > (max_solutions-1)/2 before the minimum is found.
		min_phase_set=False
		while (instant_chi2r<max_chi2r or (instant_chi2r>max_chi2r and (previous_chi2r-instant_chi2r)>=0)) and min_phase-i <= (max_solutions-1)/2: #While we haven't reached the max chi2r value, or we are going down in chi2r value, or we don't have the max amount of wanted solutions on one side.
			previous_chi2r=instant_chi2r
			(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,i,max_chi2r)
			if exists==False:
				print("Tempo2 can't fit a solution beyond this point.")
				break
			if (previous_chi2r-instant_chi2r)<0 and min_phase_set==False: #This means we have found the minimum!
				min_phase=i+1
				min_phase_set=True
				print("")
				print("Minimum found at phase turn",min_phase)
			phase.insert(0,i)
			chi2r.insert(0,instant_chi2r)
			i=i-1

		# Walk forward until we hit the chi2r wall once again.
		i=starting_phase
		instant_chi2r=right_chi2r
		previous_chi2r=middle_chi2r
		while (instant_chi2r<max_chi2r or (instant_chi2r>max_chi2r and (instant_chi2r-previous_chi2r)<=0)) and i-min_phase <= (max_solutions-1)/2:
			previous_chi2r=instant_chi2r
			(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,i,max_chi2r)
			if exists==False:
				print("Tempo2 can't fit a solution beyond this point.")
				break	
			phase.append(i)
			chi2r.append(instant_chi2r)
			i=i+1

	if direction==2: #We look for the minimum of the parabola on the right. 

		print("Looking for the minimum at positive phase turns.")

		# Walk forward until we hit the chi2r or the max solutions wall.
		i=starting_phase
		instant_chi2r=right_chi2r
		previous_chi2r=middle_chi2r
		min_phase=9999999999
		min_phase_set=False
		while (instant_chi2r<max_chi2r or (instant_chi2r>max_chi2r and (instant_chi2r-previous_chi2r)<=0)) and i-min_phase <= (max_solutions-1)/2:
			previous_chi2r=instant_chi2r
			(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,i,max_chi2r)
			if exists==False:
				print("Tempo2 can't fit a solution beyond this point.")
				break
			if (previous_chi2r-instant_chi2r)<0 and min_phase_set==False:
				min_phase=i-1
				min_phase_set=True
				print("")
				print("Minimum found at phase turn",min_phase)
			phase.append(i)
			chi2r.append(instant_chi2r)
			i=i+1

		# Walk forward until we hit the chi2r wall once again.
		i=-starting_phase
		instant_chi2r=left_chi2r
		previous_chi2r=middle_chi2r
		while (instant_chi2r<max_chi2r or (instant_chi2r>max_chi2r and (previous_chi2r-instant_chi2r)>=0)) and min_phase-i <= (max_solutions-1)/2:
			previous_chi2r=instant_chi2r
			(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,i,max_chi2r)
			if exists==False:
				print("Tempo2 can't fit a solution beyond this point.")
				break
			phase.insert(0,i)
			chi2r.insert(0,instant_chi2r)
			i=i-1

	if direction==1: #The minimum of the oparabola is already known.

		print("Phase turn 0 (no turns) is the minimum.")

		# Walk forward until we hit the chi2r or the max solutions wall.
		i=starting_phase
		instant_chi2r=right_chi2r
		previous_chi2r=middle_chi2r
		while (instant_chi2r<max_chi2r or (instant_chi2r>max_chi2r and (instant_chi2r-previous_chi2r)<=0)) and i <= (max_solutions-1)/2:
			previous_chi2r=instant_chi2r
			(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,i,max_chi2r)
			if exists==False:
				print("Tempo2 can't fit a solution beyond this point.")
				break
			phase.append(i)
			chi2r.append(instant_chi2r)
			i=i+1

		# Walk forward until we hit the chi2r wall once again.
		i=-starting_phase
		instant_chi2r=left_chi2r
		previous_chi2r=middle_chi2r
		while (instant_chi2r<max_chi2r or (instant_chi2r>max_chi2r and (previous_chi2r-instant_chi2r)>=0)) and -i <= (max_solutions-1)/2:
			previous_chi2r=instant_chi2r
			(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,i,max_chi2r)
			if exists==False:
				print("Tempo2 can't fit a solution beyond this point.")
				break
			phase.insert(0,i)
			chi2r.insert(0,instant_chi2r)
			i=i-1

	#This is left-over code from when the maximum amount of solutions wasn't accounted for in the explorer. This could potentially produce hundreths of solutions that were later deleted at the expense of computing time.

	# Walk forward to make sure we are within the range.
#	i=2
#	instant_chi2r=chi2r[2]
#	previous_chi2r=chi2r[1]
#	while instant_chi2r<max_chi2r or (instant_chi2r>max_chi2r and (instant_chi2r-previous_chi2r)<0):
#		previous_chi2r=instant_chi2r
#		(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,i,max_chi2r)
#		if exists==False:
#			print("Tempo2 can't fit a solution beyond this point.")
#			break	
#		phase.append(i)
#		chi2r.append(instant_chi2r)
#		i=i+1

	# Walk back to make sure we are within the range.
#	i=-2
#	instant_chi2r=chi2r[0]
#	previous_chi2r=chi2r[1]
#	while instant_chi2r<max_chi2r or (instant_chi2r>max_chi2r and (previous_chi2r-instant_chi2r)>0):
#		previous_chi2r=instant_chi2r
#		(instant_chi2r,exists)=remove_jump_add_phase(parFile,phase_jump_times,jump_index,i,max_chi2r)
#		if exists==False:
#			print("Tempo2 can't fit a solution beyond this point.")
#			break
#		phase.insert(0,i)
#		chi2r.insert(0,instant_chi2r)
#		i=i-1

	chi2r=np.array(chi2r)
	phase=np.array(phase)

	print("")
	print("At most "+str(max_solutions)+" solutions with CHI2R<"+str(max_chi2r)+" are kept and sent to the next jump removal.")

	# Remove the remaining unnecessary files.
	sorting=np.argsort(chi2r)
	phase=phase[sorting]
	phases_to_delete=phase[max_solutions:]
	i=0
	for phase_d in phases_to_delete:
		if phase_d>=0:
			subprocess.run(["rm",parFile.split(".")[0]+"_"+str(jump_index)+"+"+str(phase_d)+".par"],stdout=subprocess.DEVNULL)
		if phase_d<0:
			subprocess.run(["rm",parFile.split(".")[0]+"_"+str(jump_index)+str(phase_d)+".par"],stdout=subprocess.DEVNULL)
		i=i+1

	print("Phase jumps:",phase[:max_solutions])
	print("CHI2R values:",chi2r[sorting[:max_solutions]])

	return 1

def order_and_trim_files(max_files,parFiles):

	chi2r_list=[]

	for file in parFiles:

		chi2r_list.append(read_chi2r_from_par(file))

	parFiles=np.array(parFiles)
	chi2r_list=np.array(chi2r_list)

	parFiles=parFiles[np.argsort(chi2r_list)]
	parFiles=parFiles[0:max_files]

	return parFiles

def remove_by_position(parFiles,ref_pos,tolerance):

	distance_list=[]
	uncertainty_list=[]

	for file in parFiles:

		(distance,uncertainty)=compute_distance(file,ref_pos)
		distance_list.append(distance)
		uncertainty_list.append(uncertainty)

		if distance-uncertainty>tolerance:
			subprocess.run(["rm",file],stdout=subprocess.DEVNULL)

	parFiles=np.array(parFiles)
	distance_list=np.array(distance_list)
	uncertainty_list=np.array(uncertainty_list)

	parFiles=parFiles[ distance-uncertainty<=tolerance ]

	return parFiles

parser=argparse.ArgumentParser(description="Take in a tempo2 parameter file and a tim file, and attempt to find a phase connection with JUMP and PHASE statements. It requires an installation of tempo2 and numpy.")
parser.add_argument("-p","--parameter",help="Tempo2 parameter file WITHOUT 'JUMP MJD' or 'PHASE' statements. There can be other kinds of jumps, but the last observation MUST either NOT be jumped, or have its jump value FIXED. For instance, if you have backend jumps, the backend with the last observation should be the non-jumped one. Only parameters with 1 will be fit.")
parser.add_argument("-t","--tim",help="Tempo2 tim file. It requires: observation name in the 1st column, and ToA in the third columns.")
parser.add_argument("--max_chi2r",type=float,help="Largest acceptable chi2r value for a solution. Default: 2.0",default=2.0)
parser.add_argument("--max_solutions",type=int,help="Largest amount of solutions that are taken from each jump removal attempt. The code will padd this number a bit is ambiguities are found in some jump removals. Default: 5",default=5)
parser.add_argument("--max_files",type=int,help="Largest global amount of files taken from each jump removal attempt, always with the lowest CHI2R value. The code will padd this number a bit is ambiguities are found in some jump removals. Default: all of them.")
parser.add_argument("--n_gulp",type=int,help="Number of jumps to remove at the same time (multithreading). The on-screen outputs become funky. Default: a single thread (serial).")
parser.add_argument("--pre_fits",type=int,help="Number of fits done to the initial file once jumps are added.",default=1)
parser.add_argument("--par_with_jumps",type=bool,help="If set, then jumps are assumed to be added manually and they are are not added by dracula2. Make sure that they are in the correct format!",default=False)
parser.add_argument("--up_to_jump",type=int,help="If set, then skip the run will end after removing this jump (counted from 0). Useful for splicing runs.")
parser.add_argument("--skip_jumps",type=int,help="If set, then skip the removal of this many jumps. Useful for continuing a broken run. MUST BE USED WITH --par_with_jumps",default=0)
parser.add_argument("--position_prior",help="Prior position of the pulsar. Solutions with a different position will be deleted. It must be used --position_tolerance")
parser.add_argument("--position_tolerance",type=float,help="Positional tolerance of found positions. It must be used with --position_prior.")

args = parser.parse_args()

parFile=args.parameter
timFile=args.tim
root=parFile.split(".")[0]

# Create a PAR file with jumps and fit for them.
(n_jumps,chi2r,time_intervals,phase_jumps_times)=add_jumps_and_fit(parFile,timFile,args.par_with_jumps,args.pre_fits)

max_jump=n_jumps-1

if args.up_to_jump:
	max_jump=args.up_to_jump

# Order the resulting time intervals to know where to start removing.
ordering=np.argsort(time_intervals)

print("")
print("Jumps will be removed in this order:",ordering)
print(phase_jumps_times[ordering])

#Loop over jumps.
i=0+args.skip_jumps
first=True
phases=[]
chi2r=[]
# At each jump removal, your directory will become VERY cluttered.
# However, all files end up nicelly collected in new folders at the end of the script.
# Ideally, only 2 extra files should be left at in your direcotory:
# - The original one with all the JUMP MJD statements.
# - The final, phase-connected solution with no JUMP MJD statements.
# If there are more than 1 possible solutions, they will be left in your folder as well.
while i<=max_jump:
	# Loop over phases
	if first==True:
		parFile=parFile.split(".")[0]+"_jumps.par"
		dummy=find_chi2r_interval(parFile,phase_jumps_times,ordering[i],args.max_chi2r,args.max_solutions)
		first=False
		i=i+1
	else:

		parFile=parFile.split(".")[0]+"_*.par"
		parFiles=glob.glob(parFile)
		print("")
		print("Surviving PAR files from removing the previous jump:",parFiles)

		if (args.position_prior and args.position_tolerance):

			for file in parFiles:

				if compatible(file,args.position_prior,args.position_tolerance)==False:

					#subprocess.run(["rm",file],stdout=subprocess.DEVNULL)
					print("This one would be removed.")

			parFile=parFile.split(".")[0]+"_*.par"

		if args.max_files:

			if len(parFiles) > args.max_files:

				print("")
				print("Selecting only the",args.max_files,"files with the lowest CHI2R value.")

				parFiles_solving=order_and_trim_files(args.max_files,parFiles)

			else:

				parFiles_solving=parFiles
		else:

			parFiles_solving=parFiles

		if args.n_gulp:

			nFiles=len(parFiles_solving)
			j=0

			while j < nFiles:

				multiprocesses=multi.Pool(processes=args.n_gulp)
				dummy_array=multiprocesses.map(partial(find_chi2r_interval,phase_jump_times=phase_jumps_times,jump_index=ordering[i],max_chi2r=args.max_chi2r,max_solutions=args.max_solutions),parFiles_solving[j:j+args.n_gulp])
				j=j+args.n_gulp

				multiprocesses.close()
				multiprocesses.join()
				j=j+1

		else:

			for file in parFiles_solving:
				dummy=find_chi2r_interval(file,phase_jumps_times,ordering[i],args.max_chi2r,args.max_solutions)

		print("")
		print("Moving files:",parFiles)
		subprocess.run(["mkdir",root+"_JUMP"+str(ordering[i-1])],stdout=subprocess.DEVNULL)
		for file in parFiles:
			subprocess.run(["mv",file,root+"_JUMP"+str(ordering[i-1])],stdout=subprocess.DEVNULL)

		i=i+1