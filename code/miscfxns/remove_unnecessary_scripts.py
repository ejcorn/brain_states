import os
import re
import glob
import numpy as np
import itertools

homedir = '/Users/Eli/Dropbox/Cornblath_Bassett_Projects/BrainStateTransitions/brain_states/'
pipeline_file = 'finalmain.sh'
fin = open(homedir+pipeline_file) # load pipeline
pipeline=fin.readlines() # into python list

def unlist(x):
	return(list(itertools.chain.from_iterable(x)))
def setdiff(list1,list2):
	return list(set(list1) - set(list2))
def contains(l,pattern):
	# analagous to matlabs contains() function or R grepl() function. finds pattern in each element of a list
	return([pattern in str(x) for x in l])
# find all .sh scripts called in pipeline
used_jobs = [re.findall('[a-zA-Z0-9_]*\.sh',s) for s in pipeline]
used_jobs = np.squeeze([x for x in used_jobs if len(x) > 0])
used_jobs = [x for x in used_jobs if os.path.exists(homedir+'jobs/' + x)] # remove if not in jobs folder
all_jobs = [os.path.basename(x) for x in glob.glob(homedir+'jobs/*')]
unused_jobs = setdiff(all_jobs,used_jobs)
[os.system('mv ' + homedir + 'jobs/'+f + ' '+ homedir + 'all_jobs_unused/'+f) for f in unused_jobs] # delete unused job scripts

# load those .sh scripts
job_files = [open(homedir+'jobs/'+j).readlines() for j in used_jobs]
# loop through those .sh scripts and find the .m, .R and .py files called
used_scripts = unlist([[re.findall('[a-zA-Z0-9_-]*\.[mRpy]',s) for s in jf] for jf in job_files])
used_scripts = np.squeeze([x for x in used_scripts if len(x)>0])
all_scripts = glob.glob(homedir + 'code/*/*.*')
all_scripts_unused_mask = np.array([os.path.basename(x) not in used_scripts for x in all_scripts]) # check if each of the scripts is used

unused_scripts_full = [all_scripts[x] for x in np.where(all_scripts_unused_mask)[0]]
used_scripts_full = [all_scripts[x] for x in np.where(~all_scripts_unused_mask)[0]]

# next need to load every used script and check if the unused scripts are called within them 
# (in the case of a script being a function)

def find_functions_in_scripts(paths_to_functions,paths_to_scripts):
	# INPUTS:
	# paths_to_functions: list of paths to scripts that might be called as functions
	# within other scripts
	# paths_to_scripts: list of paths to scripts that might contain these functions
	#
	# OUTPUTS:
	# (paths_to_used,paths_to_unused): split paths_to_functions into those called and not called 
	# in the files found at paths_to_scripts

	used_script_files = [open(x).readlines() for x in paths_to_scripts]
	# list element for each script, whose elements contain a nested list that should have somdething in it
	# if that script 
	used_script_files_script_search  = list()
	for script in paths_to_functions:
		# in MATLAB, scripts are called using their basename, i.e. addpaths.m called as addpaths, so search for that
		potential_fxn = os.path.basename(script).split('.')[0]
		# search through first layer of used scripts that are called in .sh file
		mentions_in_used_scripts = unlist([[re.findall(potential_fxn,s) for s in sf] for sf in used_script_files])	
		used_script_files_script_search.append(len([x for x in mentions_in_used_scripts if len(x) > 0]))
	used_script_files_script_search = np.array(used_script_files_script_search)

	# return paths to scripts which are called as functions at least once across all the scripts
	paths_to_used = [paths_to_functions[x] for x in np.where(used_script_files_script_search >0)[0]]
	# return paths to scripts which are NOT called as functions at least once across all the scripts
	paths_to_unused = [paths_to_functions[x] for x in np.where(used_script_files_script_search ==0)[0]]
	return paths_to_used,paths_to_unused

paths_to_used,paths_to_unused = find_functions_in_scripts(unused_scripts_full,used_scripts_full)
# add the initially unused scripts that are called as functions back to used script list
used_scripts_full = used_scripts_full + paths_to_used 
# run again to see if there are any functions within functions
paths_to_used,paths_to_unused = find_functions_in_scripts(paths_to_unused,used_scripts_full)
# rinse and repeat until no more
used_scripts_full = used_scripts_full + paths_to_used 
paths_to_used,paths_to_unused = find_functions_in_scripts(paths_to_unused,used_scripts_full)
# now paths_to_used is empty
# move all files in paths_to_unused

[os.system('mkdir -p ' + homedir + 'all_code_unused/' + x) for x in os.listdir('code')]
[os.system('mv ' +f + ' '+ homedir + 'all_code_unused/'+f.split('/')[-2] + '/'+f.split('/')[-1]) for f in paths_to_unused] # move unused job scripts

