import os 
import numpy as np
import sys
from termcolor import colored
sys.path.append('../')
from pythonlibs.process_outputs_tamcmc_library import get_aj_inc_star, get_aj_inc_star_IDL, get_ajAlm_inc_star
from pythonlibs.compute_stats import make_stats
from pythonlibs.nice_hist import nice_hist_matrix
from pythonlibs.a2CF_a2AR import do_a2CF_a2AR
from analyse_evidence import translate_modeltype2dir

def show_ajAlm_pdfs(dir_mcmc, keep_aj=None, file_out='test.jpg', binning=40, extra_gaussian=[], 
				 idlfiles=False, show_inc=False,
				 ignore_slope=True,
				process_name=None,
				phase="A", chain=0,
				first_index=0, last_index=-1, period=1,
				cpp_path="cpp_prg/", outtmpdir="tmp/", confidence=[2.25,16,50,84,97.75]):
	'''
		Basic function to show the pdfs for the activity parameters: epsilon_nl, theta0, delta
		dir_mcmc: directory that specify where to look at the data. 
		file_out: Name of the output image file
		binning: Binning of the pdf (number of classes)
		extra_gaussian: Can be used to show additional data on the top of some of the projected pdf
		idlfiles: If true, will expect reading idlfiles (gradually depreciated)
	'''
	labels=[r'$a_1$ (nHz)', r'$a_3$ (nHz)', r'$a_5$ (nHz)', r"$\epsilon\,(10^3)$", r"$\theta_0$", r"$\delta$"]
	labels_txtonly=["a1", "a3", "a5", "epsilon", "theta0", "delta"]
	units_txtonly= ["nHz","nHz","nHz","10^-3","deg","deg"]
	unit_nHz=1000 # To convert into nHz
	if idlfiles == True:
		raise ValueError("ERROR : The case idlfiles = True is not implemented in show_ajAlm_pdfs(). You have to add it or use directly the TAMCMC samples.")
	else:
		aj_stats, inc_stats, aj_samples, inc_samples, Nsamples=get_ajAlm_inc_star(dir_mcmc, confidence=confidence,				 
				ignore_slope=ignore_slope,
				process_name=process_name,
				phase=phase, chain=chain,
				first_index=first_index, last_index=last_index, period=period,
				cpp_path=cpp_path, outtmpdir=outtmpdir)
	Nparams=len(aj_stats[:,0])
	if keep_aj==None:
		keep_aj=np.repeat(True, Nparams)
		
	#
	# -- Filter according to keep_aj --
	poskeep=np.where(np.asarray(keep_aj, dtype=bool) == True)[0]
	Nparams_keep=len(poskeep)
	aj_stats_keep=np.zeros((Nparams_keep, len(confidence)))
	aj_samples_keep=np.zeros((Nparams_keep, Nsamples))
	labels_keep=[]
	# -- Apply unit conversion and filter --
	cpt=0
	for i in range(Nparams):
		if keep_aj[i] == True:
			if labels[i] != r"$\theta_0$" and labels[i] != r"$\delta$":
				aj_stats_keep[cpt,:]=aj_stats[i,:]*unit_nHz
				aj_samples_keep[cpt,:]=aj_samples[i,:]*unit_nHz
			else:
				aj_stats_keep[cpt,:]=aj_stats[i,:]
				aj_samples_keep[cpt,:]=aj_samples[i,:]
			labels_keep.append(labels[i])
			cpt=cpt+1

	if show_inc == True:
		inc_stats = inc_stats
		inc_samples = inc_samples
		aj_stats_keep = np.vstack((aj_stats_keep, inc_stats))
		aj_samples_keep = np.vstack((aj_samples_keep, inc_samples[:,0]))
		labels_keep.append("Inclination")
		Nparams_keep=Nparams_keep+1

	#### AD HOC MOD. REMOVE IF UNEXPECTED
	'''
	i = 0
	pos_keep=np.where(aj_samples_keep[i,:] > 400)[0]
	aj_samples_keep_2=np.zeros((Nparams_keep, len(pos_keep)))
	aj_stats_keep=np.zeros((Nparams_keep, len(confidence)))
	for i in range(Nparams_keep):
		aj_samples_keep_2[i,:] = aj_samples_keep[i,pos_keep]
		aj_stats_keep[i,:]=np.percentile(aj_samples_keep_2[i,:], confidence)
	aj_samples_keep=aj_samples_keep_2
	'''
	####
	# 		
	extra_data=[] # Contain all of the pdfs info that are going to be plotted. Minimalistically, extra_data[0]=extra_a2CF
	# -- Perform the plotting --
	nice_hist_matrix(aj_samples_keep, aj_stats_keep, labels_keep, binning=binning, posplots=None, file_out=file_out, extra_data=extra_data)
	#print('Saved files with syntax: ', file_out)
	for i in range(Nparams):
		if labels[i] != r"$\theta_0$" and labels[i] != r"$\delta$":
			aj_stats[i,:]=aj_stats[i,:]*unit_nHz
	return labels_txtonly, units_txtonly, aj_stats, inc_stats # We return all what was read... because it is supposed to be written later
	
def show_aj_pdfs(dir_mcmc, keep_aj=None, file_out='test.jpg', binning=40, extra_gaussian=[], 
				 idlfiles=False, show_a2CF=False, show_inc=False,
				 ignore_slope=True,
				process_name=None,
				phase="A", chain=0,
				first_index=0, last_index=-1, period=1,
				cpp_path="cpp_prg/", outtmpdir="tmp/", confidence=[2.25,16,50,84,97.75]):
	'''
		Basic function to show the pdfs for the activity parameters: epsilon_nl, theta0, delta
		dir_mcmc: directory that specify where to look at the data. 
		file_out: Name of the output image file
		binning: Binning of the pdf (number of classes)
		extra_gaussian: Can be used to show additional data on the top of some of the projected pdf
		idlfiles: If true, will expect reading idlfiles (gradually depreciated)
	'''
	labels=[r'$a_1$ (nHz)', r'$a_2$ (nHz)', r'$a_3$ (nHz)', r'$a_4$ (nHz)', r'$a_5$ (nHz)', r'$a_6$ (nHz)']
	labels_txtonly=["a1", "a2", "a3", "a4", "a5", "a6"]
	units_txtonly= ["nHz","nHz","nHz","nHz","nHz","nHz"]
	unit_nHz=1000 # To convert into nHz
	#
	if idlfiles == True:
		aj_stats, inc_stats, aj_samples, inc_samples, Nsamples=get_aj_inc_star_IDL(dir_mcmc, confidence=confidence)
	else:
		aj_stats, inc_stats, aj_samples, inc_samples, Nsamples=get_aj_inc_star(dir_mcmc, confidence=confidence,				 
				ignore_slope=ignore_slope,
				process_name=process_name,
				phase=phase, chain=chain,
				first_index=first_index, last_index=last_index, period=period,
				cpp_path=cpp_path, outtmpdir=outtmpdir)
	Nparams=len(aj_stats[:,0])
	if keep_aj==None:
		keep_aj=np.repeat(True, Nparams)
		
	#
	# -- Filter according to keep_aj --
	poskeep=np.where(np.asarray(keep_aj, dtype=bool) == True)[0]
	Nparams_keep=len(poskeep)
	aj_stats_keep=np.zeros((Nparams_keep, len(confidence)))
	aj_samples_keep=np.zeros((Nparams_keep, Nsamples))
	labels_keep=[]
	# -- Apply unit conversion and filter --
	cpt=0
	for i in range(Nparams):
		if keep_aj[i] == True:
			aj_stats_keep[cpt,:]=aj_stats[i,:]*unit_nHz
			aj_samples_keep[cpt,:]=aj_samples[i,:]*unit_nHz
			#### AD HOC MOD. REMOVE IF UNEXPECTED
			#if cpt == 0:
			#	aj_samples_keep[cpt,:] = np.where(aj_samples_keep[cpt,:] >= 400, aj_samples_keep[cpt,:], np.nan)
			####
			labels_keep.append(labels[i])
			cpt=cpt+1
	
	if show_inc == True:
		inc_stats = inc_stats
		inc_samples = inc_samples
		aj_stats_keep = np.vstack((aj_stats_keep, inc_stats))
		aj_samples_keep = np.vstack((aj_samples_keep, inc_samples[:,0]))
		labels_keep.append("Inclination")
		
	extra_data=[] # Contain all of the pdfs info that are going to be plotted. Minimalistically, extra_data[0]=extra_a2CF
	height_extra=[-1]
	# Add the a2_CF (centrifugal term to a2) for reference
	if show_a2CF == True:
		a2_CF_all, a2_CF_l, a2_CF_mean12, a2_AR_mean=do_a2CF_a2AR(dir_mcmc, step_factor=10, first_sample=0, use_Anl_weight=True)
		a2_CF_stats=make_stats(a2_CF_mean12, confidence=[2.25,16,50,84,97.75])
		if extra_gaussian == []:
			Npdfs=1
			pdf_extra=np.zeros((1, len(a2_CF_mean12)))
			stats_extra=np.zeros((1, len(a2_CF_stats)))
		else:
			Npdfs=1 + len(extra_gaussian) # a2_CF + extra_gaussians
			pdf_extra=np.zeros((Npdfs, len(a2_CF_mean12)))
			stats_extra=np.zeros((Npdfs, len(a2_CF_stats)))		
		# Setting a2_CF
		pos_extra=[1] # a2 position
		pdf_extra[0,:]=a2_CF_mean12
		stats_extra[0,:]=a2_CF_stats
		# Setting extra_pdfs (user-requested, eg. a4 multiple solutions)
		for k in range(len(extra_gaussian)):
			pos_extra.append([extra_gaussian[k][0]]) # Position given by the user
			samples=np.random.normal(extra_gaussian[k][2], extra_gaussian[k][3], len(a2_CF_mean12)) # Generate random gaussian serie of same size as the a2_CF serie
			pdf_extra[k+1,:]=samples
			stats_extra[k+1,:]=[extra_gaussian[k][2] - 2*extra_gaussian[k][3], extra_gaussian[k][2] - extra_gaussian[k][3], extra_gaussian[k][2], extra_gaussian[k][2] + extra_gaussian[k][3], extra_gaussian[k][2] + 2*extra_gaussian[k][3]]
			height_extra.append(extra_gaussian[k][1])
		#
		extra_data.append(pos_extra) # Index for a2
		extra_data.append(pdf_extra)
		extra_data.append(stats_extra)
		extra_data.append(height_extra)

	# -- Perform the plotting --
	nice_hist_matrix(aj_samples_keep, aj_stats_keep, labels_keep, binning=binning, posplots=None, file_out=file_out, extra_data=extra_data)
	#print('Saved files with syntax: ', file_out)
	if show_a2CF == True:
		return labels_txtonly, units_txtonly, aj_stats*unit_nHz, inc_stats, a2_CF_stats  # We return all what was read... because it is supposed to be written later
	else:
		return labels_txtonly, units_txtonly, aj_stats*unit_nHz, inc_stats, []  # We return all what was read... because it is supposed to be written later
	
def set_plots_ajAlm(modeltypes):
	'''
		Set the proper set of aj and activity plots to show
		modeltypes: Defines the shortname to identify the model scenario. Either "1001","1101", "1011", "1111", "1201", "1211"
		keep_aj : A boolean vector that specifies what aj are plotted
	'''
	passed=False
	if str(modeltypes) == "1001":
		keep_aj=[True, False, False, False, False, False]
		passed=True
	if str(modeltypes) == "1011":
		keep_aj=[True, False, True, True, False, False]
		passed=True
	if str(modeltypes) == "1101":
		keep_aj=[True, True, False, True, False, False]
		passed=True
	if str(modeltypes) == "1111":
		keep_aj=[True, True, True, True, False, False]
		passed=True
	if str(modeltypes) == "1201":
		keep_aj=[True, False, False, True, True, True]
		passed=True
	if str(modeltypes) == "1211":
		keep_aj=[True, True, False, True, True, True]
		passed=True
	if passed != True:
		raise ValueError("ERROR : We coud not identify the modeltypes associated to the fit. Set it either to 1001, 1011, 1101, 1111, 1201 or 1211")
	return keep_aj


def set_dir_out(outdir_core, modeltypes, datasource, subtype1="", subtype2=""):
	#
	#	Set the output directory for the post processed data depending on the analysed scenario. 
	#	dir_root: Top directory of all the data
	#	modeltypes: Defines the shortname to identify the model scenario. Either "1101", "1111", "1201", "1211"
	#	subtype1: Defines subtype1 if required. Required for models "1201" or "1211". Eg. Gate or Triangle
	#	subtype2: Defines subtype1 if required. Required for models "1201" or "1211". Eg. decompose_-1
	#
	passed=False
	if str(modeltypes) == "1001":
		outdir=outdir_core + "/" + datasource + "/" + str(modeltypes) + "/" 
		passed=True
	if str(modeltypes) == "1011":
		outdir=outdir_core + "/" + datasource + "/" + str(modeltypes) + "/" 
		passed=True
	if str(modeltypes) == "1101":
		outdir=outdir_core + "/" + datasource + "/" + str(modeltypes) + "/" 
		passed=True
	if str(modeltypes) == "1111":
		outdir=outdir_core + "/" + datasource + "/" + str(modeltypes) + "/" 
		passed=True
	if str(modeltypes) == "1201":
		outdir_core=outdir_core + "/" + datasource + "/" + str(modeltypes) + "/" 
		if subtype1 != "" and subtype2 != "":
			outdir= outdir_core + str(subtype1) + "/" + str(subtype2)+ "/"
		else:
			raise ValueError("Error : You must specify a subtype1 subtype2 when modeltypes is 1201 or 1211")
		passed=True
	if str(modeltypes) == "1211":
		outdir_core=outdir_core + "/" + datasource + "/" + str(modeltypes) + "/" 
		if subtype1 != "" and subtype2 != "":
			outdir= outdir_core + str(subtype1) + "/" + str(subtype2) + "/"
		else:
			raise ValueError("Error : You must specify a subtype1 subtype2 when modeltypes is 1201 or 1211")
		passed=True
	if passed != True:
		raise ValueError("ERROR : We coud not identify the modeltypes associated to the fit. Set it either to 1001, 1101, 1111, 1201 or 1211")
	return outdir

def do_analyse_ajAlm(dir_tamcmc_root, dir_out_root, process_name, data_source, modeltype, subtype1, subtype2, 
			   show_inc=True,
			   cpp_path="/var/services/homes/obenomar/TAMCMC_bin_1.85.1_AMD/", tmpdir="../tmp/"):
	'''
	Main function to perform the plots of aj and/or the Alm terms
	'''
	modeltype_full=[modeltype, subtype1, subtype2]
	#keep_aj, dir_mcmc=set_dir_data_branch(dir_tamcmc_root, modeltypes, data_source, subtype1=subtype1, subtype2=subtype2)
	dir_mcmc=translate_modeltype2dir(dir_tamcmc_root+"/" + data_source, modeltype_full, check_dir_exist=True)
	keep_aj=set_plots_ajAlm(modeltype)
	outdir=set_dir_out(dir_out_root, modeltype, data_source, subtype1=subtype1, subtype2=subtype2)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	outfilename=outdir + process_name + "_" + modeltype + "_" + subtype1 + "_" + subtype2 + ".jpg"
	# Run the correct plotting function in function of the modeltype
	if str(modeltype) == "1001" or str(modeltype) == "1101" or str(modeltype) == "1011" or str(modeltype) == "1111":
		aj_stats, inc_stats,a2CF_stats=show_aj_pdfs(dir_mcmc, keep_aj=keep_aj, file_out=outfilename, binning=40, 
					ignore_slope=True, show_a2CF=True, idlfiles=False, show_inc=show_inc,
					process_name=process_name,
					phase="A", chain=0,
					first_index=0, last_index=-1, period=20,
					cpp_path=cpp_path, outtmpdir=tmpdir)#, extra_gaussian=extra_gaussian)
	else:
		ajAlm_stats, inc_stats=show_ajAlm_pdfs(dir_mcmc, keep_aj=keep_aj, file_out=outfilename, binning=40, 
					ignore_slope=True, idlfiles=False, show_inc=show_inc,
					process_name=process_name,
					phase="A", chain=0,
					first_index=0, last_index=-1, period=20,
					cpp_path=cpp_path, outtmpdir=tmpdir)#, extra_gaussian=extra_gaussian)		

def do_all_analyse_ajAlm(dir_tamcmc_root, dir_out_root, data_source, prefix="kplr",                       
							allowed_modeltype=[1001, 1101,1111,1201,1211],
                            allowed_modelsubtype1=["Gate", "Triangle"],
                            allowed_modelsubtype2=["decompose_-1", "decompose_1", "decompose_2"],
							show_inc=True, show_a2CF=False, 
			   				cpp_path="/var/services/homes/obenomar/TAMCMC_bin_1.85.1_AMD/",
							phase="A", period=20,
							confidence=[2.25,16,50,84,97.75],
							append=False,
							tmpdir="../tmp/"):
	'''
	Calls all possible combination of models that were requested and for all of the existing combinations of data_source and subtypes
	that are requested by the user
	'''   
	# Define allowed types of models by generating combinations of all these type and subtypes
	modeltypes = []
	if any("1001" in str(model) for model in allowed_modeltype):
		modeltypes.append([1001])
	if any("1101" in str(model) for model in allowed_modeltype):
		modeltypes.append([1101])
	if any("1011" in str(model) for model in allowed_modeltype):
		modeltypes.append([1011])
	if any("1111" in str(model) for model in allowed_modeltype):
		modeltypes.append([1111])
	for modeltype in allowed_modeltype:
		for subtype1 in allowed_modelsubtype1:
			for subtype2 in allowed_modelsubtype2:
				if modeltype not in [1001, 1011,1101, 1111]:
					modeltypes.append([modeltype, subtype1, subtype2])

	Nconfidence=len(confidence)
	Nvars=12 # a1, a2, a2CF, a3, a4, a5, a6, epsilon, theta0, delta, inclination
	NcharMax=20*Nvars + 50*3 # Maximum number of characters that are stored within the 3D chararrays
	data_all=[]
	odir_list=[] # To keep record of where this was written
	# Get a list of all of the processes within a given branch of the tree
	for i in range(len(modeltypes)):
		# Identify the type of plot we should do: aj or ajAlm and with what free parameters?
		#print( " modeltypes[i][0] = " , modeltypes[i][0])
		keep_aj=set_plots_ajAlm(modeltypes[i][0])
		# Define the model using the combination of information from modeltypes[i])
		dir_model=translate_modeltype2dir(dir_tamcmc_root +"/" + data_source, modeltypes[i], check_dir_exist=True)
		# Get the names of all of the models that have to be processed within that dir_model
		process_names = [name for name in os.listdir(dir_model) if os.path.isdir(os.path.join(dir_model, name)) and (prefix is None or name.startswith(prefix))]
		# The containers for the data of each models of the confidence array
		Nprocesses=len(process_names)
		data_models=np.chararray((Nconfidence, Nprocesses), itemsize=NcharMax)
		# Create the tree structure for directories, if necessary
		if len(modeltypes[i]) != 1:
			outdir=set_dir_out(dir_out_root, modeltypes[i][0], data_source, subtype1=modeltypes[i][1], subtype2=modeltypes[i][2])
		else:
			outdir=set_dir_out(dir_out_root, modeltypes[i][0], data_source, subtype1="", subtype2="")
		odir_list.append(outdir)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		for j in range(len(process_names)):
			# Define the output filename
			outfile=outdir + "/" + process_names[j]
			outfile_txt=outfile + ".rot"
			# Run the correct plotting function in function of the modeltype
			print(colored("        >>>>> Processing {}...".format(outfile), "green"), flush=True)
			if str(modeltypes[i][0]) == "1000" or str(modeltypes[i][0]) == "1001" or str(modeltypes[i][0]) == "1101" or str(modeltypes[i][0]) == "1111":
				modelcode=str(modeltypes[i][0])
				labels, units, aj_stats, inc_stats,a2CF_stats=show_aj_pdfs(dir_model, keep_aj=keep_aj, file_out=outfile, binning=40, 
							ignore_slope=True, show_a2CF=show_a2CF, idlfiles=False, show_inc=show_inc,
							process_name=process_names[j],
							phase=phase, chain=0,
							first_index=0, last_index=-1, period=period,
							cpp_path=cpp_path, outtmpdir=tmpdir, confidence=confidence)
				write_ajAlm(outfile_txt, labels, units, aj_stats, inc_stats, a2CF_stats, confidence, modeltypes[i][0], stats_type="aj")
				# Prepare everything to be written in a single file that is model-agnostic in Nconfidence files. Files will be taged using "confidence[k]"
				for c in range(len(confidence)):
					d=format_all_ajAlm(labels, units, aj_stats, inc_stats, a2CF_stats, process_names[j], modelcode, "aj", c)
					data_models[c,j]=d
			else:
				# Create a modelcode
				for k in range(len(modeltypes[i])):
					if k == 0:
						modelcode=str(modeltypes[i][k])
					else:
						#modelcode=modelcode+ "_" + str(modeltypes[i][k])
						modelcode=modelcode + str(modeltypes[i][k])
				# Process the MCMC content
				labels, units, ajAlm_stats, inc_stats=show_ajAlm_pdfs(dir_model, keep_aj=keep_aj, file_out=outfile, binning=40, 
							ignore_slope=True, idlfiles=False, show_inc=show_inc,
							process_name=process_names[j],
							phase=phase, chain=0,
							first_index=0, last_index=-1, period=period,
							cpp_path=cpp_path, outtmpdir=tmpdir, confidence=confidence)	
				# Write a summary in the directory of the pdf
				write_ajAlm(outfile_txt, labels, units, ajAlm_stats, inc_stats, [], confidence, modelcode, stats_type="ajAlm")
				# Prepare everything to be written in a single file that is model-agnostic in Nconfidence files. Files will be taged using "confidence[k]"
				for c in range(len(confidence)):
					d=format_all_ajAlm(labels, units, ajAlm_stats, inc_stats, [], process_names[j], modelcode, "ajAlm", c)
					data_models[c,j]=d
		data_all.append(data_models) # because we do not know firsthand the final number of processes within a model, we have to do this list thing
	# Open new files (erasing old ones if any) and setting their headers
	labels_all=["ID", "modelcode", "a1", "a2", "a2CF", "a3", "a4", "a5", "a6", "epsilon", "theta0", "delta", "Inclination"]
	units_all =["N/A", "N/A", "nHz","nHz","nHz","nHz", "nHz","nHz","nHz", "10^-3","deg" , "deg", "deg"]
	for c in range(Nconfidence):
		filename=dir_out_root + "rotation_stats_{:.2f}.summary".format(confidence[c])
		if append == False:
			header="# This file is generated by show_pdf_aj_obs.py::do_all_analyses_ajAlm() version {}\n".format(show_version(verbose=False))
			header=header + "# Condensed statistical summary for all stars and for the rotation-related information\n"
			header=header + "# Confidence = {}%\n".format(confidence[c])
			header=header + "{:<50}".format(labels_all[0])
			header=header + "{:<40}".format(labels_all[1])
			for l in labels_all[2:]:
				header=header + "{:<20}".format(l)
			header=header + "\n"
			header=header + "{:<50}".format(units_all[0])
			header=header + "{:<40}".format(units_all[1])
			for u in units_all[2:]:
				header=header + "{:<20}".format(u)
			header=header + "\n"
			with open(dir_out_root + "rotation_stats_{:.2f}.summary".format(confidence[c]), "w") as f:
				f.write(header)
	# Open the files on which we will write the full summary of the median, etc...
	for c in range(Nconfidence):
		for d in data_all:
			with open(dir_out_root + "rotation_stats_{:.2f}.summary".format(confidence[c]), "a") as f:
				for j in range(len(d[c,:])):
					f.write(d[c,j].decode('utf-8'))
					f.write("\n")
	print(colored("Successfully written all summary files at:","green"), flush=True)
	print(colored("   {}","green").format(dir_out_root), flush=True)
	print(colored("   Individual summary files and plots are at:","green"), flush=True)
	for odir in odir_list:
		print(colored("   {}","green").format(odir), flush=True)
	

def write_ajAlm(outfile, labels, units, aj_stats, inc_stats, a2CF_stats, confidence, modelcode, stats_type="aj"):
	'''
	Write the summary of all of the aj coefficients 
	Only information relevant for the model will be written in a single file per star and model
	See format_all_ajAlm() to see the function that format the full summary in before writting in a single file instead.
	'''
	Naj=len(aj_stats[:,0])
	Nconfidence=len(aj_stats[0,:])
	if stats_type == "aj":
		Ng=Naj+2
		grouped=np.zeros((Ng,Nconfidence))
		labels_all=labels[0:2]
		labels_all.append("a2CF")
		for l in labels[2:]:
			labels_all.append(l)
		units_all=units[0:2]
		units_all.append("nHz")
		for u in units[2:]:
			units_all.append(u)
		grouped[0:2,:]=aj_stats[0:2,:] # a1, a2
		if a2CF_stats != []:
			grouped[2,:]=a2CF_stats
		else:
			grouped[2,:]=np.repeat(0,Nconfidence)
		grouped[3:Naj+1,:]=aj_stats[2:,:] # a3,a4,a5,a6
		grouped[Naj+1,:]=inc_stats
	elif stats_type == "ajAlm":
		labels_all=labels
		units_all=units
		Ng=Naj+1
		grouped=np.zeros((Ng,Nconfidence))
		grouped[0:Naj,:]=aj_stats[:,:] # a1, a3, a5, epsilon ,theta0, delta
		grouped[Naj,:]=inc_stats
	else:
		raise ValueError("stats_type invalid. You must specify either 'aj' or 'ajAlm'")
	header="#ModelCode={}\n".format(modelcode)
	header=header + "#{:<20}".format("Confidence")
	for l in labels_all:
		header=header + "{:<20}".format(l)
	header=header + "{:<20}\n".format("Inclination")
	header=header +"#" + "{:<20}".format("%")
	for u in units_all:
		header=header + "{:<20}".format(u)
	header=header + "{:<20}\n".format("deg")
	data=""
	for c in range(Nconfidence):
		data=data + "{:<20.4f}".format(confidence[c])
		for j in range(Ng):
			data=data + "{:<20.4f}".format(grouped[j,c])
		data=data + "\n"
	f=open(outfile, "w")
	f.write(header)
	f.write(data)
	f.close()
	print("written on :", outfile, flush=True)


def format_all_ajAlm(labels, units, aj_stats, inc_stats, a2CF_stats, starID, modelcode, stats_type, confidence_index):
	'''
	Formating outputs adequately in order to write the full summary of all of the aj coefficients for all of the models 
	All the information will have be written, independently of the model
	data are returned to be written on file later
	'''
	Naj=len(aj_stats[:,0])
	Nconfidence=len(aj_stats[0,:])
	Ng=11 # a1, a2, a2CF, a3, a4, a5, a6, epsilon, theta0, delta, Inclination
	labels_all=["a1", "a2", "a2CF", "a3", "a4", "a5", "a6", "epsilon", "theta0", "delta", "Inclination"]
	grouped=np.zeros((Ng,Nconfidence))
	if stats_type == "aj":
		grouped[0:2,:]=aj_stats[0:2,:] # a1, a2
		if a2CF_stats != []:
			grouped[2,:]=a2CF_stats
		else:
			grouped[2,:]=np.repeat(0,Nconfidence)
		grouped[3:Naj+1,:]=aj_stats[2:,:] # a3,a4,a5,a6
		
	elif stats_type == "ajAlm":
		for i in range(len(labels)):
			if labels[i] in labels_all:
				index = labels_all.index(labels[i])
				grouped[index,:]=aj_stats[i,:]
			else:
				raise ValueError("Could not find the model label within the list of allowed labels. Do you use the correct MCMC model? Debug required.")
	else:
		raise ValueError("stats_type invalid. You must specify either 'aj' or 'ajAlm'")
	grouped[-1,:]=inc_stats # inclination at the end
	data="{:<50}{:<40}".format(starID, modelcode)
	for j in range(Ng):
		data=data + "{:<20.4f}".format(grouped[j,confidence_index])
	data=data + "\n" 
	return data



def test1():
	dir_tamcmc_root='/var/services/homes/obenomar/Temporary/ajAlm_project/data/TRANSFORMED_RUN/Outputs/Legacy/'
	dir_out_root="/var/services/homes/obenomar/Temporary/ajAlm_project/Products/"
	data_source="Legacy"
	allowed_types=[1101,1111,1201,1211]
	allowed_subtypes1=["Gate", "Triangle"]
	allowed_subtypes2=["decompose_-1", "decompose_1", "decompose_2"]
	do_all_analyse_ajAlm(dir_tamcmc_root, dir_out_root, data_source, prefix="kplr",                       
								allowed_modeltype=allowed_types,
								allowed_modelsubtype1=allowed_subtypes1,
								allowed_modelsubtype2=allowed_subtypes2,
								show_inc=True, 
								cpp_path="/var/services/homes/obenomar/TAMCMC_bin_1.85.1_AMD/")


def show_version(verbose=True):
	version="1.31"
	if verbose == True:
		print("show_pdf_aj_obs version {}".format(version), flush=True)
		print("Updated on 1 May 2024: ")
		print("   - Writting the modelCode without underscore in the statistical summary files. This to be consistent with Odds ratio summary file. e.g. 1201_Gate_decompose_-1 -> 1201Gatedecompose-1")	
		print("Updated on 09 Apr 2024: ", flush=True)
		print("   - Adding flush=True for all print to ensure that log files get updated in real time", flush=True)
		print("   - adding user-defined tmpdir for mcmc products using bin2txt etc...")
		print("Updated on 02 Apr 2024: ", flush=True)
		print("   - Adding support for model 1011", flush=True)
		print("Updated on 31 Jan 2024: ", flush=True)
		print("   - Significant refactoring of many functions in order to process large ensemble of Kepler stars", flush=True)
		print("     The main changes are on do_all_analyse_ajAlm() which was extended to create summary files of rotation parameters", flush=True)
		print("Updated on 12 Oct 2023: ", flush=True)
		print("   - Adding the capability of using directly the binary data from TAMCMC instead of the IDL postprocessed data", flush=True)
	return version
show_version()
