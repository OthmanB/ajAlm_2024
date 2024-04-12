from eval_aj_range import read_matrix_tab
from eval_aj_range import read_id_list
from eval_aj_range import match_id_within_lists
import matplotlib.pyplot as plt
import numpy as np
import time
from termcolor import colored
from a2CF_a2AR import do_a2CF_a2AR
from a2CF_a2AR import getdir_idlist
from read_aj_mcmcdata import *
from whisker import Draw_BoxAndWiskers_vertical

def get_pastdata(dir_core='/Users/obenomar/Work/tmp/paper2_asphericity/', IDs_file='/Users/obenomar/Work/tmp/paper2_asphericity/RealData/id_list.txt'):
    '''
    This function regroups the information about Mass, Radius, Teff, a1 and a3 from previous publications for the list of stars
    defined in IDs_file
    '''
    R_sun=6.96342e5 #in km
    M_sun=1.98855e30 #in kg
    V=np.pi*R_sun**3 * 4./3
    #rho_sun=M_sun/V * 1e-12 # Conversion into g/cm-3
    # Reading the IDs within the provided file
    print('Reading the IDs within the provied file: ', IDs_file)
    myID_list=read_id_list(IDs_file)
    # Table of M, R derived from Bellinger+2019 (https://www.aanda.org/articles/aa/pdf/2019/02/aa34461-18.pdf)
    Bellinger_tabA1=dir_core+'/External_data/Legacy_params_Bellinger2019.txt'
    print('Reading Bellinger et al. 2019 table A1...')
    data_bellinger, header_bellinger=read_matrix_tab(Bellinger_tabA1)
    i_M=2
    i_R=3
    # Teff are comming from Pinsonault+2012
    print('Reading Pinsonnault et al. 2012 table 7...')
    Pinsonnault_tab7=dir_core+'/External_data/Pinsonnault2012_table7.dat'
    data_pinsonnault, header_pinsonnault=read_matrix_tab(Pinsonnault_tab7)
    i_Teff=1
    # Benomar et al 2018 for a1 and a3 table
    print('Reading Benomar et al. 2018 summary table of a1 and a3...')
    Benomar_tab=dir_core+'/External_data/Level2-Outputs/a1a3_summary.txt'
    data_benomar, header_benomar=read_matrix_tab(Benomar_tab)
    i_a1=1
    i_a3=4
    # Cross referencing
    print('Attempting to match ID of the published tables with the one provided by the user...')
    print('   -  Verifying that all ID from the user exist...')
    for myID in myID_list:
        id_pos, label=match_id_within_lists(myID, data_bellinger, data_pinsonnault, data_benomar, verbose=False)
    # Gathering data
    print('Gathering data...')
    M=[]
    R=[]
    Teff=[]
    a1_benomar2018=[]
    a3_benomar2018=[]
    for myID in myID_list:
        id_pos, label=match_id_within_lists(myID, data_bellinger, data_pinsonnault, data_benomar, verbose=False)        
        M.append(data_bellinger[id_pos[0], i_M])
        R.append(data_bellinger[id_pos[0],i_R])
        Teff.append(data_pinsonnault[id_pos[1], i_Teff])
        a1_benomar2018.append(data_benomar[id_pos[2], i_a1])
        a3_benomar2018.append(data_benomar[id_pos[2], i_a3])
        #print('Star:', myID, '    Mass: ', M,  '   Radius: ', R,   '    Teff: ', Teff, '     a1: ', a1_benomar2018, '     a3: ', a3_benomar2018)
    return myID_list, M, R, Teff, a1_benomar2018, a3_benomar2018

def propag_err_theta_min_max(theta0, delta, err_theta0, err_delta):
    '''
    Use the error propagation law for gaussian errors in order to compute uncertainties 
    on theta_min and thetea_max coming from (theta0, delta)
    theta0: Value of the mean latitude of activity
    delta: Extension of the active band
    err_theta0: A 2-list with the minus and plus uncertainty on theta0
    err_delta : A 2-list with the minus and plus uncertainty on delta
    '''
    theta_max=theta0 + delta/2
    theta_min=theta0- delta/2
    err_theta_max= [np.sqrt(err_theta0[0]**2 + err_delta[0]**2/4),  np.sqrt(err_theta0[1]**2 + err_delta[1]**2/4)]
    err_theta_min= [np.sqrt(err_theta0[0]**2 + err_delta[0]**2/4),  np.sqrt(err_theta0[1]**2 + err_delta[1]**2/4)]
    return theta_min, theta_max, err_theta_min, err_theta_max

def data_Benomar2022():
    # epsilon, err_epsi_m, err_epsi_p, theta0, theta0_m, theta0_p, delta, delta_m, delta_p, theta_min, err_theta_min_m, err_theta_min_p, theta_max, err_theta_max_m, err_theta_max_p
    # WARNING: ERROR PROPAGATION WAS USED FOR THETA_MIN AND THETA_MAX COMPUTATION
    #          IT IS BETTER LATER TO USE THE FULL PDF TO ACCOUNT FOR CORRELATIONS
    # ACTIVE SUN:
    theta_min, theta_max, err_theta_min, err_theta_max= propag_err_theta_min_max(76, 7, [7,8], [5,16])
    Active_Sun={'Teff':5777, 
                'a1':[421, 10, 10],
                'epsilon':[7.6*1e-4, 5.3*1e-4, 25.5*1e-4], 
                'theta0': [76, 7, 8], 
                'delta': [7, 5, 16], 
                'theta_min':[theta_min, err_theta_min[0], err_theta_min[1]], 
                'theta_max':[theta_max, err_theta_max[0], err_theta_max[1]] 
                }
    # 16 Cyg A:
    theta_min, theta_max, err_theta_min, err_theta_max= propag_err_theta_min_max(71, 4, [14,13], [3,15])
    #
    CygA={"Teff":5825, 
            'a1':[614, 37, 37],
            'epsilon':[5.3*1e-4, 4.1*1e-4, 19.6*1e-4], 
            'theta0': [71, 14, 13], 
            'delta':[4, 3, 15],  
            'theta_min':[theta_min, err_theta_min[0], err_theta_min[1]], 
            'theta_max':[theta_max, err_theta_max[0], err_theta_max[1]]
            }
    # 16 Cyg B (Sol with largest probability and largest uncertainty):
    theta_min, theta_max, err_theta_min, err_theta_max= propag_err_theta_min_max(67, 3, [25,14], [2,13])
    CygB_higha4={'Teff':5750, 
                    'a1':[607, 78, 78],
                    'epsilon':[4.6*1e-4, 3.7*1e-4, 20.8*1e-4], 
                    'theta0': [67, 25, 14], 
                    'delta':[3, 2, 13],  
                    'theta_min':[theta_min, err_theta_min[0], err_theta_min[1]], 
                    'theta_max':[theta_max, err_theta_max[0], err_theta_max[1]] 
                    }
    # 16 Cyg B (Sol with lower probability and lowest uncertainty):
    theta_min, theta_max, err_theta_min, err_theta_max= propag_err_theta_min_max(67, 3, [25,14], [2,13])
    CygB_lowa4={'Teff':0, 
                    'a1':[607, 78, 78],
                    'epsilon':[12.5*1e-4, 7.5*1e-4, 31.9*1e-4], 
                    'theta0': [58, 3, 3], 
                    'delta':[8, 6, 17],  
                    'theta_min':[theta_min, err_theta_min[0], err_theta_min[1]], 
                    'theta_max':[theta_max, err_theta_max[0], err_theta_max[1]] 
                    }
    return Active_Sun, CygA, CygB_lowa4, CygB_higha4

def match_id_newVSpast(IDpastdata, ajdata, ajclassdata, ajdataIDcol=0, ajclassdataIDcol=0, verbose=True):
    '''
        A function that takes the table of parameters of Mass, Radius, Teff, a1_benomar2018, a3_benomar2018 
        and cross reference it with the latest user-provided ID from the new analysis.
        IDpastdata: The table of data  for ID associated to the table with M, R, Teff, a1_benomar2018, a3_benomar2018
        ajdata: The table of data for aj and inclination from the latest analysis
        ajdataIDcol : Column where ID is expected inside ajdata
        ajclassdataIDcol: Column where ID is expected inside the ajclassdata table
        verbose: Show all message (True) or not (False)
        sort: If True, return a sorted pastdata table. If false, return just a pair of indexes in a list + labels that explain the meaning for the indexes
        If it fails, WARNS about the missing stars and continue WITHOUT adding the star to the final list
    '''
    # Number of stars within the new data
    Najdata=len(ajdata[:,0]) 
    i_id=[]
    for i in range(Najdata):
        myID=np.asarray(ajdata[i,ajdataIDcol], dtype=int)
        # Searching a match for the ID between new and old data
        if myID !=0:
            print('Searching ID:', myID)
            pos_pastdata=np.where(myID == np.asarray(IDpastdata, dtype=int))
            pos_classdata=np.where(myID == np.asarray(ajclassdata[:, ajclassdataIDcol], dtype=int))
            err=False
            if len(pos_pastdata[0]) != 0:
                if verbose == True:
                    print('       ID found at (DATA)', pos_pastdata)
                    if len(pos_classdata[0]) !=0:
                        print('       ID found at (CLASS)', pos_classdata)
                        i_id.append([np.squeeze(pos_pastdata), i,np.squeeze(pos_classdata) ])
                    else:
                        i_id.append([np.squeeze(pos_pastdata), i, -1])
                    labels=['M,R,Teff,a1_benomar2018,a3_benomar2018 past data', 'new data for aj and inclination', 'score-based classification (quality) of the new data for aj and inclination']
            else:
                if verbose == True:
                    print(colored('       Error while searching in the past data table: could not find the requested ID', 'red'))
                err=True
    return i_id, labels

def scoreclass2color(score, limits=[-1, 2.5, 3.5]):
    '''
        Receives the score-based table and in function of thresholds return a color:
        The classification is as follow:
            - if score < 0: Unclassified ==> RED
            - if 0 <= score <= 2.5: Bad ==> Light Gray
            - if 2.5 < score <= 3.5: Fair ==> Blue
            - if score > 3.5: Good  ==> Dark Green
        The score is expected to be on a 5-point basis
    '''
    if score <limits[0]:
        color='red'
    if 0 <= score <= limits[1]:
        color='lightgray'
    if limits[1] < score <= limits[2]:
        color='blue'
    if score > limits[2]:
        color='darkgreen'
    return color

def find_indices(lst, element):
    '''
    Returns all element of a list that match a given value or string
    lst: The list within which we search for values
    element: the value or string that is looked inside the list lst
    '''
    result = []
    offset = -1
    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return result
        result.append(offset)


def main_a2CF_a2AR_group(dir_stars='/Users/obenomar/Work/tmp/paper2_asphericity/RealData/products/a1a2a3a4_l0123_a1Gauss/', ids_file='/Users/obenomar/Work/tmp/paper2_asphericity/RealData/id_list.txt', outfile='stars.summary.aj'):
    '''
    Function that measure a2CF and calculate a2AR=a2-a2CF for an ensemble of stars
    Then it writes the outputs into a structure file
    '''
    step_factor=3
    first_sample=0
    use_Anl_weight=True #False
    #
    found_dirs, found_id, missing_id=getdir_idlist(dir_stars, ids_file, ids_col=0)
    if len(missing_id) !=0:
        print('Warning: Some stars were not found in the root directory:')
        for m in missing_id:
            print(m)
        time.sleep(2)
    #
    i=0
    outputs='#{0:10} {1:10} {2:10} {3:10}'.format('id' ,  'a1'  ,  'err_a1_m',  'err_a1_p')
    outputs=outputs +'{0:10} {1:10} {2:10} {3:10} {4:10} {5:10} {6:10} {7:10} {8:10}'.format('a2','err_a2_m','err_a2_p','a2_CF','err_a2_CF_m','err_a2_CF_p','a2_AR','err_AR_m','err_AR_p')
    outputs=outputs + '{0:10} {1:10} {2:10}'.format('a3','err_a3_m', 'err_a3_p')
    outputs=outputs + '{0:10} {1:10} {2:10}'.format('a4','err_a4_m', 'err_a4_p')
    outputs=outputs + '{0:10} {1:10} {2:10}\n'.format('inc','err_inc_m', 'err_inc_p')
    for d in found_dirs:
        a2_CF_all, a2_CF_l, a2_CF_mean12_samples, a2_AR_mean_samples=do_a2CF_a2AR(d, step_factor=step_factor, first_sample=first_sample, use_Anl_weight=use_Anl_weight)
        a2_CF_stats=make_stats(a2_CF_mean12_samples, confidence=[2.25,16,50,84,97.75])
        a2_AR_stats=make_stats(a2_AR_mean_samples, confidence=[2.25,16,50,84,97.75])
        a2_CF_err=make_error_from_stats(a2_CF_stats)
        a2_AR_err=make_error_from_stats(a2_AR_stats)
        aj_stats, inc_stats, aj_samples, inc_samples, Nsamples=get_aj_inc_star(d)
        aj_err=make_error_from_stats(aj_stats)
        inc_err=make_error_from_stats(inc_stats)
        #    id   a1    err_a1_m  err_a1_p
        outputs=outputs + '{0:10} {1:10.4f} {2:10.4f} {3:10.4f}'.format(found_id[i], aj_stats[0,2], aj_err[0,0], aj_err[1,0])
        #    a2  err_a2_m  err_a2_p     a2_CF   err_a2_CF_m   err_a2_CF_p    a2_AR  err_AR_m err_AR_p 
        outputs=outputs +  '{0:10.4f} {1:10.4f} {2:10.4f} {3:10.4f} {4:10.4f} {5:10.4f} {6:10.4f} {7:10.4f} {8:10.4f}'.format(aj_stats[1,2], aj_err[0,1], aj_err[1,1], a2_CF_stats[2], a2_CF_err[0], a2_CF_err[1], a2_AR_stats[2], a2_AR_err[0], a2_AR_err[1])
        #    a3  err_a3_m  err_a3_p 
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}'.format(aj_stats[2,2], aj_err[0,2], aj_err[1,2])
        #    a4  err_a4_m  err_a4_p
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}'.format(aj_stats[3,2], aj_err[0,3], aj_err[1,3])
        #    inc  err_inc_m  err_inc_p
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}\n'.format(inc_stats[2], inc_err[0], inc_err[1])
        i=i+1
    #
    f=open(outfile, 'w')
    f.write(outputs)
    f.close() 
    print("All data are saved in:", outfile)

def main_ajAlm_group(dir_stars='/Users/obenomar/Work/tmp/paper2_asphericity/RealData/products/ajAlm/', ids_file='/Users/obenomar/Work/tmp/paper2_asphericity/RealData/id_list.txt', outfile='stars.summary.ajAlm'):
    '''
    Function that measure a2CF and calculate a2AR=a2-a2CF for an ensemble of stars
    Then it writes the outputs into a structure file
    '''
    #
    found_dirs, found_id, missing_id=getdir_idlist(dir_stars, ids_file, ids_col=0)
    if len(missing_id) !=0:
        print('Warning: Some stars were not found in the root directory:')
        for m in missing_id:
            print(m)
        time.sleep(2)
    #
    i=0
    outputs='#{0:10} {1:10} {2:10} {3:10}'.format('id' ,  'a1'  ,  'e_a1_m',  'e_a1_p')
    outputs=outputs + '{0:10} {1:10} {2:10}'.format('a3','e_a3_m', 'e_a3_p')
    outputs=outputs + '{0:10} {1:12} {2:12}'.format('epsi','e_epsi_m', 'e_epsi_p')
    outputs=outputs + '{0:10} {1:10} {2:10}'.format('theta','e_theta_m', 'e_theta_p')
    outputs=outputs + '{0:10} {1:10} {2:10}'.format('delta','e_delta_m', 'e_delta_p')
    outputs=outputs + '{0:10} {1:10} {2:10}'.format('inc','e_inc_m', 'e_inc_p')
    outputs=outputs + '{0:12} {1:15} {2:15}'.format('theta_min','e_theta_min_m', 'e_theta_min_p')
    outputs=outputs + '{0:12} {1:15} {2:15}\n'.format('theta_max','e_theta_max_m', 'e_theta_max_p')
    param_file='plength.txt'
    for d in found_dirs:
        # We first need to know i0_ajAlm:
        plength=read_parameters_length(d + '/Files/', param_file) # Read the parameters_length file and retrieves plength
        i0_ajAlm=sum(plength[0:6]) # First Position after Nf_el list
        i0_inc=sum(plength[0:-2]) # The last parameter is before the extra parameters that are at the end ([-1] position)
        # Then we read the aj and Alm parameters 
        aj_params, epsi_params, theta0_params, delta_params=read_rot_ajAlm(d+'/Files/', i0_ajAlm) # Extract all of the samples for a1,a3, epsilon, theta0, delta in a ajAlm model
        # Compose theta_min and theta_max from theta0 and delta (assuming a window-function for the active zone)
        theta_min=theta0_params - delta_params/2
        theta_max=theta0_params + delta_params/2
        # Read the inclination
        inc_samples=read_inc(d+'/Files/', i0_inc)
        # Then we need to organise them in simple statistics
        a1_stats=make_stats(aj_params[0,:], confidence=[2.25,16,50,84,97.75])
        a3_stats=make_stats(aj_params[1,:], confidence=[2.25,16,50,84,97.75])
        #a5_stats=make_stats(ajAlm_params[2,:], confidence=[2.25,16,50,84,97.75])
        epsilon_stats=make_stats(epsi_params, confidence=[2.25,16,50,84,97.75])
        theta0_stats=make_stats(theta0_params, confidence=[2.25,16,50,84,97.75])
        delta_stats=make_stats(delta_params, confidence=[2.25,16,50,84,97.75])
        theta_min_stats=make_stats(theta_min, confidence=[2.25,16,50,84,97.75])
        theta_max_stats=make_stats(theta_max, confidence=[2.25,16,50,84,97.75])
        inc_stats=make_stats(inc_samples)
        a1_err=make_error_from_stats(a1_stats)
        a3_err=make_error_from_stats(a3_stats)
        epsilon_err=make_error_from_stats(epsilon_stats)
        theta0_err=make_error_from_stats(theta0_stats)
        theta_min_err=make_error_from_stats(theta_min_stats)
        theta_max_err=make_error_from_stats(theta_max_stats)
        delta_err=make_error_from_stats(delta_stats)
        inc_err=make_error_from_stats(inc_stats)
        #    id   a1    err_a1_m  err_a1_p
        outputs=outputs + '{0:10} {1:10.4f} {2:10.4f} {3:10.4f}'.format(found_id[i], a1_stats[2], a1_err[0], a1_err[1])
        #    a3  err_a3_m  err_a3_p 
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}'.format(a3_stats[2], a3_err[0], a3_err[1])
        #    epsilon  err_epsi_m  err_epsi_p
        outputs=outputs + '{0:12.6f} {1:12.6f} {2:12.6f}'.format(epsilon_stats[2], epsilon_err[0], epsilon_err[1])
        #    theta0  err_theta0_m  err_theta0_p
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}'.format(theta0_stats[2], theta0_err[0], theta0_err[1])
        #    delta  err_delta_m  err_delta_p
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}'.format(delta_stats[2], delta_err[0], delta_err[1])
        #    inc  err_inc_m  err_inc_p
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}'.format(inc_stats[2], inc_err[0], inc_err[1])
        #    theta_min  err_theta_min_m  err_theta_min_p
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}'.format(theta_min_stats[2], theta_min_err[0], theta_min_err[1])
        #    theta_max  err_theta_max_m  err_theta_max_p
        outputs=outputs + '{0:10.4f} {1:10.4f} {2:10.4f}\n'.format(theta_max_stats[2], theta_max_err[0], theta_max_err[1])
        i=i+1
    #
    f=open(outfile, 'w')
    f.write(outputs)
    f.close() 
    print("All data are saved in:", outfile)

def main_a2CF_a2AR():
    '''
    Function to quickly call  main_ajAlm_group in a proper manner
    '''
    #dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/19992002_incfix_fast_Priorevalrange/'
    #dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/20062009_incfix_fast/'
    #dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/kplr012069424_kasoc-wpsd_slc_v1_a2a3a4_nol3/'
    #dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/kplr012069449_kasoc-wpsd_slc_v1_a2a3a4_nol3/'
    dir_core='/Users/obenomar/tmp/test_a2AR/tmp/Realdata/products/kplr012069449_kasoc-wpsd_slc_v1_a2a3a4_nol3_GU2/'
    step_factor=3
    first_sample=0
    use_Anl_weight=True #False
    a2_CF_all, a2_CF_l, a2_CF_mean12, a2_AR_mean=do_a2CF_a2AR(dir_core, step_factor=step_factor, first_sample=first_sample, use_Anl_weight=use_Anl_weight)
    #
    exit()


def show_aj(dir_core, IDs_file, aj_summary_file, aj_classification_file, activity_summary_file=None, dir_out=''):
    '''
    Function that plots the values of a2 and a4 in function of
        - new a1
        - new a3
        - new inclination
        - previous a1
        - previous a3
        - Teff
        - M
        - R
    If an activity_summary_file (with epsilon, theta0, delta) is also provided, shows its quantities the same quantities listed above
    '''

    print('Reading latest aj analysis summary from ', aj_summary_file)
    data_ajnew, header_ajnew=read_matrix_tab(aj_summary_file, delete0=True)
    i_a1=1 # err_m and err_p follow
    i_a2=4 # err_m and err_p follow
    i_a2CF=7 # err_m and err_p follow
    i_a2AR=10 # err_m and err_p follow
    i_a3=13 # err_m and err_p follow
    i_a4=16 # err_m and err_p follow
    i_inc=19 # err_m and err_p follow
    if activity_summary_file !=None:
        print('Reading activity analysis summary from ', activity_summary_file)
        data_ajAlm, header_ajnew=read_matrix_tab(activity_summary_file, delete0=True)
        i_a1_act=1 # err_m and err_p follow
        i_a3_act=4 # err_m and err_p follow
        i_epsi=7 # err_m and err_p follow
        i_theta0=10 # err_m and err_p follow
        i_delta=13 # err_m and err_p follow
        i_inc_act=16 # err_m and err_p follow
        i_theta_min=19 
        i_theta_max=22 
    else:
        print('No activity summary file provided... skipping some plots')
    #
    print('Reading the score-based classification file...')
    data_aj_class, header_class=read_matrix_tab(aj_classification_file, delete0=True)
    print('Get past published data...')
    myID_list, M, R, Teff, a1_benomar2018, a3_benomar2018=get_pastdata(dir_core, IDs_file)
    # Getting the Sun and 16CygA/B values from past data
    Active_Sun, CygA, CygB_lowa4, CygB_higha4=data_Benomar2022()
    print('Verifying that all new data has corresponding past/published data and cross match the index between tables...')
    # Verifying that all new data has corresponding past/published data and cross match the index between tables
    indexes, labels=match_id_newVSpast(myID_list, data_ajnew, data_aj_class, verbose=True)
    indexes_ajAlmVSpast, labels=match_id_newVSpast(myID_list, data_ajAlm, data_aj_class, verbose=True)
    # Sorting the pastdata to be sure that it matches the data_ajnew
    NajNew_lines=len(data_ajnew[:,0])
    NajNew_cols=len(data_ajnew[0,:])
    M_sorted=[]
    R_sorted=[]
    Teff_sorted=[]
    a1past_sorted=[]
    a3past_sorted=[] #np.zeros((NajNew_lines, NajNew_cols))
    colors=[]
    M_sorted_Alm=[]
    R_sorted_Alm=[]
    Teff_sorted_Alm=[]
    a1past_sorted_Alm=[]
    a3past_sorted_Alm=[] #np.zeros((NajNew_lines, NajNew_cols))
    colors_Alm=[]
    #print(indexes)
    for ind in indexes:
        M_sorted.append(M[ind[0]][0])
        R_sorted.append(R[ind[0]][0])
        Teff_sorted.append(Teff[ind[0]][0])
        a1past_sorted.append(a1_benomar2018[ind[0]][0])
        a3past_sorted.append(a3_benomar2018[ind[0]][0])
        colors.append(scoreclass2color(data_aj_class[ind[2],1]))
        #print(data_ajnew[ind[1],0], data_aj_class[ind[2],0],  data_aj_class[ind[2],1], scoreclass2color(data_aj_class[ind[2],1]))
    M_sorted=np.array(M_sorted)
    R_sorted=np.array(R_sorted)
    Teff_sorted=np.array(Teff_sorted)
    a1past_sorted=np.array(a1past_sorted)
    a3past_sorted=np.array(a3past_sorted)
    #
    cols=['red', 'lightgray', 'blue', 'darkgreen']
    pos_unclassified=find_indices(colors, cols[0])
    pos_bad=find_indices(colors, cols[1])
    pos_fair=find_indices(colors, cols[2])
    pos_good=find_indices(colors, cols[3])
    #
    if activity_summary_file !=None:
        for ind in indexes_ajAlmVSpast:
            M_sorted_Alm.append(M[ind[0]][0])
            R_sorted_Alm.append(R[ind[0]][0])
            Teff_sorted_Alm.append(Teff[ind[0]][0])
            a1past_sorted_Alm.append(a1_benomar2018[ind[0]][0])
            a3past_sorted_Alm.append(a3_benomar2018[ind[0]][0])
            colors_Alm.append(scoreclass2color(data_aj_class[ind[2],1]))
    M_sorted_Alm=np.array(M_sorted_Alm)
    R_sorted_Alm=np.array(R_sorted_Alm)
    Teff_sorted_Alm=np.array(Teff_sorted_Alm)
    a1past_sorted_Alm=np.array(a1past_sorted_Alm)
    a3past_sorted_Alm=np.array(a3past_sorted_Alm)
    #
    pos_unclassified_Alm=find_indices(colors_Alm, cols[0])
    pos_bad_Alm=find_indices(colors_Alm, cols[1])
    pos_fair_Alm=find_indices(colors_Alm, cols[2])
    pos_good_Alm=find_indices(colors_Alm, cols[3])
    #
    # a2AR and a4 function of M
    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    #ax[0].errorbar(M_sorted[pos_bad], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
    ax[0].errorbar(M_sorted[pos_fair], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[0].errorbar(M_sorted[pos_good], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[0].axhline(0, linestyle='--', color='black')
    ax[0].set_ylabel('a2AR (nHz)')
    #ax[1].errorbar(M_sorted[pos_bad], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
    ax[1].errorbar(M_sorted[pos_fair], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[1].errorbar(M_sorted[pos_good], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[1].set_ylabel('a4 (nHz)')
    ax[1].axhline(0, linestyle='--', color='black')
    ax[1].set_xlabel('Mass')
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper left')	
    fig_1d.savefig(dir_out+ 'summary_a2a4_Mass.jpg', dpi=300)
    #
    # a2AR and a4 function of R
    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    #ax[0].errorbar(R_sorted[pos_bad], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
    ax[0].errorbar(R_sorted[pos_fair], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[0].errorbar(R_sorted[pos_good], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[0].axhline(0, linestyle='--', color='black')
    ax[0].set_ylabel('a2 (nHz)')
    #ax[1].errorbar(R_sorted[pos_bad], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
    ax[1].errorbar(R_sorted[pos_fair], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[1].errorbar(R_sorted[pos_good], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[1].axhline(0, linestyle='--', color='black')
    ax[1].set_ylabel('a4 (nHz)')
    ax[1].set_xlabel('Radius')
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper left')	
    fig_1d.savefig(dir_out+ 'summary_a2a4_Radius.jpg', dpi=300)
    #
    # a2AR and a4 function of Teff
    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    #ax[0].errorbar(Teff_sorted[pos_bad], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
    ax[0].errorbar(Teff_sorted[pos_fair], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[0].errorbar(Teff_sorted[pos_good], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[0].axhline(0, linestyle='--', color='black')
    ax[0].set_ylabel('a2AR (nHz)')
    #ax[1].errorbar(Teff_sorted[pos_bad], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
    ax[1].errorbar(Teff_sorted[pos_fair], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[1].errorbar(Teff_sorted[pos_good], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[1].axhline(0, linestyle='--', color='black')
    ax[1].set_ylabel('a4 (nHz)')
    ax[1].set_xlabel('Teff (K)')
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper left')	
    fig_1d.savefig(dir_out+ 'summary_a2a4_Teff.jpg', dpi=300)
    #
    # a2AR and a4 function of past a1
    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    #ax[0].errorbar(a1past_sorted[pos_bad], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
    ax[0].errorbar(a1past_sorted[pos_fair], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[0].errorbar(a1past_sorted[pos_good], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[0].axhline(0, linestyle='--', color='black')
    ax[0].set_ylabel('a2AR (nHz)')
    #ax[1].errorbar(a1past_sorted[pos_bad], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
    ax[1].errorbar(a1past_sorted[pos_fair], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[1].errorbar(a1past_sorted[pos_good], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[1].axhline(0, linestyle='--', color='black')
    ax[1].set_ylabel('a4 (nHz)')
    ax[1].set_xlabel('Past a1 (nHz)')
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper left')	
    fig_1d.savefig(dir_out+ 'summary_a2a4_pasta1.jpg', dpi=300)
   #
    # a2AR and a4 function of past a3
    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    #ax[0].errorbar(a3past_sorted[pos_bad], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
    ax[0].errorbar(a3past_sorted[pos_fair], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[0].errorbar(a3past_sorted[pos_good], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[0].axhline(0, linestyle='--', color='black')
    ax[0].set_ylabel('a2AR (nHz)')
    #ax[1].errorbar(a3past_sorted[pos_bad], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
    ax[1].errorbar(a3past_sorted[pos_fair], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[1].errorbar(a3past_sorted[pos_good], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[1].axhline(0, linestyle='--', color='black')
    ax[1].set_ylabel('a4 (nHz)')
    ax[1].set_xlabel('Past a3 (nHz)')
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper left')	
    fig_1d.savefig(dir_out+ 'summary_a2a4_pasta3.jpg', dpi=300)

    # a2AR and a4 function of new a1
    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    #ax[0].errorbar(data_ajnew[pos_bad, i_a1], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
    ax[0].errorbar(data_ajnew[pos_fair, i_a1], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[0].errorbar(data_ajnew[pos_good, i_a1], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[0].axhline(0, linestyle='--', color='black')
    ax[0].set_ylabel('a2AR (nHz)')
    #ax[1].errorbar(data_ajnew[pos_bad, i_a1], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
    ax[1].errorbar(data_ajnew[pos_fair, i_a1], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[1].errorbar(data_ajnew[pos_good, i_a1], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[1].axhline(0, linestyle='--', color='black')
    ax[1].set_ylabel('a4 (nHz)')
    ax[1].set_xlabel('New a1')
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper left')	
    fig_1d.savefig(dir_out+ 'summary_a2a4_newa1.jpg', dpi=300)
   #
    # a2AR and a4 function of new a3
    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    #ax[0].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
    ax[0].errorbar(data_ajnew[pos_fair, i_a3], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[0].errorbar(data_ajnew[pos_good, i_a3], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[0].axhline(0, linestyle='--', color='black')
    ax[0].set_ylabel('a2AR (nHz)')
    #ax[1].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
    ax[1].errorbar(data_ajnew[pos_fair, i_a3], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[1].errorbar(data_ajnew[pos_good, i_a3], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[1].axhline(0, linestyle='--', color='black')
    ax[1].set_ylabel('a4 (nHz)')
    ax[1].set_xlabel('New a3')
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper left')	
    fig_1d.savefig(dir_out+ 'summary_a2a4_newa3.jpg', dpi=300)

    # a2AR and a4 function of inclination
    fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
    #ax[0].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
    ax[0].errorbar(data_ajnew[pos_fair, i_inc], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[0].errorbar(data_ajnew[pos_good, i_inc], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[0].axhline(0, linestyle='--', color='black')
    ax[0].set_ylabel('a2AR (nHz)')
    #ax[1].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
    ax[1].errorbar(data_ajnew[pos_fair, i_inc], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
    ax[1].errorbar(data_ajnew[pos_good, i_inc], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
    ax[1].axhline(0, linestyle='--', color='black')
    ax[1].set_ylabel('a4 (nHz)')
    ax[1].set_xlabel('Inclination')
    # Handling legends
    ax[0].legend(fontsize=10, loc='upper left')	
    fig_1d.savefig(dir_out+ 'summary_a2a4_inc.jpg', dpi=300)

    ###############################################################
    ###############################################################
    ########################### ajAlm #############################
    ###############################################################
    ###############################################################
    if activity_summary_file != None:
        # a2AR and a4 function of theta0 FROM ajAlm
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        #ax[0].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
        ax[0].errorbar(data_ajAlm[pos_fair, i_theta0], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(data_ajAlm[pos_good, i_theta0], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('a2AR (nHz)')
        ax[0].set_xlim(0, 90)
        #ax[1].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
        ax[1].errorbar(data_ajAlm[pos_fair, i_theta0], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(data_ajAlm[pos_good, i_theta0], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].axhline(0, linestyle='--', color='black')
        ax[1].set_ylabel('a4 (nHz)')
        ax[1].set_xlabel('theta0')
        ax[1].set_xlim(0, 90)
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_a2a4_theta0.jpg', dpi=300)

        # a2AR and a4 function of delta FROM ajAlm
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        #ax[0].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a2AR], yerr=[data_ajnew[pos_bad, i_a2AR+1],data_ajnew[pos_bad, i_a2AR+2]], marker='o',linestyle='', color=cols[1])
        ax[0].errorbar(data_ajAlm[pos_fair, i_delta], data_ajnew[pos_fair, i_a2AR], yerr=np.squeeze([data_ajnew[pos_fair, i_a2AR+1],data_ajnew[pos_fair, i_a2AR+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(data_ajAlm[pos_good, i_delta], data_ajnew[pos_good, i_a2AR], yerr=np.squeeze([data_ajnew[pos_good, i_a2AR+1],data_ajnew[pos_good, i_a2AR+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('a2AR (nHz)')
        ax[0].set_xlim(0, 45)
        #ax[1].errorbar(data_ajnew[pos_bad, i_a3], data_ajnew[pos_bad, i_a4], yerr=np.squeeze([data_ajnew[pos_bad, i_a4+1],data_ajnew[pos_bad, i_a4+2]]), marker='o',linestyle='', color=cols[1])
        ax[1].errorbar(data_ajAlm[pos_fair, i_delta], data_ajnew[pos_fair, i_a4], yerr=np.squeeze([data_ajnew[pos_fair, i_a4+1],data_ajnew[pos_fair, i_a4+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(data_ajAlm[pos_good, i_delta], data_ajnew[pos_good, i_a4], yerr=np.squeeze([data_ajnew[pos_good, i_a4+1],data_ajnew[pos_good, i_a4+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].axhline(0, linestyle='--', color='black')
        ax[1].set_ylabel('a4 (nHz)')
        ax[1].set_xlabel('delta')
        ax[1].set_xlim(0, 45)
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_a2a4_delta.jpg', dpi=300)

        # theta0 and delta and epsilon function of M
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(M_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(M_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[0].set_yscale('log')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[1].errorbar(M_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(M_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(M_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(M_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[2].set_xlabel('Mass')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_Mass.jpg', dpi=300)

        # theta0 and delta and epsilon function of R
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(R_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(R_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[0].set_yscale('log')
        ax[1].errorbar(R_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(R_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(R_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(R_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[2].set_xlabel('Radius')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_Radius.jpg', dpi=300)

        # theta0 and delta and epsilon function of Teff
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[0].set_yscale('log')
        ax[1].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[2].set_xlabel('Teff (K)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_Teff.jpg', dpi=300)

        # theta0 and delta and epsilon function of past a1
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[0].set_yscale('log')
        ax[1].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(45, linestyle='--', color='black')
        ax[2].set_xlabel('Past a1 (nHz)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_pasta1.jpg', dpi=300)

        # theta0 and delta and epsilon function of past a3
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(a3past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(a3past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel('epsilon')
        ax[0].set_yscale('log')
        ax[1].errorbar(a3past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta0+1],data_ajAlm[pos_fair_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(a3past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta0], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta0+1],data_ajAlm[pos_good_Alm, i_theta0+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel('theta0')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(90, linestyle='--', color='black')
        ax[2].errorbar(a3past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_delta+1],data_ajAlm[pos_fair_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(a3past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_delta], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_delta+1],data_ajAlm[pos_good_Alm, i_delta+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel('delta')
        ax[2].set_ylim(0, 45)
        #ax[2].axhline(45, linestyle='--', color='black')
        ax[2].set_xlabel('Past a3 (nHz)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_pasta3.jpg', dpi=300)

        # ------
        # Added on 14 Dec 2022
        # ------
        # theta_min and theta_max function of Teff
        fig_1d, ax = plt.subplots(3,1, figsize=(12, 6))
        ax[0].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
        ax[1].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta_min], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta_min+1],data_ajAlm[pos_fair_Alm, i_theta_min+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[1].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta_min], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta_min+1],data_ajAlm[pos_good_Alm, i_theta_min+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[1].set_ylabel(r'$\theta_{min}$')
        ax[1].set_ylim(0, 90)
        #ax[1].axhline(0, linestyle='--', color='black')
        ax[2].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_theta_max], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_theta_max+1],data_ajAlm[pos_fair_Alm, i_theta_max+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[2].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_theta_max], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_theta_max+1],data_ajAlm[pos_good_Alm, i_theta_max+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        ax[2].set_ylabel(r'$\theta_{max}$')
        ax[2].set_ylim(0, 90)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[2].set_xlabel('Teff (K)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetamin_thetamax_Teff.jpg', dpi=300)
        # theta_min , theta_max, theta0 in a single plot ONLY FAIR AND GOOD STARS
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        ax[0].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        # The Sun, 16 CygA, 16 Cyg B 
        extra_a1=[Active_Sun['Teff'], CygA['Teff'], CygB_higha4['Teff']]
        extra_epsi=[Active_Sun['epsilon'][0], CygA['epsilon'][0] , CygB_higha4['epsilon'][0]]
        extra_epsi_err=[Active_Sun['epsilon'][1], CygA['epsilon'][1] , CygB_higha4['epsilon'][1]]
        extra_epsi_err=[Active_Sun['epsilon'][2], CygA['epsilon'][2] , CygB_higha4['epsilon'][2]]
        ax[0].errorbar(extra_a1, extra_epsi, yerr=extra_epsi_err, marker='o',linestyle='', color=cols[2])
#
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
#        for p in pos_bad_Alm:
#            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+2] ]
#            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['lightgray', 'lightgray', 'gray'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_fair_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['blue', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_good_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # Add the Sun, 16 Cyg A and B from Benomar 2022
        # The Active Sun
        stats=[Active_Sun['theta_min'][0] - Active_Sun['theta_min'][1] , Active_Sun['theta_min'][0], Active_Sun['theta0'][0], Active_Sun['theta_max'][0], Active_Sun['theta_max'][0] + Active_Sun['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(Active_Sun['Teff'], stats, ax=ax[1], width=15, color=['gold', 'gold', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg A
        stats=[CygA['theta_min'][0] - CygA['theta_min'][1] , CygA['theta_min'][0], CygA['theta0'][0], CygA['theta_max'][0], CygA['theta_max'][0] + CygA['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygA['Teff'], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg B HIGH A4 (HIGHEST UNCERTAINTY)
        stats=[CygB_higha4['theta_min'][0] - CygB_higha4['theta_min'][1] , CygB_higha4['theta_min'][0], CygB_higha4['theta0'][0], CygB_higha4['theta_max'][0], CygB_higha4['theta_max'][0] + CygB_higha4['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygB_higha4['Teff'], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        #ax[1].set_xlim(5200, 7200)
        ax[1].set_ylabel(r'Activity band $\Delta\theta$')
        ax[1].set_ylim(0, 91)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[1].set_xlabel('Teff (K)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetaBand_Teff_GoodFair.jpg', dpi=300)
        #
        # theta_min , theta_max, theta0 in a single plot WITH gray BAD STARS
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        ax[0].errorbar(Teff_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(Teff_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        # The Sun, 16 CygA, 16 Cyg B 
        extra_a1=[Active_Sun['Teff'], CygA['Teff'], CygB_higha4['Teff']]
        extra_epsi=[Active_Sun['epsilon'][0], CygA['epsilon'][0] , CygB_higha4['epsilon'][0]]
        extra_epsi_err=[Active_Sun['epsilon'][1], CygA['epsilon'][1] , CygB_higha4['epsilon'][1]]
        extra_epsi_err=[Active_Sun['epsilon'][2], CygA['epsilon'][2] , CygB_higha4['epsilon'][2]]
        ax[0].errorbar(extra_a1, extra_epsi, yerr=extra_epsi_err, marker='o',linestyle='', color=cols[2])
#
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
        for p in pos_bad_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['lightgray', 'lightgray', 'gray'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_fair_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['blue', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_good_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # The Active Sun
        stats=[Active_Sun['theta_min'][0] - Active_Sun['theta_min'][1] , Active_Sun['theta_min'][0], Active_Sun['theta0'][0], Active_Sun['theta_max'][0], Active_Sun['theta_max'][0] + Active_Sun['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(Active_Sun['Teff'], stats, ax=ax[1], width=15, color=['gold', 'gold', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg A
        stats=[CygA['theta_min'][0] - CygA['theta_min'][1] , CygA['theta_min'][0], CygA['theta0'][0], CygA['theta_max'][0], CygA['theta_max'][0] + CygA['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygA['Teff'], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg B HIGH A4 (HIGHEST UNCERTAINTY)
        stats=[CygB_higha4['theta_min'][0] - CygB_higha4['theta_min'][1] , CygB_higha4['theta_min'][0], CygB_higha4['theta0'][0], CygB_higha4['theta_max'][0], CygB_higha4['theta_max'][0] + CygB_higha4['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygB_higha4['Teff'], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        #ax[1].set_xlim(5200, 7200)
        ax[1].set_ylabel(r'Activity band $\Delta\theta$')
        ax[1].set_ylim(0, 91)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[1].set_xlabel('Teff (K)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetaBand_Teff_GoodFairBad.jpg', dpi=300)

        # ---------
        # ---------

        # theta_min , theta_max, theta0 in a single plot ONLY FAIR AND GOOD STARS
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        ax[0].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        # The Sun, 16 CygA, 16 Cyg B 
        extra_a1=[Active_Sun['a1'][0], CygA['a1'][0], CygB_higha4['a1'][0]]
        extra_epsi=[Active_Sun['epsilon'][0], CygA['epsilon'][0] , CygB_higha4['epsilon'][0]]
        extra_epsi_err=[Active_Sun['epsilon'][1], CygA['epsilon'][1] , CygB_higha4['epsilon'][1]]
        extra_epsi_err=[Active_Sun['epsilon'][2], CygA['epsilon'][2] , CygB_higha4['epsilon'][2]]
        ax[0].errorbar(extra_a1, extra_epsi, yerr=extra_epsi_err, marker='o',linestyle='', color=cols[2])

        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
#        for p in pos_bad_Alm:
#            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+2] ]
#            Draw_BoxAndWiskers_vertical(Teff_sorted_Alm[p], stats, ax=ax[1], width=15, color=['lightgray', 'lightgray', 'gray'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_fair_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['blue', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_good_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # Add the Sun, 16 Cyg A and B from Benomar 2022
        # The Active Sun
        stats=[Active_Sun['theta_min'][0] - Active_Sun['theta_min'][1] , Active_Sun['theta_min'][0], Active_Sun['theta0'][0], Active_Sun['theta_max'][0], Active_Sun['theta_max'][0] + Active_Sun['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(Active_Sun['a1'][0], stats, ax=ax[1], width=15, color=['gold', 'gold', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg A
        stats=[CygA['theta_min'][0] - CygA['theta_min'][1] , CygA['theta_min'][0], CygA['theta0'][0], CygA['theta_max'][0], CygA['theta_max'][0] + CygA['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygA['a1'][0], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg B HIGH A4 (HIGHEST UNCERTAINTY)
        stats=[CygB_higha4['theta_min'][0] - CygB_higha4['theta_min'][1] , CygB_higha4['theta_min'][0], CygB_higha4['theta0'][0], CygB_higha4['theta_max'][0], CygB_higha4['theta_max'][0] + CygB_higha4['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygB_higha4['a1'][0], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        #ax[1].set_xlim(5200, 7200)
        ax[1].set_ylabel(r'Activity band $\Delta\theta$')
        ax[1].set_ylim(0, 91)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[1].set_xlabel(r'$a_1$ (nHz)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetaBand_a1_GoodFair.jpg', dpi=300)
        #
        # theta_min , theta_max, theta0 in a single plot WITH gray BAD STARS
        fig_1d, ax = plt.subplots(2,1, figsize=(12, 6))
        ax[0].errorbar(a1past_sorted_Alm[pos_fair_Alm], data_ajAlm[pos_fair_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_fair_Alm, i_epsi+1],data_ajAlm[pos_fair_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[2], label='Fair (2.5<s<3.5)')
        ax[0].errorbar(a1past_sorted_Alm[pos_good_Alm], data_ajAlm[pos_good_Alm, i_epsi], yerr=np.squeeze([data_ajAlm[pos_good_Alm, i_epsi+1],data_ajAlm[pos_good_Alm, i_epsi+2]]), marker='o',linestyle='', color=cols[3], label='Good (s>3.5)')
        # The Sun, 16 CygA, 16 Cyg B 
        extra_a1=[Active_Sun['a1'][0], CygA['a1'][0], CygB_higha4['a1'][0]]
        extra_epsi=[Active_Sun['epsilon'][0], CygA['epsilon'][0] , CygB_higha4['epsilon'][0]]
        extra_epsi_err=[Active_Sun['epsilon'][1], CygA['epsilon'][1] , CygB_higha4['epsilon'][1]]
        extra_epsi_err=[Active_Sun['epsilon'][2], CygA['epsilon'][2] , CygB_higha4['epsilon'][2]]
        ax[0].errorbar(extra_a1, extra_epsi, yerr=extra_epsi_err, marker='o',linestyle='', color=cols[2])
        #ax[0].axhline(0, linestyle='--', color='black')
        ax[0].set_ylabel(r'$\epsilon$')
        ax[0].set_yscale('log')
        for p in pos_bad_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['lightgray', 'lightgray', 'gray'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_fair_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['blue', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        for p in pos_good_Alm:
            stats=[data_ajAlm[p, i_theta_min] - data_ajAlm[p, i_theta_min+1], data_ajAlm[p, i_theta_min], data_ajAlm[p, i_theta0], data_ajAlm[p, i_theta_max], data_ajAlm[p, i_theta_max] + data_ajAlm[p, i_theta_max+1] ]
            Draw_BoxAndWiskers_vertical(a1past_sorted_Alm[p], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # The Active Sun
        stats=[Active_Sun['theta_min'][0] - Active_Sun['theta_min'][1] , Active_Sun['theta_min'][0], Active_Sun['theta0'][0], Active_Sun['theta_max'][0], Active_Sun['theta_max'][0] + Active_Sun['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(Active_Sun['a1'][0], stats, ax=ax[1], width=15, color=['gold', 'gold', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg A
        stats=[CygA['theta_min'][0] - CygA['theta_min'][1] , CygA['theta_min'][0], CygA['theta0'][0], CygA['theta_max'][0], CygA['theta_max'][0] + CygA['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygA['a1'][0], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        # 16 Cyg B HIGH A4 (HIGHEST UNCERTAINTY)
        stats=[CygB_higha4['theta_min'][0] - CygB_higha4['theta_min'][1] , CygB_higha4['theta_min'][0], CygB_higha4['theta0'][0], CygB_higha4['theta_max'][0], CygB_higha4['theta_max'][0] + CygB_higha4['theta_max'][1] ]
        Draw_BoxAndWiskers_vertical(CygB_higha4['a1'][0], stats, ax=ax[1], width=15, color=['green', 'blue', 'Black'], fill=True)#, color=color, fill=fill, extra=extra, show_stats=True, ax=ax[1])
        #ax[1].set_xlim(5200, 7200)
        ax[1].set_ylabel(r'Activity band $\Delta\theta$')
        ax[1].set_ylim(0, 91)
        #ax[2].axhline(0, linestyle='--', color='black')
        ax[1].set_xlabel(r'$a_1$ (nHz)')
        # Handling legends
        ax[0].legend(fontsize=10, loc='upper left')	
        fig_1d.savefig(dir_out+ 'summary_Alm_thetaBand_a1_GoodFairBad.jpg', dpi=300)
    else:
        print('WARNING: No activity summary file provided. Plots skipped!')

    print('Image files saved in: ', dir_out)

# Main part of the code
dir_core='/Users/obenomar/Work/tmp/paper2_asphericity/'
IDs_file=dir_core+ '/RealData/id_list.txt'
aj_summary_file=dir_core + '/RealData/stars.summary.ajfit'
#activity_summary_file='/Users/obenomar/Work/trash/stars.summary.ajAlm'
activity_summary_file='/Users/obenomar/Work/trash/stars.summary.14Dec2022.ajAlm'
activity_summary_file='/Users/obenomar/Work/trash/stars.summary.15Dec2022-MORESAMPLES.ajAlm'
aj_classification_file=dir_core + '/RealData/stars.classification.ajfit'
show_aj(dir_core, IDs_file, aj_summary_file, aj_classification_file, activity_summary_file=activity_summary_file)

print('All tasks done')