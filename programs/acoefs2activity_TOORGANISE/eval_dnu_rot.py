import numpy as np
from fit_a2sig import nu_CF, nu_AR
from fit_a2sig import read_mcmcobs
from acoefs import Pslm
file='/Users/obenomar/Work/tmp/test_a2AR/tmp/Realdata/products/kplr012069424_kasoc-wpsd_slc_v1_a2a3a4_nol3/aj_raw.txt'

en, el, nu_nl, a1, a2, a3, a4, a5, a6, sig_a1, sig_a2, sig_a3, sig_a4, sig_a5, sig_a6=read_mcmcobs(file)
a1=a1[0]
a3=a3[0]

Dnu=103.35
epsilon_nl=5e-4
theta0=71*np.pi/180.
delta=4*np.pi/180.
dnuCF=[]
dnuAR=[]
dnuROT=[]
dnu_nl=[]

print('a1 = ', a1)
print('a3 = ', a3)
print('nu_nl = ', nu_nl)
print(' ----- ')
print('')

for i in range(len(en)):
	# Centrifugal effect
	for m in range(-el[i],el[i]+1):
		nu=nu_CF(nu_nl[i], Dnu, a1, el[i], m, a1_unit='nHz')
		if m !=0:
			dnuCF.append((nu)/m)
	# Activity effect
	nu=nu_AR(nu_nl[i], epsilon_nl, theta0, delta, 'gate', el[i])
	for m in range(-el[i],el[i]+1):
		if m !=0:
			dnuAR.append((nu)/m)
	# Rotation
	for m in range(-el[i],el[i]+1):
		if m !=0:
			if el[i]==1:
				dnuROT.append(a1)
			else:
				dnuROT.append(a1 + a3*Pslm(3,el[i],m)/m)				

# The sum of all dnu_nl/m
for j in range(len(dnuROT)):
	dnu_nl.append(dnuCF[j]+dnuAR[j]+dnuROT[j])

print('dnuAR:', dnuAR)
print(' ----- ')
print('')

print('dnuCF:', dnuCF)
print(' ----- ')
print('')
print('dnuROT:', dnuROT)
print(' ----- ')
print('')

print('dnu_nl:', dnu_nl)
print(' ----- ')
print('')

# The mean of all dnu_nl:
print("Mean dnu_nl:", np.mean(dnu_nl))