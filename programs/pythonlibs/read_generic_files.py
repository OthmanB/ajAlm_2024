# Regroup functions to read generic file format here
# 
import numpy as np

def read_matrix_tab(filein, separator=" ", first_data_is_label=True, second_data_is_unit=False, datatype="string"):
	'''
		Reader of tables in a matrix format with a header marked by '#' as first character
		file: filename of the file to read
		separator: Identifies how data are separated
		first_data_is_label: If True, consider the first data line as providing the labels for the following lines
		second_data_is_unit: If True, consider the second data line as providing the units for the following lines
		datatype: If "string", it will return a list of lists with the upper level list index the line number and  lower level index the columns
                  Otherwise, it will return a numpy array of values for the main data.
	'''
	if datatype == "string":
		return read_matrix_tab_string(filein, separator=separator, first_data_is_label=first_data_is_label, second_data_is_unit=second_data_is_unit)
	else:
		return read_matrix_tab_float(filein, separator=separator, first_data_is_label=first_data_is_label, second_data_is_unit=second_data_is_unit)
	
def read_matrix_tab_float(file, separator=" ", first_data_is_label=True, second_data_is_unit=False):
	'''
		Reader of tables in a matrix format with a header marked by '#' as first character
		It will return a numpy array of float
		file: filename of the file to read
		separator: Identifies how data are separated
		first_data_is_label: If True, consider the first data line as providing the labels for the following lines
	'''
	f=open(file, encoding='utf-8-sig') # 'utf-8-sig is used instead of utf-8 because Excel seems to save with BOM as first character and this messes up everything
	raw=f.read()
	f.close()
	raw=raw.splitlines()
	Nlines=len(raw)
	header=[]
	label=[]
	unit=[]
	label_done=not first_data_is_label
	unit_done=not second_data_is_unit
	i=0
	for line in raw:
		if line != '':
			if line[0] == '#':
				header.append(line)
			else: 
				s=line.split(separator)
				if first_data_is_label == True and label_done == False:
					label=s
					label_done=True
					continue  # skip the rest of the loop 
				if second_data_is_unit == True and unit_done == False:
					unit=s
					unit_done=True
					continue # skip the rest of the loop
				if i==0:
					data=np.zeros((Nlines, len(s)), dtype=float) 
					data[i,:]=np.asarray(s, dtype=float)
				else:
					data[i,:]=np.asarray(s, dtype=float)
				i=i+1
	# Remove potential empty lines at the end
	if i-1 < Nlines:
		data=data[0:i]
	return data, header,label, unit

def read_matrix_tab_string(file, separator=" ", first_data_is_label=True, second_data_is_unit=False):
    '''
        Reader of tables in a matrix format with a header marked by '#' as first character
        Returns list of lists. The upper level list index the line number. The lower level index the columns
        file: filename of the file to read
        separator: Identifies how data are separated
        first_data_is_label: If True, consider the first data line as providing the labels for the following lines
    '''
    f=open(file, encoding='utf-8-sig') # 'utf-8-sig is used instead of utf-8 because Excel seems to save with BOM as first character and this messes up everything
    raw=f.read()
    f.close()
    raw=raw.splitlines()
    data=[]
    header=[]
    label=[]
    unit=[]
    label_done=not first_data_is_label
    unit_done=not second_data_is_unit
    for line in raw:
        line=line.strip()
        if line != '':
            if line[0] == '#':
                header.append(line)
            else:
                s=line.split(separator)
                if first_data_is_label == True and label_done == False:
                    label=s
                    label_done=True
                    continue    
                if second_data_is_unit == True and unit_done == False:
                    unit=s
                    unit_done=True
                    continue
                if label_done == True and unit_done == True:
                    data.append(s)
    return data, header,label, unit