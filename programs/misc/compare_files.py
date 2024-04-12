def compare_files(file1, file2, column1, column2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = [line for line in f1 if not line.startswith('#')]
        lines2 = [line for line in f2 if not line.startswith('#')]

        delimiter1 = '\t' if '\t' in lines1[0] else ' '  # Check if first line in file1 has a tabulation delimiter
        delimiter2 = '\t' if '\t' in lines2[0] else ' '  # Check if first line in file2 has a tabulation delimiter

        column1_file1 = [line.split(delimiter1)[column1] for line in lines1]
        column2_file2 = [line.split(delimiter2)[column2] for line in lines2]

        differences = set(column1_file1) - set(column2_file2)
        return differences

file1 = '../../data/stars_list_Kamiaka2018.txt'
file2 = '../../data/stars_list_done.txt'
column1 = 0  # Specify the column in file1 to compare (starting from 0)
column2 = 0  # Specify the column in file2 to compare (starting from 0)
differences = compare_files(file1, file2, column1, column2)
print(differences)
