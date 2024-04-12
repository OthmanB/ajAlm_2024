from compress_bin2tar import bin2targz

def relative_dir_tree4process():
    dir=[]
    dir.append("1001/")
    dir.append("1101/")
    dir.append("1111/")
    dir.append("1201/Gate/decompose_-1/")
    dir.append("1201/Gate/decompose_1/")
    dir.append("1201/Gate/decompose_2/")
    dir.append("1201/Triangle/decompose_-1/")
    dir.append("1201/Triangle/decompose_1/")
    dir.append("1201/Triangle/decompose_2/")
    dir.append("1211/Gate/decompose_-1/")
    dir.append("1211/Gate/decompose_1/")
    dir.append("1211/Gate/decompose_2/")
    dir.append("1211/Triangle/decompose_-1/")
    dir.append("1211/Triangle/decompose_1/")
    dir.append("1211/Triangle/decompose_2/")
    return dir
		
def compress_All(dir_core):
    #print(" I. Processing Legacy stars...")
    #dir0=dir_core + "/Legacy/"
    #relative_dirs=relative_dir_tree4process()
    #for d in relative_dirs:
    #    dir=dir0 + d
    #    bin2targz(dir)

    print(" II. Processing Kamiaka stars...")
    dir0=dir_core + "/Kamiaka2018/"
    relative_dirs=relative_dir_tree4process()
    for d in relative_dirs:
        dir=dir0 + d
        bin2targz(dir)


dir_core="/var/services/homes/obenomar/Temporary/ajAlm_project/data/TRANSFORMED_RUN/Outputs/"
print(" ------------ ", flush=True)
print(" Compressing all .bin files into individual .tar.gz files to save space...", flush=True)
print(" ------------ ", flush=True)
	
compress_All(dir_core)
