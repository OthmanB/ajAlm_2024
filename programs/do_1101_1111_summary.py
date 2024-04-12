from do_main_summary import *

cpp_path="/var/services/homes/obenomar/TAMCMC_bin_1.86.78_x86/"
recompute_evidence=[False, True, cpp_path]
data_dir="/var/services/homes/dataonly/Kepler_FullScale_Analysis/Sun-like/Outputs/Tcoef_1.7_1.86.77/"
tmpdir="tmp/"
# QUESTION: DO WE HAVE SIGNIFICANT a3 irrespected of the shape of the activity?
dir_out="/var/services/homes/dataonly/Kepler_FullScale_Analysis/Sun-like/Products/Tcoef_1.7_1.86.77/statistical_summary/1101_vs_1111/"
done_modeltype=[1101,1111]
do_products_database(data_dir, dir_out, 
                    data_sources=["Legacy", "Kamiaka2018"],
                    done_modeltype=done_modeltype,
                    done_modelsubtype1=["Gate", "Triangle"],
                    done_modelsubtype2=["decompose_-1", "decompose_1", "decompose_2"],
                    cpp_path=cpp_path,
                    recompute_evidence=recompute_evidence, phase="A", period=60, tmpdir=tmpdir)
