#author: siayouyang
#refer to Albert et al. (2007)Nature, Albert et al. (2008)Bioinformatics, Bhardwaj et al. (2020)Nat. Commun.
from optparse import OptionParser
import numpy

##OptionParser
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="data with position shift and p-value")
parser.add_option("-t", "--threshold", dest="threshold", type="float", default=0.01, help="significant threshold. Default %default.")
(options, args) = parser.parse_args()
input = options.input
threshold = float(options.threshold)

def get_significant_data():
    input_file = open(input, 'r')
    input_file_readlines = input_file.readlines()
    significant_file = open(input+"_sig.bed", 'w')
    shift1_list = []
    shift2_list = []
    shift3_list = []
    shift4_list = []
    num = 0
    for a in range(0, len(input_file_readlines)):
        shift1 = int(input_file_readlines[a].split()[6])
        shift2 = int(input_file_readlines[a].split()[8])
        shift3 = int(input_file_readlines[a].split()[10])
        shift4 = int(input_file_readlines[a].split()[12])
        p_val1 = float(input_file_readlines[a].split()[7])
        p_val2 = float(input_file_readlines[a].split()[9])
        p_val3 = float(input_file_readlines[a].split()[11])
        p_val4 = float(input_file_readlines[a].split()[13])
        if (p_val1<=threshold) or (p_val2<=threshold) or (p_val3<=threshold) or (p_val4<=threshold):
            significant_file.write(input_file_readlines[a])
            shift1_list.append(shift1)
            shift2_list.append(shift2)
            shift3_list.append(shift3)
            shift4_list.append(shift4)
            num += 1
        else:
            pass
    shift1_median = numpy.median(shift1_list)
    shift2_median = numpy.median(shift1_list)
    shift3_median = numpy.median(shift1_list)
    shift4_median = numpy.median(shift1_list)
    print(f"number of genes with significant shift: {num}\n+1 shift median= {shift1_median}\n+2 shift median= {shift2_median}\n+3 shift median= {shift3_median}\n+4 shift median= {shift4_median}\n")
    input_file.close()
    significant_file.close()

get_significant_data()

print("Done!")

