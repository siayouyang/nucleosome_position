#author: siayouyang
#refer to Albert et al. (2007)Nature, Albert et al. (2008)Bioinformatics, Bhardwaj et al. (2020)Nat. Commun.
#output format:[0]chr,[1]chr_start,[2]chr_end,[3]gene,[4]score,[5]strand,[6-7]+1_shift,p-val,[8-9]+2_shift,p-val,[10-11]+3_shift,p-val,[12-13]+4_shift,p-val
import os
from optparse import OptionParser
from scipy.stats import ttest_ind_from_stats

##OptionParser
parser = OptionParser()
parser.add_option("--input1", dest="input1", help="WT data with +1 to +4 nucleosome positions")
parser.add_option("--input2", dest="input2", help="mutant data with +1 to +4 nucleosome positions")
parser.add_option("--nuc1", dest="nucleosome_pos1", help="WT nucleosome positions data")
parser.add_option("--nuc2", dest="nucleosome_pos2", help="mutant nucleosome positions data")
(options, args) = parser.parse_args()

input1 = options.input1
input2 = options.input2
nucleosome_pos1 = options.nucleosome_pos1
nucleosome_pos2 = options.nucleosome_pos2

#select common genes
def select_common():
    input1_file = open(input1, 'r')
    input1_file_readlines = input1_file.readlines()
    input2_file = open(input2, 'r')
    input2_file_readlines = input2_file.readlines()
    common1_file = open(input1+input2+"_common.txt", 'w')
    common2_file = open(input2+input1+"_common.txt", 'w')
    for a in range(0, len(input1_file_readlines)):
        for b in range(0, len(input2_file_readlines)):
            chr1 = input1_file_readlines[a].split()[0]
            chr_start1 = input1_file_readlines[a].split()[1]
            chr_end1 = input1_file_readlines[a].split()[2]
            gene1 = input1_file_readlines[a].split()[3]
            chr2 = input2_file_readlines[b].split()[0]
            chr_start2 = input2_file_readlines[b].split()[1]
            chr_end2 = input2_file_readlines[b].split()[2]
            gene2 = input2_file_readlines[b].split()[3]
            if (str(chr1) == str(chr2)) and (int(chr_start1) == int(chr_start2)) and (int(chr_end1) == int(chr_end2)) and (str(gene1) == str(gene2)):
                common1_file.write(input1_file_readlines[a])
                common2_file.write(input2_file_readlines[b])
            else:
                pass
    input1_file.close()
    input2_file.close()
    common1_file.close()
    common2_file.close()


select_common()

def find_shift_pval():
    common1_file = open(input1+input2+"_common.txt", 'r')
    common1_file_readlines = common1_file.readlines()
    common2_file = open(input2+input1+"_common.txt", 'r')
    common2_file_readlines = common2_file.readlines()
    nucpos1_file = open(nucleosome_pos1, 'r')
    nucpos1_file_readlines = nucpos1_file.readlines()
    nucpos2_file = open(nucleosome_pos2, 'r')
    nucpos2_file_readlines = nucpos2_file.readlines()
    shift_file = open(input2+"_shift.txt", 'w')
    for a in range(0, len(common1_file_readlines)):
        common1_chr = common1_file_readlines[a].split()[0]
        common1_plus1 = common1_file_readlines[a].split()[6]
        common1_plus2 = common1_file_readlines[a].split()[8]
        common1_plus3 = common1_file_readlines[a].split()[10]
        common1_plus4 = common1_file_readlines[a].split()[12]
        common2_chr = common2_file_readlines[a].split()[0]
        common2_plus1 = common2_file_readlines[a].split()[6]
        common2_plus2 = common2_file_readlines[a].split()[8]
        common2_plus3 = common2_file_readlines[a].split()[10]
        common2_plus4 = common2_file_readlines[a].split()[12]
        common_strand = common1_file_readlines[a].split()[5]
        for b in range(0, len(nucpos1_file_readlines)):
            nuc1_chr = nucpos1_file_readlines[b].split()[0]
            nuc1_start = nucpos1_file_readlines[b].split()[1]
            nuc1_std = nucpos1_file_readlines[b].split()[4]
            nuc1_mean = nucpos1_file_readlines[b].split()[5]
            nuc1_readcount = nucpos1_file_readlines[b].split()[6]
            if (str(common1_chr) == str(nuc1_chr)) and (int(common1_plus1) == int(nuc1_start)):
                common1_plus1_std = nuc1_std
                common1_plus1_mean = nuc1_mean
                common1_plus1_readcount = nuc1_readcount
            elif (str(common1_chr) == str(nuc1_chr)) and (int(common1_plus2) == int(nuc1_start)):
                common1_plus2_std = nuc1_std
                common1_plus2_mean = nuc1_mean
                common1_plus2_readcount = nuc1_readcount
            elif (str(common1_chr) == str(nuc1_chr)) and (int(common1_plus3) == int(nuc1_start)):
                common1_plus3_std = nuc1_std
                common1_plus3_mean = nuc1_mean
                common1_plus3_readcount = nuc1_readcount
            elif (str(common1_chr) == str(nuc1_chr)) and (int(common1_plus4) == int(nuc1_start)):
                common1_plus4_std = nuc1_std
                common1_plus4_mean = nuc1_mean
                common1_plus4_readcount = nuc1_readcount
            else:
                pass
        for c in range(0, len(nucpos2_file_readlines)):
            nuc2_chr = nucpos2_file_readlines[c].split()[0]
            nuc2_start = nucpos2_file_readlines[c].split()[1]
            nuc2_std = nucpos2_file_readlines[c].split()[4]
            nuc2_mean = nucpos2_file_readlines[c].split()[5]
            nuc2_readcount = nucpos2_file_readlines[c].split()[6]
            if (str(common2_chr) == str(nuc2_chr)) and (int(common2_plus1) == int(nuc2_start)):
                common2_plus1_std = nuc2_std
                common2_plus1_mean = nuc2_mean
                common2_plus1_readcount = nuc2_readcount
            elif (str(common2_chr) == str(nuc2_chr)) and (int(common2_plus2) == int(nuc2_start)):
                common2_plus2_std = nuc2_std
                common2_plus2_mean = nuc2_mean
                common2_plus2_readcount = nuc2_readcount
            elif (str(common2_chr) == str(nuc2_chr)) and (int(common2_plus3) == int(nuc2_start)):
                common2_plus3_std = nuc2_std
                common2_plus3_mean = nuc2_mean
                common2_plus3_readcount = nuc2_readcount
            elif (str(common2_chr) == str(nuc2_chr)) and (int(common2_plus4) == int(nuc2_start)):
                common2_plus4_std = nuc2_std
                common2_plus4_mean = nuc2_mean
                common2_plus4_readcount = nuc2_readcount
            else:
                pass
        if str(common_strand) == "+":
            shift_plus1 = int(common2_plus1) - int(common1_plus1)
            shift_plus2 = int(common2_plus2) - int(common1_plus2)
            shift_plus3 = int(common2_plus3) - int(common1_plus3)
            shift_plus4 = int(common2_plus4) - int(common1_plus4)
        elif str(common_strand) == "-":
            shift_plus1 = int(common1_plus1) - int(common2_plus1)
            shift_plus2 = int(common1_plus2) - int(common2_plus2)
            shift_plus3 = int(common1_plus3) - int(common2_plus3)
            shift_plus4 = int(common1_plus4) - int(common2_plus4)
        t_plus1 = ttest_ind_from_stats(float(common1_plus1_mean), float(common1_plus1_std), int(common1_plus1_readcount), float(common2_plus1_mean), float(common2_plus1_std), int(common2_plus1_readcount), equal_var=False, alternative='two-sided')
        t_plus2 = ttest_ind_from_stats(float(common1_plus2_mean), float(common1_plus2_std), int(common1_plus2_readcount), float(common2_plus2_mean), float(common2_plus2_std), int(common2_plus2_readcount), equal_var=False, alternative='two-sided')
        t_plus3 = ttest_ind_from_stats(float(common1_plus3_mean), float(common1_plus3_std), int(common1_plus3_readcount), float(common2_plus3_mean), float(common2_plus3_std), int(common2_plus3_readcount), equal_var=False, alternative='two-sided')
        t_plus4 = ttest_ind_from_stats(float(common1_plus4_mean), float(common1_plus4_std), int(common1_plus4_readcount), float(common2_plus4_mean), float(common2_plus4_std), int(common2_plus4_readcount), equal_var=False, alternative='two-sided')
        p_value1 = float(t_plus1[1])
        p_value2 = float(t_plus2[1])
        p_value3 = float(t_plus3[1])
        p_value4 = float(t_plus4[1])
        shift_file.write(f'{common1_file_readlines[a].split()[0]}\t{common1_file_readlines[a].split()[1]}\t{common1_file_readlines[a].split()[2]}\t{common1_file_readlines[a].split()[3]}\t{common1_file_readlines[a].split()[4]}\t{common1_file_readlines[a].split()[5]}\t{shift_plus1}\t{p_value1}\t{shift_plus2}\t{p_value2}\t{shift_plus3}\t{p_value3}\t{shift_plus4}\t{p_value4}\n')
    common1_file.close()
    common2_file.close()
    nucpos1_file.close()
    nucpos2_file.close()
    shift_file.close()


find_shift_pval()
os.remove(input2+input1+"_common.txt")
os.remove(input1+input2+"_common.txt")
print("Done")