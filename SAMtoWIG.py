#author: siayouyang
#for sacCer3 genome
#for center 1bp read mapping
from optparse import OptionParser
from collections import Counter
import os
import re

##OptionParser
parser = OptionParser()
parser.add_option("-s", "--sam",dest="sam", help="input data in SAM format")
(options, args) = parser.parse_args()

sam = options.sam

def find_center():
    sam_file = open(sam, 'r')
    sam_file_readlines = sam_file.readlines()
    temp_file = open(sam + ".temp" , 'w')
    for a in range(0, len(sam_file_readlines)):
        if (len(sam_file_readlines[a].split()) > 5):
            if (str(sam_file_readlines[a].split()[0])) == str("@PG"):
                continue
            else:
                chr = sam_file_readlines[a].split()[2]
                pos = int(sam_file_readlines[a].split()[3])
                length = int(sam_file_readlines[a].split()[8])
                if length > 0:
                    temp_file.write(f'{chr}\t{pos+(length//2)-1}\t{pos+(length//2)}\t{chr}{pos+(length//2)}\n')
                else:
                    continue
        else:
            continue
    sam_file.close()
    temp_file.close()

find_center()

def count_reads():
    temp_file = open(sam + ".temp" , 'r')
    temp_file_readlines = temp_file.readlines()
    temp2_file = open(sam + "2.temp" , 'w')
    alias_list = []
    for a in range(0, len(temp_file_readlines)):
        alias_list.append(temp_file_readlines[a].split()[3])
    count_dict = Counter(alias_list)
    keys = list(dict.keys(count_dict))
    values = list(dict.values(count_dict))
    for a in range(0, len(keys)):
        split_keys = re.split('(\d+)',keys[a])
        chr = split_keys[0]
        start = int(split_keys[1])-1
        end = int(split_keys[1])
        temp2_file.write(f'{chr}\t{start}\t{end}\t{values[a]}\n')
    temp_file.close()
    temp2_file.close()
    os.remove(sam + ".temp")

count_reads()

def create_wig():
    temp2_file = open(sam + "2.temp", 'r')
    temp2_file_readlines = temp2_file.readlines()
    wig_file =  open(sam + ".wig", 'w')
    chr = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII",
           "chrXIII", "chrXIV", "chrXV", "chrXVI", "chrM"]
    chr_len = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816,1078177, 924431, 784333,1091291, 948066, 85779]
    for a in range(0, len(chr)):
        # create dictionary (from temp2)
        list_start = []
        list_val = []
        for b in range(0, len(temp2_file_readlines)):
            if str(chr[a]) == str(temp2_file_readlines[b].split()[0]):
                list_start.append(int(temp2_file_readlines[b].split()[1]))
                list_val.append(float(temp2_file_readlines[b].split()[3]))   #reads
            else:
                pass
        for c in range(0, int(chr_len[a])):
            if int(c) in list_start:
                idx = list_start.index(int(c))
                start = int(c)
                end = int(c)+1
                val = list_val[idx]
                wig_file.write(f'{chr[a]}\t{start}\t{end}\t{val}\n')
            else:
                wig_file.write(f'{chr[a]}\t{int(c)}\t{int(c)+1}\t0\n')
    temp2_file.close()
    wig_file.close()
    os.remove(sam + "2.temp")

create_wig()


print("Done!")
