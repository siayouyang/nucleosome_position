#author: siayouyang
#refer to Albert et al. (2007)Nature, Albert et al. (2008)Bioinformatics, Bhardwaj et al. (2020)Nat. Commun., Bakel et al., (2013)PLoS Genet.
#output format:[0]chr,[1]chr_start,[2]chr_end,[3]gene,[4]score,[5]strand,[6-7]+1_start,end,[8-9]+2_start,end,[10-11]+3_start,end,[12-13]+4_start,end
import os
from optparse import OptionParser

##OptionParser
parser = OptionParser()
parser.add_option("-n", "--nucleosome", dest="nucleosome_position", help="nucleosome positions data in BED format")
parser.add_option("-t", "--tss", dest="TSS", help="list of TSS coordinates derived from van Bakel et al., 2013")
(options, args) = parser.parse_args()

nucleosome_position = options.nucleosome_position
TSS = options.TSS

#genes >= 560bp
def filter_gene_length():
    TSS_raw_file = open(TSS, 'r')
    TSS_raw_file_readlines = TSS_raw_file.readlines()
    TSS_file = open(TSS+"_clean.txt" , 'w')
    for a in range(0, len(TSS_raw_file_readlines)):
        if abs(int(TSS_raw_file_readlines[a].split()[4])-int(TSS_raw_file_readlines[a].split()[3])) >= 560:
            TSS_file.write(TSS_raw_file_readlines[a])
        else:
            pass
    TSS_raw_file.close()
    TSS_file.close()

filter_gene_length()

def assign_nucleosome_position():
    TSS_file = open(TSS+"_clean.txt" , 'r')
    TSS_file_readlines = TSS_file.readlines()
    TSS_file2 = open(nucleosome_position + "_4Nuc.txt", 'w')
    nucleosome_position_file = open(nucleosome_position, 'r')
    nucleosome_position_file_readlines = nucleosome_position_file.readlines()

    for a in range(0, len(TSS_file_readlines)):
        if TSS_file_readlines[a].split()[6] == "+":
            TSS_position = int(TSS_file_readlines[a].split()[3])
        elif TSS_file_readlines[a].split()[6] == "-":
            TSS_position = int(TSS_file_readlines[a].split()[4])
        for b in range(0, len(nucleosome_position_file_readlines)):
            if str(nucleosome_position_file_readlines[b].split()[0]) == str(TSS_file_readlines[a].split()[0]):
                distance0 = 10000000000
                distance1 = abs(int(nucleosome_position_file_readlines[b].split()[1]) - TSS_position)
                c = 1
                while distance1 < distance0:
                    if (b+c) < int(len(nucleosome_position_file_readlines)):
                        if str(nucleosome_position_file_readlines[b+c].split()[0]) == str(TSS_file_readlines[a].split()[0]):
                            distance0 = distance1
                            distance1 = abs(int(nucleosome_position_file_readlines[b+c].split()[1]) - TSS_position)
                            c += 1
                        else:
                            break
                    else:
                        break
                if TSS_file_readlines[a].split()[6] == "+":
                    if abs(TSS_position - int(nucleosome_position_file_readlines[b+c-2].split()[1])) <= 150:
                        plus_1_pos = int(nucleosome_position_file_readlines[b+c-2].split()[1])
                        if (abs(int(nucleosome_position_file_readlines[b+c-2].split()[1])-int(nucleosome_position_file_readlines[b+c-1].split()[1])) <= 300) and (str(nucleosome_position_file_readlines[b+c-1].split()[0])==str(TSS_file_readlines[a].split()[0])):
                            plus_2_pos = int(nucleosome_position_file_readlines[b + c - 1].split()[1])
                            if (abs(int(nucleosome_position_file_readlines[b + c - 1].split()[1]) - int(nucleosome_position_file_readlines[b + c].split()[1])) <= 300) and (str(nucleosome_position_file_readlines[b + c].split()[0]) == str(TSS_file_readlines[a].split()[0])):
                                plus_3_pos = int(nucleosome_position_file_readlines[b + c].split()[1])
                                if (abs(int(nucleosome_position_file_readlines[b + c].split()[1]) - int(nucleosome_position_file_readlines[b + c + 1].split()[1])) <= 300) and (str(nucleosome_position_file_readlines[b + c + 1].split()[0]) == str(TSS_file_readlines[a].split()[0])):
                                    plus_4_pos = int(nucleosome_position_file_readlines[b + c + 1].split()[1])
                                    chr = TSS_file_readlines[a].split()[0]
                                    chr_start = TSS_file_readlines[a].split()[3]
                                    chr_end = TSS_file_readlines[a].split()[4]
                                    gene = TSS_file_readlines[a].split()[8]
                                    score = 0
                                    strand = TSS_file_readlines[a].split()[6]
                                    TSS_file2.write(f'{chr}\t{chr_start}\t{chr_end}\t{gene}\t{score}\t{strand}\t{plus_1_pos}\t{plus_1_pos+1}\t{plus_2_pos}\t{plus_2_pos+1}\t{plus_3_pos}\t{plus_3_pos+1}\t{plus_4_pos}\t{plus_4_pos+1}\n')
                                    break
                                else:
                                    pass
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass
                elif TSS_file_readlines[a].split()[6] == "-":
                    if abs(TSS_position - int(nucleosome_position_file_readlines[b+c-2].split()[1])) <= 150:
                        plus_1_pos = int(nucleosome_position_file_readlines[b+c-2].split()[1])
                        if (abs(int(nucleosome_position_file_readlines[b+c-2].split()[1])-int(nucleosome_position_file_readlines[b+c-3].split()[1])) <= 300) and (str(nucleosome_position_file_readlines[b+c-3].split()[0])==str(TSS_file_readlines[a].split()[0])):
                            plus_2_pos = int(nucleosome_position_file_readlines[b + c - 3].split()[1])
                            if (abs(int(nucleosome_position_file_readlines[b + c - 3].split()[1]) - int(nucleosome_position_file_readlines[b + c - 4].split()[1])) <= 300) and (str(nucleosome_position_file_readlines[b + c -4].split()[0]) == str(TSS_file_readlines[a].split()[0])):
                                plus_3_pos = int(nucleosome_position_file_readlines[b + c - 4].split()[1])
                                if (abs(int(nucleosome_position_file_readlines[b + c - 4].split()[1]) - int(nucleosome_position_file_readlines[b + c - 5].split()[1])) <= 300) and (str(nucleosome_position_file_readlines[b + c - 5].split()[0]) == str(TSS_file_readlines[a].split()[0])):
                                    plus_4_pos = int(nucleosome_position_file_readlines[b + c - 5].split()[1])
                                    chr = TSS_file_readlines[a].split()[0]
                                    chr_start = TSS_file_readlines[a].split()[3]
                                    chr_end = TSS_file_readlines[a].split()[4]
                                    gene = TSS_file_readlines[a].split()[8]
                                    score = 0
                                    strand = TSS_file_readlines[a].split()[6]
                                    TSS_file2.write(f'{chr}\t{chr_start}\t{chr_end}\t{gene}\t{score}\t{strand}\t{plus_1_pos}\t{plus_1_pos+1}\t{plus_2_pos}\t{plus_2_pos+1}\t{plus_3_pos}\t{plus_3_pos+1}\t{plus_4_pos}\t{plus_4_pos+1}\n')
                                    break
                                else:
                                    pass
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                pass
    TSS_file.close()
    TSS_file2.close()
    nucleosome_position_file.close()
    os.remove(TSS+"_clean.txt" )

assign_nucleosome_position()

print("Done!")

