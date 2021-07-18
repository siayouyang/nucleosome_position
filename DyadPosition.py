#author: siayouyang
#refer to Albert et al. (2007)Nature, Albert et al. (2008)Bioinformatics, Bhardwaj et al. (2020)Nat. Commun.
#WARNING: if the program has been interrupted, remove all the .temp data manually before rerun.
import os
from optparse import OptionParser
import math
import numpy

##OptionParser
parser = OptionParser()
parser.add_option("-i", "--input",dest="input", help="centered input data in wiggle format(bamCoverage --MNase --> bigWigToWig)")
parser.add_option("-s","--sigma", dest="sigma", type="float", default=20,help="Sigma to use when smoothing reads to call peaks. Default %default")
parser.add_option("-e", "--exclusion", dest="exclusion", type="int", default=147, help="Exclusion zone around each peak that prevents others from being called. Default %default.")
parser.add_option("-F", "--filter", dest="filter", type="float", default=0.0001, help="outputs only peaks with larger peak height. Default %default.")
(options, args) = parser.parse_args()

input = options.input
sigma = options.sigma
exclusion = options.exclusion
filter = options.filter

def clean_header(input):
    #remove lines start with "#"
    wig_file_temp = open(input, 'r')
    wig_file_temp_readlines = wig_file_temp.readlines()
    wig_file = open(input + "_nohead.temp", 'w')

    for r in wig_file_temp_readlines:
        if r.startswith("#"):
            pass
        else:
            wig_file.write(r)
    wig_file_temp.close()
    wig_file.close()

clean_header(input)


def split_chr_bp():
    # split lines to 1bp
    wig_file = open(input + "_nohead.temp", 'r')
    wig_file_readlines = wig_file.readlines()
    wig_file = open(input + "_split_chr.temp", 'w')
    for a in range(0, len(wig_file_readlines)):
        minus = int((wig_file_readlines[a]).split()[2]) - int((wig_file_readlines[a]).split()[1])
        b = int((wig_file_readlines[a]).split()[1])
        if minus > 1:
            for c in range(b, b + minus):
                wig_file.write(
                    f'{(wig_file_readlines[a]).split()[0]}\t{c}\t{c + 1}\t{(wig_file_readlines[a]).split()[3]}\n')
        elif minus == 1:
            wig_file.write(wig_file_readlines[a])
    wig_file.close()
    os.remove(input + "_nohead.temp")

split_chr_bp()


chr = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII","chrXIII", "chrXIV", "chrXV", "chrXVI"]
def split_chr():
    wig_file = open(input + "_split_chr.temp", 'r')
    wig_file_readlines = wig_file.readlines()
    for a in range(0, len(wig_file_readlines)):
        for b in range(0, len(chr)):
            if str(wig_file_readlines[a].split()[0]) == str(chr[b]):
                chr_file = open(input + "_split_"+ str(chr[b]) +".temp", 'a')
                chr_file.write(wig_file_readlines[a])
            else:
                pass
        chr_file.close()
    wig_file.close()
    os.remove(input + "_split_chr.temp")

split_chr()


#smoothing
def gaussian_smoothing(sigma, input, chr_a):
    width=(4 * sigma)
    wig_file = open(input + "_split_"+ str(chr_a) +".temp", 'r')
    wig_file_readlines = wig_file.readlines()
    smooth_file = open(input + "_smooth_"+ str(chr_a) +".temp", 'w')

    def normal_func(x):
        return math.exp(-x * x / (2 * (sigma**2)))

    gaussian = list(map(normal_func, range(-width, width)))
    gaussian = numpy.array(gaussian, numpy.float)
    # normalization
    gaussian_array = 1.0 / math.sqrt(2 * numpy.pi * (sigma**2)) * gaussian

    for a in range(0, len(wig_file_readlines)):
        raw_array=[]
        for b in range((a-width),(a+width)):
            if b < 0:
                raw_array.append(float(0))
            elif b >= (len(wig_file_readlines)):
                raw_array.append(float(0))
            else:
                raw_array.append(float(wig_file_readlines[b].split()[3]))
        smooth_array=[]
        for c in range(0, len(gaussian_array)):
            d = float(gaussian_array[c])*float(raw_array[c])
            smooth_array.append(d)
        smooth_value = sum(smooth_array)
        smooth_file.write(f'{wig_file_readlines[a].split()[0]}\t{wig_file_readlines[a].split()[1]}\t{wig_file_readlines[a].split()[2]}\t{smooth_value}\n')
    wig_file.close()
    smooth_file.close()

def choose_chr_to_smooth():
    for a in range(0, len(chr)):
        if os.path.isfile(input + "_split_"+ chr[a] +".temp"):
            gaussian_smoothing(sigma, input, chr[a])
        else:
            pass

choose_chr_to_smooth()


#find peaks
def find_peaks(chr_a):
    smooth_file = open(input + "_smooth_" + str(chr_a) + ".temp", 'r')
    smooth_file_readlines = smooth_file.readlines()
    peaks_file = open(input + "_peaks_" + str(chr_a) + ".temp", 'w')
    for a in range(0, len(smooth_file_readlines)):
        if a == 0:
            pass
        elif a == (int(len(smooth_file_readlines))-1):
            pass
        else:
            if (float(smooth_file_readlines[a].split()[3]) > float(smooth_file_readlines[a-1].split()[3])) and (float(smooth_file_readlines[a].split()[3]) > float(smooth_file_readlines[a+1].split()[3])):
                if float(smooth_file_readlines[a].split()[3]) > filter:
                    peaks_file.write(smooth_file_readlines[a])
                else:
                    pass
            else:
                pass
    smooth_file.close()
    peaks_file.close()

def choose_chr_to_find_peaks():
    for a in range(0, len(chr)):
        if os.path.isfile(input + "_smooth_" + str(chr[a]) + ".temp"):
            find_peaks(chr[a])
        else:
            pass

choose_chr_to_find_peaks()


#stddev
def stddev(chr_a):
    peaks_file = open(input + "_peaks_" + str(chr_a) + ".temp", 'r')
    peaks_file_readlines = peaks_file.readlines()
    wig_file = open(input + "_split_" + str(chr_a) + ".temp", 'r')
    wig_file_readlines = wig_file.readlines()
    std_file = open(input + "_std_" + str(chr_a) + ".temp", 'w')
    for a in range(0, len(peaks_file_readlines)):
        peak_range = numpy.arange((int(peaks_file_readlines[a].split()[1])-(2*sigma)),(int(peaks_file_readlines[a].split()[1])+(2*sigma)))
        position_array = []
        frequency_array = []
        for b in peak_range:
            if b < 0:
                pass
            elif b >= (len(wig_file_readlines)):
                pass
            else:
                position_array.append(int(wig_file_readlines[b].split()[1]))
                frequency_array.append(float(wig_file_readlines[b].split()[3]))
        flat_list = []
        for c in range(0, len(position_array)):
            position = int(position_array[c])
            frequency = int(float(frequency_array[c])*10000)
            n = 0
            while n < frequency:
                flat_list.append(position)
                n += 1
        peaks_std = numpy.std(flat_list)
        std_file.write(f'{peaks_file_readlines[a].split()[0]}\t{peaks_file_readlines[a].split()[1]}\t{peaks_file_readlines[a].split()[2]}\t{peaks_file_readlines[a].split()[3]}\t{peaks_std}\n')
    peaks_file.close()
    wig_file.close()
    std_file.close()
    os.remove(input + "_peaks_" + str(chr_a) + ".temp")
    os.remove(input + "_split_" + str(chr_a) + ".temp")

def choose_chr_stddev():
    for a in range(0, len(chr)):
        if os.path.isfile(input + "_peaks_" + str(chr[a]) + ".temp"):
            stddev(chr[a])
        else:
            pass

choose_chr_stddev()


#exclusion
def perform_exclusion(chr_a):
    std_file = open(input + "_std_" + str(chr_a) + ".temp", 'r')
    std_file_readlines = std_file.readlines()
    exclude_file = open(input + "_exclude_" + str(chr_a) + ".temp", 'w')

    #create dictionary
    list_pos = []
    list_val = []
    for a in range(0, len(std_file_readlines)):
        list_pos.append(int(std_file_readlines[a].split()[1]))
        list_val.append(float(std_file_readlines[a].split()[3]))

    dict = {}
    for b in range(0, len(list_pos)):
        dict[list_pos[b]] = list_val[b]

    safe_keys = []
    while len(dict.keys()) > 0:
        max_value = (max(dict.values()))

        def get_key(dict, value):
            return [k for k, v in dict.items() if v == value]

        max_key = get_key(dict, max_value)
        safe_keys.append(max_key[0])

        exclusion_range = numpy.arange(int(max_key[0]) - int(exclusion//2), int(max_key[0]) + int(exclusion//2))
        for r in exclusion_range:
            if r in dict.keys():
                del dict[r]
            else:
                pass

    for c in range(0, len(std_file_readlines)):
        if int(std_file_readlines[c].split()[1]) in safe_keys:
            exclude_file.write(std_file_readlines[c])
        else:
            pass
    std_file.close()
    exclude_file.close()
    #os.remove(input + "_std_" + str(chr_a) + ".temp")

def choose_chr_perform_exclusion():
    for a in range(0, len(chr)):
        if os.path.isfile(input + "_std_" + str(chr[a]) + ".temp"):
            perform_exclusion(chr[a])
        else:
            pass

choose_chr_perform_exclusion()


#merge chr files
def merging(chr_a):
    smooth_file = open(input + "_smooth_" + str(chr_a) + ".temp", 'r')
    smooth_file_readlines = smooth_file.readlines()
    exclude_file = open(input + "_exclude_" + str(chr_a) + ".temp", 'r')
    exclude_file_readlines = exclude_file.readlines()
    smooth_merged_file = open(input + "_smooth_merged_.wig", 'a')
    exclude_merged_file = open(input + "_exclude_merged_.bed", 'a')
    for s in range(0, len(smooth_file_readlines)):
        smooth_merged_file.write(smooth_file_readlines[s])
    for e in range(0, len(exclude_file_readlines)):
        exclude_merged_file.write(exclude_file_readlines[e])
    smooth_file.close()
    smooth_merged_file.close()
    exclude_file.close()
    exclude_merged_file.close()
    os.remove(input + "_smooth_" + str(chr_a) + ".temp")
    os.remove(input + "_exclude_" + str(chr_a) + ".temp")

def choose_chr_to_merge():
    for a in range(0, len(chr)):
        if os.path.isfile(input + "_exclude_" + str(chr[a]) + ".temp"):
            merging(chr[a])
        else:
            pass

choose_chr_to_merge()

print("Done!")