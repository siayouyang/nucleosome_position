#author: siayouyang
#for PE SAM format data
#Usage(for clean header and minus strand):python RandomSelect.py -s sample.sam -c
#Usage(for random select):python RandomSelect.py -s sample.sam_clean.sam -r -n 1000000
from optparse import OptionParser
import random

##OptionParser
parser = OptionParser()
parser.add_option("-s", "--sam",dest="sam", help="input data in SAM format")
parser.add_option("-c", "--clean", dest="clean", action="store_true", help="clean header and remove minus strand")
parser.add_option("-r", "--random", dest="rand", action="store_true", help="random select specific number of reads")
parser.add_option("-n", "--number", dest="number", type="int", default=0, help="reads number")
(options, args) = parser.parse_args()

sam = options.sam
clean = options.clean
rand = options.rand
number = options.number

#exclude header and minus strand
def clean_sam():
    sam_file = open(sam, 'r')
    sam_file_readlines = sam_file.readlines()
    clean_file = open(sam + "_clean.sam" , 'w')
    for a in range(0, len(sam_file_readlines)):
        if (len(sam_file_readlines[a].split()) > 5):
            if (str(sam_file_readlines[a].split()[0])) == str("@PG"):
                continue
            else:
                if int(sam_file_readlines[a].split()[8]) < 0:
                    continue
                else:
                    clean_file.write(sam_file_readlines[a])
        else:
            continue
    sam_file.close()
    clean_file.close()


if clean:
    clean_sam()


def random_select():
    sam_file = open(sam, 'r')
    sam_file_readlines = sam_file.readlines()
    selected_file = open(sam + "_selected.sam", 'w')
    line_list = []
    for a in range(0, len(sam_file_readlines)):
        line_list.append(int(a))
    random.shuffle(line_list)
    for b in range(0, number):
        line = int(line_list[b])
        selected_file.write(sam_file_readlines[int(line)])
    sam_file.close()
    selected_file.close()


if rand:
    random_select()

print("Done!")
