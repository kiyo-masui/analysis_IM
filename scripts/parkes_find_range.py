import fnmatch
import os
import codecs
import numpy as np

def find_pattern(pattern,root_dir):
    matches = []
    for root, dirnames, filenames in os.walk(root_dir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))

    return matches

matches = find_pattern("**.sdfits","/mnt/raid-project/gmrt/ycli/parkes/RAWDATA_SDF")

#print matches

foo = open('datalist_2012_yc.txt', 'w')
#for item in matches:
#    foo.write("{}\n".format(item))
foo.write('\n'.join(matches))
foo.close()

#list = open('datalist_2012_yc.txt', 'r').read().splitlines()

#print list
#thematches=np.array(matches)

#print thematches.size

#np.savetxt('datalist.txt', thematches, delimiter=',')

