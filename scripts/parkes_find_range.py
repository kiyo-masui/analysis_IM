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

matches = find_pattern("*east1*.sdfits","/mnt/raid-project/gmrt/raid-pen/pen/Parkes/2dF/DATA/p641/sdfits/")

foo = open('datalist.txt', 'w')
#for item in matches:
#    foo.write("{}\n".format(item))
foo.write('\n'.join(matches))
foo.close()

#thematches=np.array(matches)

#print thematches.size

#np.savetxt('datalist.txt', thematches, delimiter=',')

