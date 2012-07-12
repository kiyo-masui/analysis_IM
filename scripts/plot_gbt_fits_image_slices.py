# This script plots image slices from a fits file format (written for an early version of the mapmaker that had a fits output).
import pylab
import pyfits
import sys
from numpy import *

filename1 = sys.argv[1]
filename2 = sys.argv[2]

filename1 = filename1.split('.')[0]
filename2 = filename2.split('.')[0]

hdudata = pyfits.open(filename1+'.fits')
hdudata2 = pyfits.open(filename2+'.fits')
#0 is nothing, 1 is only one.

length = len(hdudata2[1].data)
print length

count = 0
count2 = 0

for i in range(0,length) :
    cfr = int(hdudata[1].data[i][12]/1000000)
    if cfr == 695 :
        if hdudata[1].data[i][15] == -5 :
            if hdudata[1].data[i][5] == 'F' :
                if count == 0 :
                    count += 1
                    image_cube = hdudata[1].data[i][0]
                    cfr2 = int(hdudata2[1].data[i][12]/1000000)
                    if cfr == cfr2 :
                        if hdudata2[1].data[i][15] == -5 :
                            if hdudata2[1].data[i][5] == 'F' :
                                if count2 == 0 :
                                    count2 += 1
			            image_cube_2 = hdudata2[1].data[i][0]
                                    freq = zeros(len(image_cube_2))
                                    for k in range(0,len(image_cube_2)) :
                                        freq[k] = 720.0 - k*0.023855
                                    pylab.plot(freq,image_cube,label='Before transition')
                                    pylab.plot(freq,image_cube_2, label='After transition')
                                    pylab.legend()
                                    pylab.savefig('comparison_3C48_XX_'+str(cfr)+'_.png')
                                    pylab.clf()      
count = 0
count2 = 0

for i in range(0,length) :
    cfr = int(hdudata[1].data[i][12]/1000000) 
    if cfr == 725 :
        if hdudata[1].data[i][15] == -5 : 
            if hdudata[1].data[i][5] == 'F' : 
                if count == 0 :
                    image_cube = hdudata[1].data[i][0]
                    for j in range(0,i) : 
                        cfr2 = int(hdudata2[1].data[j][12]/1000000)
                        if cfr == cfr2 :
                            if hdudata2[1].data[j][15] == -5 :
                                if hdudata2[1].data[j][5] == 'F' :
                                    if count2 == 0 : 
                                        count2+=1
                                        count += 1
                                        image_cube_2 = hdudata2[1].data[j][0]
                                        pylab.plot(image_cube,label='Before transition')
                                        pylab.plot(image_cube_2, label='After transition')
                                        pylab.legend()
                                        pylab.savefig('comparison_3C48_XX_'+str(cfr)+'_.png')
                                        pylab.clf()
count = 0 
count2 = 0
 
for i in range(0,length) : 
    cfr = int(hdudata[1].data[i][12]/1000000)  
    if cfr == 755 :
        if hdudata[1].data[i][15] == -5 :
            if hdudata[1].data[i][5] == 'F' :
                if count == 0 :
                    image_cube = hdudata[1].data[i][0]
                    for j in range(0,i) :
                        cfr2 = int(hdudata2[1].data[j][12]/1000000)
                        if cfr == cfr2 :
                            if hdudata2[1].data[j][15] == -5 :
                                if hdudata2[1].data[j][5] == 'F' :
                                    if count2 == 0 :
                                        count2+=1
                                        count += 1
                                        image_cube_2 = hdudata2[1].data[j][0]
                                        pylab.plot(image_cube,label='Before transition')
                                        pylab.plot(image_cube_2, label='After transition')
                                        pylab.legend()
                                        pylab.savefig('comparison_3C48_XX_'+str(cfr)+'_.png')
                                        pylab.clf()
count = 0 
count2 = 0
 
for i in range(0,length) : 
    cfr = int(hdudata[1].data[i][12]/1000000)  
    if cfr == 785 :
        if hdudata[1].data[i][15] == -5 :
            if hdudata[1].data[i][5] == 'F' :
                if count == 0 :
                    count += 1
                    image_cube = hdudata[1].data[i][0]
                    cfr2 = int(hdudata2[1].data[i][12]/1000000)
                    if cfr == cfr2 :
                        if hdudata2[1].data[i][15] == -5 :
                            if hdudata2[1].data[i][5] == 'F' :
                                if count2 == 0 :
                                    count2 += 1
                                    image_cube_2 = hdudata2[1].data[i][0]
                                    pylab.plot(image_cube,label='Before transition')
                                    pylab.plot(image_cube_2, label='After transition')
                                    pylab.legend()
                                    pylab.savefig('comparison_3C48_XX_'+str(cfr)+'_.png')
                                    pylab.clf()
count = 0 
count2 = 0
 
for i in range(0,length) : 
    cfr = int(hdudata[1].data[i][12]/1000000)  
    if cfr == 815 :
        if hdudata[1].data[i][15] == -5 :
            if hdudata[1].data[i][5] == 'F' :
                if count == 0 :
                    count += 1
                    image_cube = hdudata[1].data[i][0]
                    cfr2 = int(hdudata2[1].data[i][12]/1000000)
                    if cfr == cfr2 :
                        if hdudata2[1].data[i][15] == -5 :
                            if hdudata2[1].data[i][5] == 'F' :
                                if count2 == 0 :
                                    count2 += 1
                                    image_cube_2 = hdudata2[1].data[i][0]
                                    pylab.plot(image_cube,label='Before transition')
                                    pylab.plot(image_cube_2, label='After transition')
                                    pylab.legend()
                                    pylab.savefig('comparison_3C48_XX_'+str(cfr)+'_.png')
                                    pylab.clf()
count = 0 
count2 = 0
 
for i in range(0,length) : 
    cfr = int(hdudata[1].data[i][12]/1000000)  
    if cfr == 845 :
        if hdudata[1].data[i][15] == -5 :
            if hdudata[1].data[i][5] == 'F' :
                if count == 0 :
                    count += 1
                    image_cube = hdudata[1].data[i][0]
                    cfr2 = int(hdudata2[1].data[i][12]/1000000)
                    if cfr == cfr2 :
                        if hdudata2[1].data[i][15] == -5 :
                            if hdudata2[1].data[i][5] == 'F' :
                                if count2 == 0 :
                                    count2 += 1
                                    image_cube_2 = hdudata2[1].data[i][0]
                                    pylab.plot(image_cube,label='Before transition')
                                    pylab.plot(image_cube_2, label='After transition')
                                    pylab.legend()
                                    pylab.savefig('comparison_3C48_XX_'+str(cfr)+'_.png')
                                    pylab.clf()
count = 0 
count2 = 0
 
for i in range(0,length) : 
    cfr = int(hdudata[1].data[i][12]/1000000)  
    if cfr == 875 :
        if hdudata[1].data[i][15] == -5 :
            if hdudata[1].data[i][5] == 'F' :
                if count == 0 :
                    count += 1   
                    image_cube = hdudata[1].data[i][0]
                    cfr2 = int(hdudata2[1].data[i][12]/1000000)
                    if cfr == cfr2 :
                        if hdudata2[1].data[i][15] == -5 :
                            if hdudata2[1].data[i][5] == 'F' :
                                if count2 == 0 :
                                    count2 += 1
                                    image_cube_2 = hdudata2[1].data[i][0]
                                    pylab.plot(image_cube,label='Before transition')
                                    pylab.plot(image_cube_2, label='After transition')
                                    pylab.legend()
                                    pylab.savefig('comparison_3C48_XX_'+str(cfr)+'_.png')
                                    pylab.clf()
count = 0 
count2 = 0
 
for i in range(0,length) : 
    cfr = int(hdudata[1].data[i][12]/1000000)  
    if cfr == 905 :
        if hdudata[1].data[i][15] == -5 :
            if hdudata[1].data[i][5] == 'F' :
                if count == 0 :
                    count += 1
                    image_cube = hdudata[1].data[i][0]
                    cfr2 = int(hdudata2[1].data[i][12]/1000000)
                    if cfr == cfr2 :
                        if hdudata2[1].data[i][15] == -5 :
                            if hdudata2[1].data[i][5] == 'F' :
                                if count2 == 0 :
                                    count2 += 1
                                    image_cube_2 = hdudata2[1].data[i][0]
                                    pylab.plot(image_cube,label='Before transition')
                                    pylab.plot(image_cube_2, label='After transition')
                                    pylab.legend()
                                    pylab.savefig('comparison_3C48_XX_'+str(cfr)+'_.png')
                                    pylab.clf()

