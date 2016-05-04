## ==================================================================
## 	     Producing custom colormaps
## ==================================================================
## warrenj 20160428 Producing custom colormaps from using
## http://jdherman.github.io/colormap/ and writing them to colormaps.py
import re
import numpy as np
from matplotlib import colors

name = "velcmap"






with open('/home/warrenj/VIMOS_project/analysis/colortable','r') as f:
    cm = f.read()
if ',' in cm:
    cm = np.array(map(float,re.split(',|;\n',cm)))
else:
    cm = np.array(map(float,cm.split()))

index = np.arange(0, len(cm), 3, dtype=int)
red = np.round(cm[index]/255.0,3)
green = np.round(cm[index+1]/255.0,3)
blue = np.round(cm[index+2]/255.0,3)

cm_file = "/home/warrenj/VIMOS_project/analysis/colormaps.py"
interval = 1.0/(len(index)-1)

with open(cm_file, 'a') as f:
    f.write("_cm = {'red':(")
    for i in range(len(index)-1):
        f.write('(%s,%s,%s),' % (str(i*interval),str(red[i]),
                                     str(red[i]))+'\n              ')
    f.write('(%s,%s,%s)),' % (str((i+1)*interval),str(red[i+1]),
                              str(red[i+1]))+'\n     ')
    
    f.write("'green':(")
    for i in range(len(index)-1):
        f.write('(%s,%s,%s),' % (str(i*interval),str(green[i]),
                                     str(green[i]))+'\n              ')
    f.write('(%s,%s,%s)),' % (str((i+1)*interval),str(green[i+1]),
                              str(green[i+1]))+'\n      ')

    f.write("'blue':(")
    for i in range(len(index)-1):
        f.write('(%s,%s,%s),' % (str(i*interval),str(blue[i]),
                                     str(blue[i]))+'\n              ')
    f.write('(%s,%s,%s))' % (str((i+1)*interval),str(blue[i+1]),
                              str(blue[i+1]))+'\n              ')

    f.write('} \n \n')

    f.write("%s = colors.LinearSegmentedColormap('%s', _cm) \n\n" % (name,name))

    
