################################################################################
# basicio.py
# Sumeet Gupta
# Created on Jan 2020
# class with a collection of general methods for file/data handling
################################################################################

import os

class basicio():

    # find files within a directory (does not search sub directories)
    def findfiles(self, inputdir, suffix):
        return [os.path.join(inputdir, x) for x in os.listdir(inputdir) if x.endswith(suffix)]
    
    def readfile(self, filename):
        with open (filename, "r", encoding='utf-8') as fh:
            content = fh.read()
        return content

    # find files within a directory and it's subdirectories
    def getallfiles(self, inputdir, suffix):
        filelist = []
        [[filelist.append(os.path.join(i[0], j)) for j in i[2] if j.endswith(suffix)] for i in os.walk(inputdir)]
        return filelist

    # reformat the STAR and Qualimap objects for pandas dataframe
    def reformatinput(self, data):
        h = ['File']
        for sn in sorted(data.keys()):
            for k in data[sn].keys():
                if type(data[sn][k]) is not dict and k not in h:
                    h.append(str(k))
        rows = [h]

        for sn in sorted(data.keys()):
            rows.append([str(sn)] + [ str(data[sn].get(k, 'NA')) for k in h[1:] ])
        
        return rows