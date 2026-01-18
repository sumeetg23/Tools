################################################################################
# parse.py
# Sumeet Gupta
# Created on Jan 2020
# Looks for the final log files from the STAR aligner and output a summary excel
# Looks for the final summary files from Qualimap and output a summary excel
################################################################################


import os
from STAR import staralignment
import sys
import numpy as np
import pandas as pd
from qualimap import qualimap

def main():

    # input directory
    inputdir = sys.argv[1]

    # generate object
    parsedcontent = staralignment(inputdir)
    qualicontent = qualimap(inputdir)
 
    # convert to dataframe
    stardata = pd.DataFrame(parsedcontent.reformatfornp())
    qualidata = pd.DataFrame(qualicontent.reformatfornp())

    # set indexes
    stardata.columns = stardata.iloc[0]
    stardata = stardata[1:]
    stardata.set_index("File", inplace=True)
    
    qualidata.columns = qualidata.iloc[0]
    qualidata = qualidata[1:]
    qualidata.set_index("File", inplace=True)

    # write to excel
    stardata.to_excel("Summary-STAR.xlsx") 
    qualidata.to_excel("Summary-Qualimap.xlsx") 


if __name__ == "__main__":
    
    main()