#!/usr/bin/python3
# Need to install the Python3 anafile wrapper that can be downloaded here: http://cds.u-strasbg.fr/resources/doku.php?id=anafile

from cds.core import *
import numpy as np

if __name__ == "__main__":

    tablemaker = CDSTablesMaker(out='ReadMe')

    # Add TSV table
    table = tablemaker.addTable("stars.csv", description="stellar parameters")
    #table = tablemaker.addTable("flares.csv", description="stellar flares")
    tablemaker.writeCDSTables()

    # Customize ReadMe output
    tablemaker.title = "Flares in Clusters. II. M35, Hyades, and Ruprecht 147"
    tablemaker.author = 'Ilin+'
    tablemaker.date = 2019
    tablemaker.abstract = "HERE IS THE PLACE FOR THE ACCEPTED PAPER ABSTRACT"

    # Print ReadMe
    with open('ReadMe', 'w') as out:
            tablemaker.makeReadMe(out=out)
