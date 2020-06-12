#!/usr/bin/python3
# Need to install the Python3 anafile wrapper that can be downloaded here: http://cds.u-strasbg.fr/resources/doku.php?id=anafile

from cds.core import *
import numpy as np

if __name__ == "__main__":

    tablemaker = CDSTablesMaker(out='ReadMe')

    # Add TSV table
    table1 = tablemaker.addTable("stars.csv", description="stellar parameters")
    table2 = tablemaker.addTable("flares.csv", description="stellar flares")
    table3 = tablemaker.addTable("ffds.csv", description="flare frequency distributions and power law fit parameters")
    tablemaker.writeCDSTables()

    # Customize ReadMe output
    tablemaker.title = "Flares in Clusters. II. Pleiades, Hyades, Praesepe, Ruprecht 147, and M67"
    tablemaker.author = 'Ilin+'
    tablemaker.date = 2020
    tablemaker.abstract = "HERE IS THE PLACE FOR THE ACCEPTED PAPER ABSTRACT"

    # Print ReadMe
    with open('ReadMe', 'w') as out:
            tablemaker.makeReadMe(out=out)
