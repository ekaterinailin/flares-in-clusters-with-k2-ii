#!/usr/bin/python


from cds import *
import numpy as np
from astropy.table import Table

if __name__ == "__main__":

    # Create the ReadMe/Table maker
    #tablemaker = CDSTablesMaker(out='ReadMe')
    tablemaker = CDSTablesMaker()

    # Add astropy Table
    data = Table({'a': [1, 2, None],
              'b': [4.0, 1115.0, 6.12340],
              'f': [1.1,2,None],
              's': ['un','trois','deux'],
              's2': ['un', None, 'deux']},
              names=['a', 'b','f','s','s2'])
    table = tablemaker.addTable(data, name="data.cds")
    #table.setByteByByteTemplate("/home/landais/workspace/pyreadme/userbytebybyte.template")
    # Print bytebybyte 
    #tablemaker.printByteByByte(table)

    # Add ASCII aligned columns (fortran) table
    table = tablemaker.addTable("h_dm_com.dat",
                        name="h_dm_com.cds",
                        description="hipparcos table",
                        nullvalue="--")

    # Add TSV table
    #table = tablemaker.addTable("asu.tsv", description="my tsv table")

    # Add TSV table (ignore the 10 first lines)
    #cdstable = CDSFileTable("asu.tsv", description="my tsv table")
    #table = tablemaker.addTable(cdstable)
    
    # Add TSV table (contaning sexagesimal columns)
    #table = tablemaker.addTable("asu_sexa.tsv", description="my tsv table")
    #table.columns[5].setSexaRa(precision=1)
    #table.columns[6].setSexaDe()
    
    # Add numpy table
    ntab = np.array([(1.1,1,'test'), (2.2,2,'test2'), (3.3,3,None)],
                dtype=[('mag', np.float), ('recno', np.int32), ('comment', np.str_, 10)])
    #table = tablemaker.addTable(ntab,
    #                    name="myNumpyTable",
    #                    description="test numy table")


    # Print table index (all tables)
    #tablemaker.printTablesIndex()
    tablemaker.writeCDSTables()

    # Customize ReadMe output
    tablemaker.title = "my title"
    tablemaker.author = 'G.Landais'
    tablemaker.date = 2015
    tablemaker.abstract = "This is my abstract"
    tablemaker.putRef("II/246", "2mass")
    tablemaker.putRef("II/220", "2220 catalogue")

    # Print ReadMe
    tablemaker.makeReadMe()

    table = CDSFileTable("myNumpyTable")




