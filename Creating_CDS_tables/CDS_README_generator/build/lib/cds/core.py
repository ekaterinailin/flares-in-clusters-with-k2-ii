from astropy.io import ascii
from astropy.table import Table, MaskedColumn
from cds.CDSColumn import CDSColumn, ColumnFormatter
import numpy as np
import sys
import os
import re
from string import Template
from textwrap import wrap, fill
import datetime

"""
   G.Landais (CDS) 28 nov 2015
   Generate ReadME and CDS standardized tables

   Standardized table in CDS is a table in ASCII format 
   with blank aligned columns 
"""

MAX_SIZE_README_LINE = 80
MAX_COL_INTLIMIT = 10000000

class CDSException(Exception):
        def __init__(self, value):
                self.value = value
        def __str__(self):
                return repr(self.value)                    
                
class CDSTable():
        """Manage table
        """
        def __init__(self, table, name = None, description = None):
                """name: tablename in output
                   description: the table description
                """
                self.columns = []
                self.__bytebybytetemplate = None
                self.notes = None
                self.table = None
                self.nlines = None
                self.linewidth = None
                self.maskedcolumn = self.__getMaskedColumn
                pass

        def setByteByByteTemplate(self, name):
                """Set a User Byte-By-Byte template
                """
                self.__bytebybytetemplate = name

        def getByteByByteTemplate(self):
                """Get the Byte-By-Byte template
                """
                return self.__bytebybytetemplate
        
        def __writeTable(self, fo, fmtc):
                for rec in self.table:
                    fo.write(fmtc[0].str(rec[0]))
                    for i in range(1,len(rec)):
                        fo.write(" "+fmtc[i].str(rec[i]))
                    fo.write("\n")

        def makeCDSTable(self):
                """Make the standardized table in 
                   ASCII format with aligned columns.
                   return the astropy table (astropy.Table)
                """
                fmtc = []
                for column in self.columns:
                       fmtc.append(ColumnFormatter(column))

                fo = open(self.name,"w")
                self.__writeTable(fo, fmtc)
                fo.close()


        def __getMaskedColumn(self, column):
                return MaskedColumn(column)


        def initMetaColumn(self, nullvalue=None):
                """Assign units, Fortran format to columns
                   Return list of CDSColumn (self.columns)
                """
                if self.table is None:
                    return
                if self.columns: return self.columns

                self.columns = []
                self.linewidth = 0
                for column in self.table.columns:
                        mcol = self.maskedcolumn(self.table[column])
                        col = CDSColumn(mcol, nullvalue)
                        self.columns.append(col)
                        self.linewidth += col.size+1

                self.linewidth -= 1
                return self.columns
        

class CDSAstropyTable(CDSTable):
        """Manage astropy table
        """
        def __init__(self, table, name = None, description = None):
                if not isinstance(table, Table):
                        raise CDSException("input is not Astropy Table")

                CDSTable.__init__(self, table, description)
                self.maskedcolumn = self.__getMaskedNoneColumn

                self.__filename = "astropytable"
                self.table = table
                if name is None: 
                        self.name = self.__filename
                else: 
                        self.name = name
                self.description = description
                self.nlines = len(table)

        def __getMaskedNoneColumn(self, column):
            # set None value to null
            return MaskedColumn(column, mask=[col is None for col in column])

class CDSNumpyTable(CDSTable):
        """Manage numpy table
        """
        def __init__(self, table, name = None, description = None):
                if not isinstance(table, np.ndarray):
                        raise CDSException("input is not a Numpy Table")

                CDSTable.__init__(self, table, description)
                self.maskedcolumn = self.__getMaskedNoneColumn

                self.__filename = "numpytable"
                self.table = Table(table)
                if name is None: 
                        self.name = self.__filename
                else: 
                        self.name = name
                self.description = description
                self.nlines = len(table)

        def __getMaskedNoneColumn(self, column):
            # set None value to null
            return MaskedColumn(column, mask=[col is None for col in column])



class CDSFileTable(CDSTable):
        """Manage Table in TSV/CSV or Ascii aligned column file
        """
        def __init__(self, table, name = None, description = None, data_start=None):
                if not isinstance(table, str):
                        raise CDSException("input is not a Table name")
                CDSTable.__init__(self, table, description)
                self.__filename = table
                self.table = ascii.read(self.__filename, data_start=data_start)
                if name is None: 
                        self.name = self.__filename+".cds"
                else: 
                        self.name = name
                self.description = description
                self.nlines = len(self.table)


class CDSMRTTable(CDSTable):
    """
    Manage MRT file to add information in the ReadMe
    """
    def __init__(self, tablein, tableout=None, description = None, set_limit=False):
        """
        :param tableout: table in output
        :param tablein:  table in intput
        :param description: description
        :param set_limit: add limits to the byte-by-byte description
        """
        CDSTable.__init__(self, tableout)

        self.name = tableout

        self.__bbb = None
        self.notes = []
        self.__begin_data = -1

        self.__tablename = tablein
        if tableout is None:
            self.__tablename_cds = tablein.replace("mrt","").replace(".txt","")+".dat"
        self.__tablename_cds = tableout

        self.linewidth = -1
        self.nlines = -1
        self.columns = None  # not yet available

        self.description = description
        if self.description is None:
            self.description = ""

        self.__is_written = False
        self.__parseTable()

        if set_limit:
            try:
                self.__bytebybytes_statistics()
            except Exception as e:
                sys.stderr.write("(error) byte-by-bytes parsing : {0}\n".format(str(e)))

    def __parseTable(self):
        num = 0
        linewidth = "0"
        regbbb = re.compile("^\s*[0-9]*[ -]+([0-9]*)\s+[A-Z][0-9.]+\s+.*$")
        regsep = re.compile("^[-]+$")
        regnot = re.compile("^\s*Note *[(][0-9]+[)].*")
        regtit = re.compile("^\s*Table\s*:\s(.*)\s*$")

        status = 0
        fd = open(self.__tablename, "r")
        for line in fd:
            num += 1
            if status == 5:
                    continue

            if status == 0:
                mo = regtit.match(line)
                if mo:
                    if len(mo.group(1).strip()) >1:
                        self.description = mo.group(1)
                    continue

                if line.find("Byte-by-byte Description ") == 0:
                    status = 1
                    self.__bbb = ""
                else:
                    if line.find("========") != 0:
                        self.description += " "+line.strip()

                continue

            elif status == 1:
                # continue until byte-by-bytes begining
                mo = regbbb.match(line)
                if mo :
                    status = 2
                else:
                    self.__bbb += line
                    continue

            if status == 2:
                # byte-by-bytes
                mo = regsep.match(line)
                if mo:
                    status = 3
                mo = regbbb.match(line)
                if mo:
                    linewidth = mo.group(1)
                self.__bbb += line

            elif status == 3:
                # notes -
                mo = regnot.match(line);
                if mo:
                        self.notes.append(line)
                        continue
                elif len(self.notes) == 0:
                    status = 5
                    self.__begin_data = num
                    continue

                mo = regsep.match(line)
                if mo:
                    status = 4
                    continue

                self.notes[-1] = self.notes[-1]+line

            elif status == 4:
                # data
                status = 5
                self.__begin_data = num
                continue

        fd.close()

        #print (self.__bbb)
        self.nlines = int(num)-int(self.__begin_data)+1
        self.linewidth = int(linewidth)

    def __bytebybytes_statistics(self):
        self.table = ascii.read(self.__tablename)
        columns = self.initMetaColumn()
        regbbb = re.compile("^(\s*[\d \-]*\d\s+[A-Z][0-9.]*\s+[^\s]+\s+[^\s]+\s+)(.*$)$")

        ncol = 0
        bbb_out = []
        bbb_in = self.__bbb.split("\n")
        for line in bbb_in:
            line = line.replace('?=""', '?')

            mo = regbbb.match(line)
            if mo is None:
                bbb_out.append(line)
                continue

            if ncol >= len(columns):
                sys.stderr.write("(warning) byte-by-byte limit can't be interpreted")
                return
            if columns[ncol].min is not None and columns[ncol].max is not None:
                s = mo.group(1)+"["+str(columns[ncol].min)+"/"+str(columns[ncol].max)+"]"
                if mo.group(2)[0] != '?': s += " "
                s += mo.group(2).lstrip()
                bbb_out.append(s)
            else:
                bbb_out.append(line)
            ncol += 1

        if ncol == len(columns):
            self.__bbb = "\n".join(bbb_out)+"\n"


    def __writeTable(self):
        i = 0
        fd = open(self.__tablename, "r")
        fout = open(self.__tablename_cds, "w")

        for line in fd:
            i += 1
            if i < self.__begin_data:
                continue
            fout.write(line)

        fd.close()
        fout.close()
        self.__is_written = True

    def makeCDSTable(self):
        """Make the standardized table in
           ASCII format with aligned columns.
           return the astropy table (astropy.Table)
        """
        print(self.__is_written)
        if self.__is_written : return
        self.__writeTable()


    def setByteByByteTemplate(self, name):
        """Set a User Byte-By-Byte template
        """
        raise Exception("not available in MRT")


    def getByteByByteTemplate(self):
        """Get the Byte-By-Byte template
        """
        return os.path.dirname(sys.argv[0])+"/bytebybyte.template"

    def getByteByByte(self):
        return self.__bbb


class CDSTablesMaker():
        """Generate standardized tables and ReadMe
        """
        def __init__(self, out=None, debug=False):
                """out: the output file (default: stdout)
                   debug: True/False
                """
                self.__dir = os.path.dirname(os.path.realpath(__file__))
                self.__tables = []
                self.__ReadMetemplate =  self.__dir+"/ReadMe.template"
                
                if out != None:
                        sys.stdout = open(out, 'w')
                
                self.title = 'Title ?'
                self.author = '1st author ?'
                self.catalogue = ''
                self.date = 'Date ?'
                self.abstract = 'Abstract ?'
                self.authors = 'Authors ?'
                self.bibcode = ''
                self.keywords = ''
                self.ref = None
                self.__templatevalue = None

                if debug is True:
                        self.__debug = self.__debugfn
                else:
                        self.__debug = self.__nodebugfn

        def __debugfn(self, buf):
                sys.stderr.write("(debug)"+str(buf)+"\n")
        def __nodebugfn(self, buf):
                pass

        def addTable(self, table, name=None, description=None, nullvalue=None):
                """Add a Table, memorize the meta-data and generate the standardized CDS table
                   name: the name used in output
                   description: table description
                   nullvalue: the null value
                """
                self.__debug("add table")
                if isinstance(table, CDSTable):
                        cdstable = table
                elif isinstance(table, Table):
                        cdstable = CDSAstropyTable(table, name, description)
                elif isinstance(table, str):
                        cdstable = CDSFileTable(table, name, description)
                elif isinstance(table,np.ndarray):
                        cdstable = CDSNumpyTable(table, name, description)
                else: 
                        raise CDSException("type "+type(table)+" is not accepted (only String or astropy.Table)")

                self.__tables.append(cdstable)
                self.__debug("append table")
                cdstable.initMetaColumn(nullvalue)
                self.__debug("retrieve meta OK")

                return cdstable

        def writeCDSTables(self):
                """Write tables in ASCII format
                """
                for table in self.__tables:
                        table.makeCDSTable()
                        self.__debug("make CDS table "+table.name)

        def getTablesInfo(self):
                """get tables informations
                """
                info = []
                for table in self.__tables:
                        info.append({'name': table.name, 'lines':len(table.table), 'width': table.linewidth})
                return info
        
        def getTables(self):
                """Get the CDSTable list
                """
                return self.__tables
 
        def printTablesIndex(self, outBuffer=False):
                """Print the tables index
                   outBuffer: true to get buffer, else write on output (default: False)
                """
                sz = [14, 0, 8]
                for tab in self.__tables:
                        if len(tab.name) > sz[0]: 
                                sz[0] = len(tab.name)
                        l = len(str(tab.linewidth))
                        if l > sz[1]: 
                                sz[1] = l

                fmtt = "{0:"+str(sz[0])+"s} {1:"+str(sz[1])+"d} {2:>"+str(sz[2])+"s}  {3:s}"
                buff = fmtt.format("ReadMe", MAX_SIZE_README_LINE, ".", "this file")+"\n"
                for tab in self.__tables:
                        buff += fmtt.format(self.__strFmt(tab.name),
                                tab.linewidth, 
                                str(tab.nlines),
                                self.__strFmt(tab.description))
                        buff += "\n"

                if outBuffer: return buff
                sys.stdout.write(buff)

        def __strFmt(self,string):
                if string is None: return ""
                else: return string

        def __strFmtRa(self, column, fmtb, startb):
                ra = column.getSexaRA()
                n = startb
                nvalue = ""
                if column.hasNull: nvalue = "? "
                        
                buff  = fmtb.format(n, n+ra.RAh.size-1, "",
                                      ra.RAh.fformat, ra.RAh.unit, ra.RAh.name, nvalue+ra.RAh.description)
                n += ra.RAh.size+1
                buff += '\n'
                buff += fmtb.format(n, n+ra.RAm.size-1, "",
                                      ra.RAm.fformat, ra.RAm.unit, ra.RAm.name, nvalue+ra.RAm.description)
                n += ra.RAm.size+1
                buff += '\n'
                buff += fmtb.format(n, n+ra.RAs.size-1, "",
                                      ra.RAs.fformat, ra.RAs.unit, ra.RAs.name, nvalue+ra.RAs.description)
                return buff

        def __strFmtDe(self, column, fmtb, startb):
                de = column.getSexaDE()
                n = startb
                nvalue = ""
                if column.hasNull: nvalue = "? "
                buff  = fmtb.format(n, n+de.DEsign.size-1, "",
                                      de.DEsign.fformat, de.DEsign.unit, de.DEsign.name, nvalue+de.DEsign.description)
                n += de.DEsign.size
                buff += '\n'
                buff += fmtb.format(n, n+de.DEd.size-1, "",
                                      de.DEd.fformat, de.DEd.unit, de.DEd.name, nvalue+de.DEd.description)
                n += de.DEd.size+1
                buff += '\n'
                buff += fmtb.format(n, n+de.DEm.size-1, "",
                                      de.DEm.fformat, de.DEm.unit, de.DEm.name, nvalue+de.DEm.description)
                n += de.DEm.size+1
                buff += '\n'
                buff += fmtb.format(n, n+de.DEs.size-1, "",
                                      de.DEs.fformat, de.DEs.unit, de.DEs.name, nvalue+de.DEs.description)
                return buff

        def __splitLine(self, line, shift=0):
                """Split line 80 char
                 :param shift: add left blank
                """
                if shift > MAX_SIZE_README_LINE:
                        shift = 0
                return ("\n"+" "*shift).join(wrap(line, width=MAX_SIZE_README_LINE-shift))

        def __add_authors(self, line, shift=0):
                """
                Split the line containing the authors without separate given and surname
                :param line: authors list in a line
                :param shift: add left blank
                :return:
                """
                # Find all spaces followed by authors's initials
                space_letter_list = re.findall(" (?:[A-Z]\.-?)+",line)

                # Replace founded spaces by !
                for space_letter in space_letter_list:
                        line = line.replace(space_letter,"!" + space_letter.strip())

                # Wrap the text by using spaces as breakpoint and then replace ! by spaces so given name and surname
                # are not separate
                new_line = fill(line, width=MAX_SIZE_README_LINE, subsequent_indent=shift*" ").replace("!"," ")

                return new_line


        def __add_keywords(self, line, shift=0):
                """
                Split the line containing the authors without separate given and surname
                :param line: keywords list in a line
                :param shift: add left blank
                :return:
                """
                # Replace all spaces that are NOT precede by ; with !
                line = re.sub("(?<!;) ","!",line)

                # Wrap the text by using spaces as breakpoint and then replace ! by spaces so there is no break line
                # between two ;
                new_line = fill(line, width=MAX_SIZE_README_LINE, subsequent_indent=shift*" ").replace("!"," ")

                return new_line



        def printByteByByte(self, table, outBuffer=False):
                """Print byte-by-byte
                   table: the CDSTable
                   outBuffer: true to get buffer, else write on output (default: False)
                """
                if isinstance(table,CDSMRTTable):
                        buff = table.getByteByByte()
                        if table.notes != None:
                            buff += "-" * 80 + "\n"
                            for line in table.notes:
                                try:
                                    buff += line.strip().encode('ascii', 'replace')
                                except:
                                    sys.stderr.write("Error : unicode detected in Notes\n")
                                    buff += "?"*len(line)+"\n"
                                buff+="\n"
                            buff += "-" * 80 + "\n"

                        if outBuffer: return buff
                        sys.stdout.write(buff)
                        return


                columns = table.columns
                startb = 1
                sz = [0, 0, 1, 7]
                l = len(str(table.linewidth))
                if l>sz[0]:
                        sz[0] = l
                        sz[1] = l
                self.__debug("size sz="+str(sz))
                fmtb = "{0:"+str(sz[0])+"d}-{1:"+str(sz[1])+"d} {2:"+str(sz[2])+"s}"
                for column in columns:
                        if len(column.name) > sz[3]:
                                sz[3] = len(column.name)
                fmtb += " {3:6s} {4:6s} {5:"+str(sz[3])+"s} {6:s}"
                buff = ""
                nsplit = sz[0]+sz[1]+sz[2]+sz[3]+16

                for column in columns:
                        endb = column.size+startb-1
                        if column.ffmt[0] == 'R':
                                buff += self.__strFmtRa(column, fmtb, startb)+"\n"
                        elif column.ffmt[0] == 'D':
                                buff += self.__strFmtDe(column, fmtb, startb)+"\n"
                        else:
                                description = column.description
                                if column.hasNull:
                                        nullflag = "?"
                                else:
                                        nullflag = ""

                                if column.ffmt[0] in "IF":
                                        if abs(column.min) < MAX_COL_INTLIMIT and abs(column.max) < MAX_COL_INTLIMIT:
                                                if column.min == column.max:
                                                    description = "[{0:.1f}]{1} {2}".format(column.min, nullflag, description)
                                                else:
                                                    description = "[{0:.1f}/{1:.1f}]{2} {3}".format(column.min, column.max,
                                                                                     nullflag, description)

                                newline = fmtb.format(startb, endb, "",
                                              self.__strFmt(column.ffmt),
                                              self.__strFmt(column.unit),
                                              self.__strFmt(column.name),
                                              description)

                                if len(newline) > MAX_SIZE_README_LINE:
                                        buff += ("\n").join(wrap(newline,
                                                                 subsequent_indent=" "*nsplit,
                                                                 width=MAX_SIZE_README_LINE))
                                        buff += "\n"
                                else: buff += newline + "\n"
                        startb = endb+2
                
                if table.notes != None: 
                        buff+="-"*80+"\n"
                        for line in table.notes: buff += line+"\n"
                        buff +="-"*80+"\n"

                if outBuffer: return buff
                sys.stdout.write(buff)

        def __getByteByByteTemplate(self, table):
                templateValue={ 'file': table.name, 
                                'bytebybyte': self.printByteByByte(table, outBuffer=True)}

                templatename =  table.getByteByByteTemplate()
                if templatename is None: templatename = self.__dir+"/bytebybyte.template"
                filein = open( templatename )
                src = Template( filein.read() )
                return src.substitute(templateValue)
                
        def setReadmeTemplate(self, templatename, templateValue = None):
                """Set a user ReadMe template
                :param templateName: the template name
                :param templateValue: dictionary to fill added variable in templatename
                """
                self.__ReadMetemplate = templatename
                self.__templatevalue = templateValue

        def putRef(self, catname, title=""):
                """
                Put a reference.
                catname catalgogue name (string)
                title the title (string)
                """
                if self.ref is None: self.ref = []
                if catname is None: raise Exception("catname is required")
                self.ref.append((catname, title))

        def printRef(self, outBuffer):
                """
                See also in ReadMe
                outBuffer: true to get buffer, else write on output (default: False)
                """
                if self.ref is None or len(self.ref)==0 : return

                buf = ""
                # Find the highest string length in first references column
                max_len = len(max([i[0] for i in self.ref],key=len))

                # Write references and align the double dot symbols
                for ref in self.ref:
                        buf += self.__splitLine(" {0:<{max}} : {1}".format(ref[0],ref[1],max=max_len))+"\n"

                if outBuffer is True: return buf
                sys.stdout.write(buf)


        def makeReadMe(self, out=sys.stdout):
                """Print the ReadMe
                """
                templateValue={ 'catalogue': self.catalogue,
                                'title': self.__splitLine(self.title),
                                'author': self.author,
                                'date': self.date,
                                'abstract': self.__splitLine(self.abstract, shift=2),
                                'authors' : self.__add_authors(self.authors, shift=4),
                                'bibcode' : "="+self.bibcode,
                                'keywords' : self.__add_keywords(self.keywords, shift=len("keywords: ")),
                                'tablesIndex': self.printTablesIndex(outBuffer=True),
                                'seealso' : self.printRef(outBuffer=True),
                                'bytebybyte': '',
                                'today': datetime.datetime.now().strftime("%d-%b-%Y")}

                if self.__templatevalue is not None:
                    for key in self.__templatevalue: templateValue[key] =  self.__templatevalue[key]

                buff = ""
                for table in self.__tables:
                        buff += self.__getByteByByteTemplate(table)
                        buff += "\n"
                templateValue['bytebybyte'] = buff

                filein = open(self.__ReadMetemplate )
                src = Template( filein.read() )
                result = src.substitute(templateValue)
                out.write(result)


