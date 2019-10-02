import numpy
import re
import sys
from astropy.table import Column, MaskedColumn

class SexaMeta:
        """Sexa metadata"""
        def __init__(self, name, description, fformat, unit):
                self.name = name
                self.description = description
                self.fformat = fformat
                self.unit = unit
                self.size = int(fformat[1])

class SexaRA():
        def __init__(self, precision):
                self.RAh = SexaMeta("RAh", "Right ascension (hour)", "I2" , "h")
                self.RAm = SexaMeta("RAm", "Right ascension (minute)", "I2" , "min")
                self.RAs = SexaMeta("RAs", "Right ascension (seconds)", "I2" , "s")
                if precision != 0:
                        self.RAs.fformat = "F"+str(3+precision)+"."+str(precision)
                        self.RAs.size = (3+precision)

class SexaDE():
        def __init__(self, precision):
                self.DEsign = SexaMeta("DE-", "Declination (degree)", "A1" , "---")
                self.DEd = SexaMeta("DEd", "Declination (degree)", "I2" , "deg")
                self.DEm = SexaMeta("DEm", "Declination (minute)", "I2" , "arcmin")
                self.DEs = SexaMeta("DEs", "Declinationion (seconds)", "I2", "arcsec")
                if precision != 0:
                        self.DEs.fformat = "F"+str(3+precision)+"."+str(precision)
                        self.DEs.size = (3+precision)

class CDSColumn():
        def __init__(self, column, nullvalue=None):
                if not isinstance(column, Column):
                        raise Exception("column must be astropy.table.Column")

                self.__regfloat = re.compile("([+-]*)([^eE.]+)([.]*)([0-9]*)([eE]*-*)[0-9]*")
                self.__sexa = [None, None]


                self.name = column.name
                self.__dbname = column.name
                if column.description:
                        self.description = column.description
                else: 
                        self.description = "Description"
                if column.unit is not None:
                        self.unit = column.unit
                else:
                        self.unit = "---"
                self.ffmt = None # Fortran format
                self.fmt = None  # printf format 
                self.size = None
                self.nullvalue = nullvalue
                self.hasNull = self.__hasNull(column)
                self.min = None
                self.max = None

                type = self.__getType(column)
                if type is int:
                        self.__getIntegerFFormat(column)
                elif type is float:
                        self.__getFloatFFormat(column)
                else:
                        self.__getStringFFormat(column)


        def __getType(self, column):
                if column.dtype.name.startswith("i"):
                        return int
                elif column.dtype.name.startswith("f"):
                        return float
                elif column.dtype.name.startswith("s"):
                        return str

                # if contains null values dtype is 'object'
                if self.hasNull is False:
                        return str

                for value in column:
                        if value is None: continue
                        if isinstance(value, (int, numpy.integer)):
                                return int
                        elif numpy.isreal(value):
                                return float

                        else:
                                return str
                return str

        def getName(self):
                """ get the name in input 
                    name != self.column which can be changed 
                """
                return self.__dbname

        def isSexa(self):
                """Return True if the column is in sexagesimal"""
                if self.__sexa[0] is None and self.__sexa[1] is None:
                        return False
                return True
        def isSexaRA(self):
                """Return True if Right ascension"""
                if self.__sexa[0] is None: return False
                return True
        def isSexaDE(self):
                """Return True if Declination"""
                if self.__sexa[1] is None: return False
                return True

        def setSexaRa(self, precision = 0):
                """set column as a Sexagesimal Right Ascension
                   precision: number of seconds decimals
                """
                if precision != 0:
                        if self.size != 9+precision : 
                                raise Exception("bad sexa format or bad precision (format: hh mm ss[.ss])")
                
                if self.ffmt[0] != 'A' : raise Exception("bad sexa format")
                self.ffmt = self.ffmt.replace('A', 'R')
                self.__sexa[0] = SexaRA(precision)

        def setSexaDe(self, precision = 0):
                """set column as a Sexagesimal Declination
                   precision: number of seconds decimals
                """
                if precision != 0:
                        if self.size < 9+precision or self.size > 10+precision: 
                                raise Exception("bad sexa format or bad precision (format: [+-]dd mm ss[.ss])")
                
                if self.ffmt[0] != 'A' : raise Exception("bad sexa format")
                self.ffmt = self.ffmt.replace('A', 'D')
                self.__sexa[1] = SexaDE(precision)

        def getSexaRA(self):
                """get Sexagesimal  Right ascension
                """
                return self.__sexa[0];

        def getSexaDE(self):
                """get Sexagesimal DEclination
                """
                return self.__sexa[1];
                
        def __hasNull(self, t):
                """Test null values in a numpy array
                """
                try:
                        n = len(t[t.argmin(fill_value=0)])
                except :
                        return True
                return False

        def __getStringFFormat(self, column):
                reg = re.compile("^[^0-9]+(\d+)$")
                mo = reg.match(column.dtype.str)
                if mo is None:
                        mo = reg.match(column.dtype.str)
                        if mo is None:
                                self.size = 50
                        else:
                                self.size = int(mo.group(1))
                else:
                        self.size = int(mo.group(1))

                try:
                        if self.hasNull:
                                mcol = column#MaskedColumn(column, mask=[col is None for col in column])
                                mcol.fill_value = ""
                                coltmp = Column(mcol.filled(), dtype=str)
                                self.size = int(re.sub(r'^[^0-9]+(\d+)$', r'\1', coltmp.dtype.str))
                        else:
                                self.size = int(re.sub(r'^[^0-9]+(\d+)$', r'\1', column.dtype.str))
                except Exception as e:
                        sys.stderr.write("(error)" + str(e) + "\n")
                        self.size = 50

                self.ffmt = "A"+str(self.size)
                self.fmt = str(self.size)+"s"

        def __getIntegerFFormat(self, column):
                if self.hasNull:
                        mcol = column#MaskedColumn(column, mask=[col is None for col in column])
                        mcol.fill_value=-999999
                        self.max = max(mcol.filled())
                        mcol.fill_value = +999999
                        self.min = min(mcol.filled())
                else:
                        self.max = max(column)
                        self.min = min(column)

                self.size = len(str(self.max))
                self.ffmt = "I"+str(self.size)
                self.fmt =  ">"+self.ffmt[1:]

        def __splitFloatFormat(self, value, fmt):
                """ get float format
                    value (IN) the float value
                    fmt (out) [size, prec, dec, ent, sign]
                    return True has scientific notation
                """
                mo = self.__regfloat.match(value)

                if mo is None: raise Exception(value+" is not a float number")

                fmt[0] = len(value)
                if mo.group(1) != "": fmt[4] = True
                else: fmt[4] = False
                fmt[2] = len(mo.group(2))
                fmt[3] = len(mo.group(4))
                fmt[1] = fmt[2]+fmt[3]

                if mo.group(5) != "":
                        # scientific notation
                        return True
                return False


        def __getFloatFFormat(self, column):

                if self.hasNull:
                        mcol = column#MaskedColumn(column)#, mask=[col is None for col in column])
                        mcol.fill_value = -999999
                        self.max = max(mcol.filled())
                        mcol.fill_value = +999999
                        self.min = min(mcol.filled())
                else:
                        self.max = max(column)
                        self.min = min(column)

                maxSize = 1
                maxprec = 0
                maxDec = 0
                maxEnt = 1
                sign = False
                reg = re.compile("^[ -]*$")
                fformat = 'F'
                fmt = [0, 0, 0, 0, False]

                for rec in column:
                        # skip null values
                        if rec is None:
                                continue
                        s = str(rec)

                        if reg.match(s): continue
                        if s == self.nullvalue:
                                self.hasNull = True
                                continue

                        if self.__splitFloatFormat(s, fmt) is True:
                                if fformat == 'F':
                                        maxSize = 1
                                        maxprec = 0
                                        maxDec = 0
                                # scientific notation
                                fformat = 'E'
                        else:
                                if fformat == 'E': continue

                        if maxprec < fmt[1]: maxprec= fmt[1]
                        if maxDec < fmt[3] : maxDec = fmt[3]
                        if maxEnt < fmt[2] : maxEnt = fmt[2]
                        if maxSize < fmt[0] : maxSize = fmt[0]
                        if fmt[4] : sign = True

                if fformat == 'E':      
                        self.size = maxSize
                        if sign: self.size += 1
                        self.ffmt = fformat+str(self.size)+"."+str(maxprec)
                        self.fmt = str(self.size)+"."+str(maxDec)+"e"
                else:
                        self.size = maxEnt+maxDec+1
                        if sign: self.size += 1
                        self.ffmt = fformat+str(self.size)+"."+str(maxDec)
                        self.fmt = self.ffmt[1:]+"f"

        def __getUnit(self): 
                s = self.name.lower()
                if s.index("magnitude"): self.unit = "mag"
                elif s.index("[ (]days[ )]"): self.unit = "d"

class ColumnFormatter():
        """Format column with apropriate type (FORTRAN type)
        """
        def __init__(self, column):
                if not isinstance(column, CDSColumn):
                        raise Exception("wrong argument type for ColumnFormatter:"+str(column))
                self.__column = column

                self.__nullvalue = column.nullvalue
                self.__fmtstring = "{0:"+str(column.size)+"s}"

                if column.ffmt[0] == 'F':
                        self.__fmtfloat = "{0:"+column.fmt+"}"
                        if column.nullvalue is None:
                                self.__str = self.__strFloat
                        else:
                                self.__str = self.__strFloatN
                elif column.ffmt[0] == 'I':
                        self.__fmtint = "{0:"+column.fmt+"}"
                        if column.nullvalue is None:
                                self.__str = self.__intString
                        else:
                                self.__str = self.__intStringN
                elif  column.ffmt[0] == 'E':
                        self.__fmtfloat = "{0:"+column.fmt+"}"
                        if column.nullvalue is None:
                                self.__str = self.__strFloat
                        else:
                                self.__str = self.__strFloatN
                else:
                        if column.isSexa(): 
                                if column.nullvalue is None:
                                        self.__str = self.__strSexa
                                else:
                                        self.__str = self.__strSexaN
                        elif column.nullvalue is None:
                                self.__str = self.__strString
                        else:
                                self.__str = self.__strStringN


        def __intString(self, value):
                if value is None:
                        return self.__fmtint.format("")
                return self.__fmtint.format(value)
        def __intStringN(self, value):
                if value is None:
                        return self.__fmtint.format("")
                if str(value).find(self.__nullvalue) > -1:
                        return self.__fmtint.format(" ")
                return self.__fmtint.format(value)

        def __strString(self, value):
                if value is None:
                        return self.__fmtstring.format(" ")
                return self.__fmtstring.format(value)
        def __strStringN(self, value):
                if value is None:
                        return self.__fmtstring.format(" ")
                elif not isinstance(value, numpy.string_) and not isinstance(value, numpy.str_):
                        return self.__fmtstring.format(" ")
                return self.__fmtstring.format(value.replace(self.__nullvalue,""))

        def __strFloat(self, value):
                if value is None:
                        return self.__fmtstring.format("")
                return self.__fmtfloat.format(float(value))
        def __strFloatN(self, value):
                if value is None:
                        return self.__fmtstring.format("")
                if str(value).find(self.__nullvalue) > -1:
                        return self.__fmtstring.format(" ")
                return self.__fmtfloat.format(float(value))

        def __strSexa(self, value):
                if not isinstance(value, numpy.string_):
                        return self.__fmtstring.format(" ")
                return self.__fmtstring.format(str(value).replace(":", " "))
        def __strSexaN(self, value):
                if not isinstance(value, numpy.string_):
                        return self.__fmtstring.format(" ")
                if str(value).find(self.__nullvalue) > -1:
                        return self.__fmtstring.format(" ")
                return self.__fmtstring.format(str(value).replace(":", " "))

        def str(self, value):
                """Return the Format"""
                return self.__str(value)
