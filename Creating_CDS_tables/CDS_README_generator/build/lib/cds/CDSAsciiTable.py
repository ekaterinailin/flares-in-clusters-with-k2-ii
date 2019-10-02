
class CDSAsciiColumn():
        def __init__(self, begin, end):
                self.begin = begin
                self.end = end

class CDSAsciiTable():
        def __init__(self, filename):
                fd = open(filename)
                lines = fd.readlines()
                fd.close()
                self.lines = len(lines)
                self.size = 0

                for line in lines:
                        n = len(line)
                        if n> self.size: self.size = n
                blanks=[False]*self.size             
                for line in lines:
                        for i in range(len(line)):
                                if line[i] == ' ':
                                        blanks[i] = True
                 
                print(blank)
