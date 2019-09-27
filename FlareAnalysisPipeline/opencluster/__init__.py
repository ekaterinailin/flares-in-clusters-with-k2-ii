import os

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

#This I found here: https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
from pandas import options
options.mode.chained_assignment = None  # default='warn'

import logging
import datetime
# create logger
logdate = datetime.datetime.now().strftime("%Y_%m_%d")
logname = 'opencluster_{}'.format(logdate)
logger = logging.getLogger(logname)
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('ancillary/{}.log'.format(logname))
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

logger.info('Started logger.')
