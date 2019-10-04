import pandas as pd
import os

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
NORM = pd.read_csv('{}/table_u0_g_col.txt'.format(PACKAGEDIR))
NORMG = pd.read_csv('{}/table_u0_g.txt'.format(PACKAGEDIR))
NORM.columns = NORM.columns.str.strip()
NORMG.columns = NORMG.columns.str.strip()
