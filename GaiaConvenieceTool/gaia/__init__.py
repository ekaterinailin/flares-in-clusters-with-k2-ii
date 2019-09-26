import pandas as pd
import os

NORM = pd.read_csv('gaia/table_u0_g_col.txt')
NORMG = pd.read_csv('gaia/table_u0_g.txt')
NORM.columns = NORM.columns.str.strip()
NORMG.columns = NORMG.columns.str.strip()