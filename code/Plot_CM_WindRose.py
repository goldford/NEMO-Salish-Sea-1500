import matplotlib.pyplot as plt
from windrose import WindroseAxes
import numpy as np

import pickle
import os
import csv

cm_p = 'D:/temp_nemo/RUN216/PROCESS/'
cm_f = 'CM_class4_SalishSea1500-RUN216_depth-below-surface.pickle'

out_p = '../data/evaluation/'
out_f = 'cm_locations.csv'

cm_class4 = pickle.load(open(os.path.join(cm_p,cm_f), 'rb'))

