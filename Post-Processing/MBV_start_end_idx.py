# filename: MBV_start_end_idx.py

import numpy as np

def MBV_start_end_idx(MBV):

    MBV_idx = np.where(MBV == True)[0]
    print(np.where(MBV == True))
    start_idx = MBV_idx[0]
    end_idx = np.where((MBV == False) & (np.arange(len(MBV)) > start_idx))[0][0]
    return start_idx, end_idx
