import csv
import pandas as pd
import random
import numpy as np

m = 0
T = []
M = []
D = []
with open('iip_98_2000.txt', 'r') as fp, open('iip.cnt', 'w', encoding='UTF-8') as fp_cnt, open('iip.data', 'w', encoding='UTF-8') as fp_data:
    lines = fp.readlines()
    # print(lines[0].split())
    for line in lines:
        fp_cnt.write(str(10) + ' ')
        line = line.split()
        if (len(line) != 20): continue
        m += 1
        type = line[6]
        T.append(type)
        melt = float(line[16])
        M.append(melt)
        days = float(line[17])
        D.append(days)
    print(m)

    maxM = max(M)
    minM = min(M)
    maxD = max(D)
    minD = min(D)

    for i in range(len(T)):
        c = 0
        if T[i] == 'R/V':
            c = 8
        elif T[i] == 'VIS':
            c = 7
        elif T[i] == 'RAD':
            c = 6
        for j in range(0, 10, 1):
            if j < c:
                fp_data.write(str(M[i]) + ' ')
                fp_data.write(str(D[i]) + '\n')
            else:
                fp_data.write(str(2*maxM) + ' ')
                fp_data.write(str(2*maxD) + '\n')

