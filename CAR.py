import csv
import pandas as pd
import random
import numpy as np

cols = ['model', 'brand', 'price', 'powerPS', 'kilometer', 'yearOfRegistration']
df = pd.read_csv('autos.csv')[cols]
df = df.dropna()
df = df[df['price'] != 0]
df = df[df['price'] < 999999]
df = df[df['powerPS'] > 0]
df = df[df['powerPS'] < 1000]
df = df[df['yearOfRegistration'] < 2023]
max_list = df.max().to_list()
min_list = df.min().to_list()
#cm = df.corr()
#print(cm)
df = df.drop_duplicates(subset=['price', 'powerPS', 'kilometer', 'yearOfRegistration'])
cars = df.groupby(['model', 'brand'])
m = 0
with open('car.cnt', 'w', encoding='UTF-8') as fp_cnt, open('car.data', 'w', encoding='UTF-8') as fp:
    for car in cars:
        m += 1
        fp_cnt.write(str(len(car[1])) + ' ')
        data = np.array(car[1]).tolist()
        for tuple in data:
            fp.write(str(tuple[2]/100) + ' ')
            fp.write(str((max_list[3] - tuple[3])/10) + ' ')
            fp.write(str(tuple[4]/10000) + ' ')
            fp.write(str(max_list[5] - tuple[5]) + '\n')
print(m, len(df))

