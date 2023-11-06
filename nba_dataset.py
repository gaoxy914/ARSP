from email import header
from operator import index
from nba_api.stats.static import players
from nba_api.stats.endpoints import playercareerstats
from nba_api.stats.endpoints import commonplayerinfo
import pandas as pd
import numpy as np
import os

nba_players = players.get_players()
print('Number of players fetched: {}'.format(len(nba_players)))
for player in nba_players:
    career = playercareerstats.PlayerCareerStats(player_id=player['id'])
    df = career.get_data_frames()[0]
    if not os.path.exists('player.csv'):
        df.to_csv('player.csv', mode='a', index=False, index_label=False)
    else:
        df.to_csv('player.csv', mode='a', index=False, index_label=False, header=False)

    #df.to_csv('player.csv', mode='a')
    #命中率, 得分, 助攻, 篮板, 抢断，盖帽
    #record = df[['PLAYER_ID', 'FG_PCT', 'PTS', 'AST', 'REB', 'STL', 'BLK']]
    #print(record)
    #record.to_csv('player.csv', index=False, mode='a', header=False)

cols = ['PLAYER_ID', 'PTS', 'AST', 'REB', 'STL', 'BLK']
df = pd.read_csv('player.csv')[cols]
df = df.dropna()
max_list = df.max().to_list()
min_list = df.min().to_list()
print(max_list)
print(min_list)
delta = []
for i in range(6):
    delta.append(max_list[i] - min_list[i])
players = df.groupby(['PLAYER_ID'])
with open('nba.cnt', 'w', encoding='UTF-8') as fp_cnt, open('nba.data', 'w', encoding='UTF-8') as fp_data:
    print(len(players))
    for player in players:
        fp_cnt.write(str(len(player[1])) + ' ')
        data = np.array(player[1]).tolist()
        for tuple in data:
            fp_data.write(str(round(1 - (tuple[1] - min_list[1])/delta[1], 3)) + ' ')
            fp_data.write(str(round(1 - (tuple[2] - min_list[2])/delta[2], 3)) + ' ')
            fp_data.write(str(round(1 - (tuple[3] - min_list[3])/delta[3], 3)) + ' ')
            fp_data.write(str(round(1 - (tuple[4] - min_list[4])/delta[4], 3)) + ' ')
            fp_data.write(str(round(1 - (tuple[5] - min_list[5])/delta[5], 3)) + '\n')
            #fp_data.write(str(tuple[2]) + ' ')

'''
#index = [1620, 982, 2516, 369, 3374]
index = [369, 1620, 982, 349, 1119, 2516, 545, 181, 99, 2579]
cols = ['PLAYER_ID', 'PTS', 'AST', 'REB', 'STL', 'BLK']
df = pd.read_csv('player.csv')[cols]
df = df.dropna()
max_list = df.max().to_list()
min_list = df.min().to_list()
print(max_list)
print(min_list)
delta = []
for i in range(6):
    delta.append(max_list[i] - min_list[i])
players = df.groupby(['PLAYER_ID'])
cnt = 0
for player in players:
    if cnt in index:
        id = player[0]
        player_info = commonplayerinfo.CommonPlayerInfo(player_id=id)
        print(cnt)
        print(player_info.get_data_frames())
    cnt += 1
'''



