from email import header
from operator import index
from nba_api.stats.static import players
from nba_api.stats.endpoints import playercareerstats
from nba_api.stats.endpoints import commonplayerinfo
from nba_api.stats.endpoints import playergamelog
from nba_api.stats.library.parameters import SeasonAll
import pandas as pd
import numpy as np
import os
import random


nba_players = players.get_players()
print('Number of players fetched: {}'.format(len(nba_players)))
for player in nba_players:
    #career = playercareerstats.PlayerCareerStats(player_id=player['id'])
    #df = career.get_data_frames()[0]
    #if not os.path.exists('player.csv'):
    #    df.to_csv('player.csv', mode='a', index=False, index_label=False)
    #else:
    #    df.to_csv('player.csv', mode='a', index=False, index_label=False, header=False)
    games_all = playergamelog.PlayerGameLog(player_id=player['id'], season=SeasonAll.all)
    df = games_all.get_data_frames()[0]
    if not os.path.exists('player_game_by_game.csv'):
        df.to_csv('player_game_by_game.csv', mode='a', index=False, index_label=False)
    else:
        df.to_csv('player_game_by_game.csv', mode='a', index=False, index_label=False, header=False)



cols = ['Player_ID', 'PTS', 'AST', 'STL', 'BLK', 'TOV', 'REB', 'MIN', 'FGM', 'GAME_DATE']
#cols = ['PLAYER_ID', 'PTS', 'REB', 'FGM', 'MIN', 'TOV', 'BLK', 'STL', 'AST']
df = pd.read_csv('player_game_by_game.csv')[cols]
#df = pd.read_csv('player.csv')[cols]
df = df.dropna()
df = df[pd.to_datetime(df['GAME_DATE'], format="%b %d, %Y") > pd.to_datetime('2003-01-01')]
df = df[pd.to_datetime(df['GAME_DATE'], format="%b %d, %Y") < pd.to_datetime('2023-01-01')]
max_list = df.max().to_list()
min_list = df.min().to_list()
print(max_list)
print(min_list)
df = df.drop_duplicates(['PTS', 'AST', 'STL', 'BLK', 'TOV', 'REB', 'MIN', 'FGM'])
players = df.groupby(['Player_ID'])
cm = df.corr()
cm.to_csv('corr_matrix.csv')
#players = df.groupby(['PLAYER_ID'])
m = 0
with open('nba.cnt', 'w', encoding='UTF-8') as fp_cnt, open('nba.data', 'w', encoding='UTF-8') as fp:
    for player in players:
        m += 1
        fp_cnt.write(str(len(player[1])) + ' ')
        data = np.array(player[1]).tolist()
        for tuple in data:
            fp.write(str(max_list[1] - tuple[1]) + ' ')
            fp.write(str(max_list[2] - tuple[2]) + ' ')
            fp.write(str(max_list[3] - tuple[3]) + ' ')
            fp.write(str(max_list[4] - tuple[4]) + ' ')
            fp.write(str(tuple[5]) + ' ')
            fp.write(str(max_list[6] - tuple[6]) + ' ')
            fp.write(str(max_list[7] - tuple[7]) + ' ')
            fp.write(str(max_list[8] - tuple[8]) + ' ')
            fp.write(str(tuple[0]) + '\n')
        
print(m, len(df))


