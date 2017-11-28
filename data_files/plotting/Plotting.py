import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from re import match

os.chdir('..')
# output_file = 'CumulativeErrorData.txt'
output_file = 'CumulativeFullAsymCircuitErrorData.txt'
if os.path.exists(output_file):
    file_mode = 'a' # append if already exists
else:
    file_mode = 'w' # make a new file if not

file_list = [f for f in os.listdir('.') if match('CumulativeFullAsymCircuitErrorData[0-9]+\.txt', f)]
with open(output_file, file_mode) as outfile:
    for fname in file_list:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
        os.remove(fname)
data_array = np.genfromtxt(output_file, delimiter=',')

df = pd.DataFrame(data_array, columns = ['data error rate', 'code size', 'logical error rate', 'n_runs'])
print(df)
redundant_row = []
for index, row in df.iterrows():
    if index not in redundant_row:
        sub_df = df[df['data error rate'] == row['data error rate']][df['code size'] == row['code size']]
        if sub_df.shape[0] != 1:
            total_error = sum(sub_df['logical error rate']*sub_df['n_runs'])
            total_runs = sum(sub_df['n_runs'])
            index_list = list(sub_df.index)
            df.loc[index_list[0],'logical error rate'] = total_error/total_runs
            df.loc[index_list[0],'n_runs'] = total_runs
            redundant_row +=index_list[1:]

df.drop(df.index[redundant_row], inplace=True)
df.sort_values(['data error rate', 'code size'], inplace = True)
df.reset_index()
print(df)
with open(output_file, 'w+') as file:
    for index, row in df.iterrows():
        for i in range(4):
            file.write(str(row[i]))
            if i != 3:
                file.write(',')
        file.write('\n')

df = df.pivot(index='data error rate', columns='code size', values='logical error rate')
df.index.name = 'data error rate'
df.columns.name = 'code size'
df.columns = df.columns.map(int)
print(df)

df = df[[12,16,20,24]]
size_list = df.columns
error_list= df.index
fig, ax = plt.subplots(1)
for size in size_list:
    ax.plot(error_list, df[size], label=size)
ax.set_xlabel('data error rate')
ax.set_ylabel('logical error rate')
# ax.set_xlim([0.001,0.005])
ax.legend()
fig.tight_layout()
plt.show(block = False)
os.chdir('./plotting')