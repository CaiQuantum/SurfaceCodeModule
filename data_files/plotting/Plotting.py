import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# filename = '../CumulativeErrorData.txt'
filename = '../FullCircuitErrorCumulativeErrorData.txt'
data_array = np.genfromtxt(filename, delimiter=',')
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
with open(filename, 'w+') as file:
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
# ax.set_xlim([0.09,0.12])
ax.legend()
fig.tight_layout()
plt.show(block = False)
