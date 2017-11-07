import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

input_mode = 1
if input_mode == 0:
    df = pd.read_table("../ErrorData_100runs.txt", delim_whitespace= True, index_col = 0, header= 0)
else:
    data_array = np.genfromtxt('../CumulativeErrorData.txt', delimiter=',')
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

    print(df)
    df.drop(df.index[redundant_row], inplace=True)
    df.reset_index()
    print(df)
#         for index2, row2 in df.iterrows():
#             if row1['data error rate'] == row1['data error rate'] and row1['code size'] == row2['code size']:
#                 row1['logical error rate'] = row1['logical error rate']*row1['n_runs'] + row2['logical error rate']*row2['n_runs']
#                 row1['logical error rate'] /= (row1['n_runs'] + row2['n runs'])
#                 row1['n_runs'] = row1['n_runs'] + row2['n runs']
#                 redundant_row.append()
#     df = df.pivot(index='data error rate', columns='code size', values='logical error rate')
# df.index.name = 'data error rate'
# df.columns.name = 'code size'
# df.columns = df.columns.map(int)
# print(df)

# size_list = df.columns
# error_list= df.index
# fig, ax = plt.subplots(1)
# for size in size_list:
#     ax.plot(error_list, df[size], label=size)
# plt.show()
