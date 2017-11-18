import numpy as np

X = np.array([[0, 1],
              [1, 0]])
Y = np.array([[0, -1],
              [1, 0]])
Z = np.array([[1, 0],
              [0, -1]])
C = np.array([[1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 0, 1],
              [0, 0, 1, 0]])
P_list = [np.eye(2), X, Y, Z]

def get_cnot_err_table():
    error_list = []
    for key1, P1 in enumerate(P_list):
        for key2, P2 in enumerate(P_list):
            error_list.append((np.kron(P1, P2), key1, key2))

    cerr_table = np.zeros((4,4,2))
    for Ei, i1, j1 in error_list:
        Ef = np.dot(np.dot(C,Ei), C)
        for E,i2, j2 in error_list:
            if np.sum(np.absolute(Ef -E)) <=1e-5:
                # print(np.sum(np.absolute(Ef -E)))
                cerr_table[i1][j1][0] = i2
                cerr_table[i1][j1][1] = j2
                # print(i1, j1, ' -> ', i2, j2, '\n')
    # print(cerr_table)
    return cerr_table.astype(int)



def get_prod_table():
    product_table = np.zeros((4,4))
    for i, P1 in enumerate(P_list):
        for j, P2 in enumerate(P_list):
            if i ==0: product_table[i][j] = j;
            elif j==0: product_table[i][j] = i;
            elif i == j: product_table[i][j] = 0;
            else:
                for k in (1,2,3):
                    if k!=i and k!=j: product_table[i][j] = k

    # print(product_table)
    return product_table.astype(int)




cnot_err_table = get_cnot_err_table()
# np.savetxt("cnot_err_table.txt", cnot_err_table, fmt ='%d')

product_table = get_prod_table()