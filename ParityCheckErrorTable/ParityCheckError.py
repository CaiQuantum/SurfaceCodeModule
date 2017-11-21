import itertools
import numpy as np
import os

from ProdCnotTableGenrators import cnot_err_table, product_table

'''
Here we introduce:
4 hadamard error for data qubits
1 preparation error for ancilla
4 PAIRS of cnot error, hence 8 errors
1 readout error
A total of 2*4^(8+4+1) = 2^(27) possibilities. (Actually possible to enumerate over all possibilities if wanted)
'''


def get_error_table(n_trials, error_prob):
    pauli_error_rate = [1-error_prob]+[error_prob/3]*3
    cnot_error_rate = [1-error_prob]+[error_prob/15]*15
    readout_error_rate = [1-error_prob, error_prob]

    result_table = np.zeros((2, 4, 4, 4, 4))

    for trial in range(n_trials):
        # prep error on ancilla and hadamard error on data qubit
        q = np.random.choice(len(pauli_error_rate), 5, p=pauli_error_rate)

        # impose cnot error one by one:
        cnot_error_list = np.random.choice(len(cnot_error_rate), 4, p=cnot_error_rate)
        for i in range(1, 5):
            # commute the current error through the PERFECT cnot
            q[i], q[0] = cnot_err_table[q[i]][q[0]]
            # adding error to cnot
            e0 = int(cnot_error_list[i - 1] / 4)
            ei = cnot_error_list[i - 1] % 4
            q[0] = product_table[q[0]][e0]
            q[i] = product_table[q[i]][ei]
        '''
        Ideally, when there is no error, we will have ancilla in |0> when the state is in even parity. Now, we add in the error,
        we then have:
        I: state:|0>,  correct
        X: state:|1>,  wrong
        Y: state:|1>,  wrong
        Z: output|0>,  correct
        then we add in error_prob chance of flipping the readout
        '''
        if q[0] in (1, 2):
            q[0] = 1
        else:
            q[0] = 0

        if np.random.choice(len(readout_error_rate), p=readout_error_rate): q[0] = not q[0]
        result_table[q[0]][q[1]][q[2]][q[3]][q[4]] += 1

    cumul_result_name = "./Data/cumulative_result%.4f.npy"%error_prob
    if os.path.exists(cumul_result_name):
        cumul_result = np.load(cumul_result_name)
        result_table +=cumul_result
    np.save(cumul_result_name, result_table)
    # print(result_table)
    # print(np.sum(result_table))
    result_table /= np.sum(result_table)

    filename = "./Data/error_table%.4f.txt"%error_prob
    file = open(filename, 'w')
    for result in range(2):
        for q1 in range(4):
            for q2 in range(4):
                for q3 in range(4):
                    for q4 in range(4):
                        file.write("%.14f %d %d %d %d %d\n"%(result_table[result][q1][q2][q3][q4], result, q1, q2, q3, q4))
    file.close()

def get_error_table_exact(error_prob):
    pauli_error_rate = [1-error_prob]+[error_prob/3]*3
    cnot_error_rate = [1-error_prob]+[error_prob/15]*15
    readout_error_rate = [1-error_prob, error_prob]

    result_table = np.zeros((2, 4, 4, 4, 4))
    '''
    h means preparation and hadamard error, 
    c means cnot error, 
    r means readout error
    '''
    for h in itertools.product(range(4), repeat=5):
        for c in itertools.product(range(16), repeat=4):
            for r in range(2):
                # prep error on ancilla and hadamard error on data qubit
                q = np.array(h)
                # impose cnot error one by one:
                for i in range(1, 5):
                    # commute the current error through the PERFECT cnot, NOTE the first bit is control!!!
                    q[i], q[0] = cnot_err_table[q[i]][q[0]]
                    # adding error to cnot
                    e0 = int(c[i - 1] / 4)
                    ei = c[i - 1] % 4
                    q[0] = product_table[q[0]][e0]
                    q[i] = product_table[q[i]][ei]

                '''
                Ideally, when there is no error, we will have ancilla in |0> when the state is in even parity. Now, we add in the error,
                we then have:
                I: state:|0>,  correct
                X: state:|1>,  wrong
                Y: state:|1>,  wrong
                Z: output|0>,  correct
                then we add in error_prob chance of flipping the readout
                '''
                if q[0] in (1, 2):
                    q[0] = 1
                else:
                    q[0] = 0

                if r:
                    q[0] = not q[0]

                prob = 1
                for i in h:
                    prob *= pauli_error_rate[i]
                for i in c:
                    prob *= cnot_error_rate[i]
                prob *= readout_error_rate[r]

                result_table[q[0]][q[1]][q[2]][q[3]][q[4]] += prob

    filename = "./Data/new_exact_error_table%.4f.txt"%error_prob
    file = open(filename, 'w')
    for result in range(2):
        for q1, q2, q3, q4 in itertools.product(range(4), repeat=4):
            file.write("%.14f %d %d %d %d %d\n"%(result_table[result][q1][q2][q3][q4], result, q1, q2, q3, q4))
    file.close()

def get_Fowler_error_table_exact(error_prob, stabiliser_type):
    pauli_error_rate = [1-error_prob]+[error_prob/3]*3
    cnot_error_rate = [1-error_prob]+[error_prob/15]*15
    readout_error_rate = [1-error_prob, error_prob]

    result_table = np.zeros((2, 4, 4, 4, 4))
    '''
    h means preparation and hadamard error,
    c means cnot error,
    r means readout error
    '''
    if stabiliser_type == 'X':
        h_repeat = 3
    else:
        h_repeat = 1
    for h in itertools.product(range(4), repeat=h_repeat):
        for c in itertools.product(range(16), repeat=4):
            for r in range(2):
                q = np.zeros(5)
                q = q.astype(int)
                # prep error on ancilla
                q[0] = h[0]
                if stabiliser_type == 'X':
                    #hadamard error on ancilla
                    q[0] = product_table[q[0]][h[1]]
                # impose cnot error one by one:
                for i in range(1, 5):
                    # commute the current error through the PERFECT cnot
                    if stabiliser_type == 'X':
                        q[0], q[i] = cnot_err_table[q[0]][q[i]]
                    else:
                        q[i], q[0] = cnot_err_table[q[i]][q[0]]
                    # adding error to cnot
                    e0 = int(c[i - 1] / 4)
                    ei = c[i - 1] % 4
                    q[0] = product_table[q[0]][e0]
                    q[i] = product_table[q[i]][ei]

                if stabiliser_type == 'X':
                    #pass the error through Hadamard.
                    if q[0] == 1:
                        q[0] = 3
                    elif q[0] == 3:
                        q[0] = 1
                    #hadamard error on ancilla
                    q[0] = product_table[q[0]][h[2]]

                '''
                Ideally, when there is no error, we will have ancilla in |0> when the state is in even parity. Now, we add in the error,
                we then have:
                I: state:|0>,  correct
                X: state:|1>,  wrong
                Y: state:|1>,  wrong
                Z: output|0>,  correct
                then we add in error_prob chance of flipping the readout
                '''
                if q[0] in (1, 2):
                    q[0] = 1
                else:
                    q[0] = 0

                if r:
                    q[0] = not q[0]

                prob = 1
                for i in h:
                    prob *= pauli_error_rate[i]
                for i in c:
                    prob *= cnot_error_rate[i]
                prob *= readout_error_rate[r]

                result_table[q[0]][q[1]][q[2]][q[3]][q[4]] += prob

    filename = "./Data/fowler_error_table_%s%.4f.txt"%(stabiliser_type, error_prob)
    file = open(filename, 'w')
    for result in range(2):
        for q1, q2, q3, q4 in itertools.product(range(4), repeat=4):
            file.write("%.14f %d %d %d %d %d\n"%(result_table[result][q1][q2][q3][q4], result, q1, q2, q3, q4))
    file.close()



# if __name__ == "__main__":
#     # error_prob = np.arange(0.005,0.015,0.001)
#     # np.random.seed(random.SystemRandom().randint(0,2**32-1))
#     # n_trials = 10000000
#     arrayID = float(sys.argv[1])
#     error_prob = (0.001*arrayID)+0.001
#     # get_error_table(1, 0.1)
#     # get_error_table(2, 0.3)
#     get_Fowler_error_table_exact(error_prob, 'X')
    # get_Fowler_error_table_exact(error_prob, 'Z')








# print(result_table)

# // //0: ancilla
# // //1-4: data qubits
# //std::array<int, 5> qubits = {0};
# //
# //std::array<double, 4> pauli_error_rate;
# //pauli_error_rate[0] = 1 - error_prob;
# //for (int i = 1; i < 4; ++i) pauli_error_rate[i] = error_prob / 3;
# //std::array<double, 16> cnot_error_rate;
# //cnot_error_rate[0] = 1 - error_prob;
# //for (int i = 1; i < 16; ++i) cnot_error_rate[i] = error_prob / 15;
# //
# //std::random_device rd;
# //std::mt19937 gen(rd());
# //std::discrete_distribution<> pauli_error_distr(pauli_error_rate.begin(), pauli_error_rate.end());
# //std::discrete_distribution<> cnot_error_distr(cnot_error_rate.begin(), cnot_error_rate.end());
# //std::discrete_distribution<> readout_error_distr({1-error_prob, error_prob});
# //
# // //Prep error on ancilla:
# //qubits[0] = pauli_error_distr(gen);
# //
# // //Hadamard error on data qubit:
# //for (int i = 1; i < 5; ++i) {
# //qubits[i] = pauli_error_distr(gen);
# //}
#
#
#
# //cnot error on qubit and ancilla:
# //}



