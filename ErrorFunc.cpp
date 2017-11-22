//
// Created by Fish on 21/11/2017.
//
#include <fstream>
#include <iostream>
#include <array>
#include "ErrorFunc.h"
#include <tuple>

bool isError(int code_value, int error){
    return code_value == error or code_value == Y_ERROR;
}

int prod_table[4][4] = {{0,1,2,3},
                        {1,0,3,2},
                        {2,3,0,1},
                        {3,2,1,0}};

int errorComposite(const int error0, const int error1) {
    return prod_table[error0][error1];
}

int cnot_err_table[4][4][2] = {{{0,0}, {0,1}, {3,2}, {3,3}},
                               {{1,1}, {1,0}, {2,3}, {2,2}},
                               {{2,1}, {2,0}, {1,3}, {1,2}},
                               {{3,0}, {3,1}, {0,2}, {0,3}}};
std::tuple<int, int> pass_through_cnot(const int cbit, const int nbit){
    return  std::make_tuple(cnot_err_table[cbit][nbit][0], cnot_err_table[cbit][nbit][1]);
}


int pass_through_H(const int error){
    if (error == 1) return 3;
    else if (error == 3) return 1;
    else return error;
}


//inline int errorComposite(int error0, int error1) {
//    DataError error_f;
//    if (error0 == NO_ERROR) error_f = (DataError)error1;
//    else if (error0 == error1) error_f = NO_ERROR;
//    else if (error0 == X_ERROR and error1 == Y_ERROR) error_f = Z_ERROR;
//    else if (error0 == X_ERROR and error1 == Z_ERROR) error_f = Y_ERROR;
//    else if (error0 == Y_ERROR and error1 == Z_ERROR) error_f = X_ERROR;
//    else return errorComposite(error1, error0);
//    return (int)error_f;
//}

//std::array<std::array<std::array<int, 2>, 4>,4> cnot_err_table;
//
//void load_cerr_table() {
//    std::ifstream inFile;
//    //If file can't be read, maybe the filename buffer is NOT LONG ENOUGH!!!
//    char filename[100];
//    sprintf(filename, "../ParityCheckErrorTable/cnot_err_table.txt");
//    inFile.open(filename);
//    if (!inFile) {
//        std::cerr << "unable to open file for reading" << std::endl;
//    }
//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 4; j++) {
//            for (int k = 0; k < 2; k++) {
//                inFile >> cnot_err_table[i][j][k];
//            }
//        }
//    }
//    inFile.close();
//}
