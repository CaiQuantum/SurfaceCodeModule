//
// Created by Fish on 16/11/2017.
//

#ifndef TORICCODE_ERRORFUNC_H
#define TORICCODE_ERRORFUNC_H
#include <fstream>
enum DataError{
    NO_ERROR,
    X_ERROR,
    Y_ERROR,
    Z_ERROR
};

bool isError(int code_value, int error);

int errorComposite(int error0, int error1);

void load_cerr_table();

std::tuple<int, int> pass_through_cnot(int cbit, int nbit);

int pass_through_H(int error);

#endif //TORICCODE_ERRORFUNC_H
