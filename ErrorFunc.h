//
// Created by Fish on 16/11/2017.
//

#ifndef TORICCODE_ERRORFUNC_H
#define TORICCODE_ERRORFUNC_H
enum DataError{
    NO_ERROR,
    X_ERROR,
    Y_ERROR,
    Z_ERROR
};

inline bool isError(int code_value, int error){
    return code_value == error or code_value == Y_ERROR;
}

inline int errorComposite(int error0, int error1){
    DataError error_f;
    if (error0 == NO_ERROR) error_f = (DataError)error1;
    else if (error0 == error1) error_f = NO_ERROR;
    else if (error1 == X_ERROR){
        if (error0 == Y_ERROR) error_f = Z_ERROR;
        else error_f = Y_ERROR;
    }
    else if (error1 == Z_ERROR){
        if (error0 == Y_ERROR) error_f = X_ERROR;
        else error_f = Y_ERROR;
    }
    else if (error1 == Y_ERROR){
        if (error0 == X_ERROR) error_f = Z_ERROR;
        else error_f = X_ERROR;
    }
    else error_f = (DataError)error0;
    return (int)error_f;
}
#endif //TORICCODE_ERRORFUNC_H
