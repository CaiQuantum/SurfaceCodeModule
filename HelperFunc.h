//
// Created by Fish on 16/11/2017.
//

#ifndef TORICCODE_HELPERFUNC_H
#define TORICCODE_HELPERFUNC_H
#include <iostream>
template <class InputType> // printMatrix can be printed without specifying InputType, this is done using implicit instantiation
void printMatrix(InputType M){
    for (auto array: M){
        for (auto element: array){
            printf("%2d ", element);
        }
        printf("\n");
    }
    printf("\n");
}

inline int getSign(int a){
    return (a>0) - (a<0);
}






#endif //TORICCODE_HELPERFUNC_H
