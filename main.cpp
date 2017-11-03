#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <array>
#include <set>
#include "PerfectMatching.h"
//To make perfectmatching.h to work, we need to first delete example.cpp in blossom_dir, then we also cannot use the
//triangle package as suggested due to lack of X11. We use import project in Clion to rewrite Cmake. Remember to exclude
//the unwanted files.

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


template <class ErrorType>
class Code{
public:
    int n_row;
    int n_col;
    std::vector< std::vector<ErrorType>> code;
public:
    Code(int n_row, int n_col): n_row(n_row),n_col(n_col){
        code.resize(n_row);
        for(int i = 0; i < n_row; i++){
            code[i].resize(n_col);
            for(int j = 0; j < n_col; j++){
                code[i][j] = (ErrorType)0;
            }
        }
    }
    ~Code() {}

    ErrorType& operator()(int row, int col){
        row = (row% n_row + n_row)%n_row;
        col = (col% n_col + n_col)%n_col;
        return code[row][col];
    }

    std::array<int, 2> idxTransform(int row, int col){
        row = (row% n_row + n_row)%n_row;
        col = (col% n_col + n_col)%n_col;
        return {row, col};
    }

    void printCode(){
        printMatrix(code);
    }
    void reset(){
        for(int i = 0; i < n_row; i++){
            code[i].resize(n_col);
            for(int j = 0; j < n_col; j++){
                code[i][j] = (ErrorType)0;
            }
        }
    }
};

enum DataError{
    NO_ERROR = 0,
    X_ERROR = 1,
    Z_ERROR = -1,
    Y_ERROR = 2
};

enum StabiliserType{
    X_STB = 1,
    Z_STB= -1
};

class Data: public Code<DataError>{
public:
    Data(int n_row, int n_col): Code(n_row, n_col){
        assert(n_row%2 == 0);//for TORIC code, the number of rows and columns of data qubit must be even
        assert(n_col%2 == 0);
    }

    void induceError(double error_prob, DataError ERROR){
        assert(error_prob <= 1);
        DataError OTHER_ERROR = (DataError) -ERROR;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::binomial_distribution<> error_occur(1, error_prob);

        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                if (error_occur(gen)){
                    if (code[i][j] == NO_ERROR) code[i][j] = ERROR;
                    else if (code[i][j] == ERROR) code[i][j] = NO_ERROR;
                    else if (code[i][j] == OTHER_ERROR) code[i][j] = Y_ERROR;
                    else code[i][j] = OTHER_ERROR;
                }
            }
        }
    }
    std::array<int,3> neighbour(const char direction, int row, int col) {
        std::array<int, 3> nb;
        if (row%2 == 0) {
            if (direction == 'N') nb = {row/2 - 1, col, (int)Z_STB};
            else if (direction == 'S') nb = {row/2, col, (int)Z_STB};
            else if (direction == 'W') nb = {row/2, col, (int)X_STB};
            else if (direction == 'E') nb = {row/2, col+1, (int)X_STB};
        }
        else{
            if (direction == 'N') nb = {(row-1)/2, col, (int)X_STB};
            else if (direction == 'S') nb = {(row+1)/2, col, (int)X_STB};
            else if (direction == 'W') nb = {row, col-1, (int)Z_STB};
            else if (direction == 'E') nb = {row, col, (int)Z_STB};
        }
        return nb;
    }
//    std::array<int,2> neighbour(const StabiliserType stabiliser_type, int row, int col) {
//        std::array<int, 2> nb;
//
//        if (row%2 == 0) {
//            if (direction == 'N') nb = {row/2 - 1, col, (int)Z_STB};
//            else if (direction == 'S') nb = {row/2, col, (int)Z_STB};
//            else if (direction == 'W') nb = {row/2, col, (int)X_STB};
//            else if (direction == 'E') nb = {row/2, col+1, (int)X_STB};
//        }
//        else{
//            if (direction == 'N') nb = {(row-1)/2, col, (int)X_STB};
//            else if (direction == 'S') nb = {(row+1)/2, col, (int)X_STB};
//            else if (direction == 'W') nb = {row, col-1, (int)Z_STB};
//            else if (direction == 'E') nb = {row, col, (int)Z_STB};
//        }
//        return nb;
//    }

};



//boolean 0 denote no error, 1 denotes error.
class Stabiliser: public Code<int>{
public:
    StabiliserType stabiliser_type;
    std::vector<std::array<int,2>> error_locations;
    
public:
//    Stabiliser(int n_row, int n_col): Code(n_row, n_col){}
    Stabiliser(int n_row, int n_col, StabiliserType stabiliser_type): Code(n_row, n_col), stabiliser_type(stabiliser_type){}
    void induceError(double error_prob){
        assert(error_prob <= 1);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::binomial_distribution<> error_occur(1, error_prob);

        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                if (error_occur(gen)){
                    code[i][j] = not code[i][j];
                }
            }
        }
    }

    PerfectMatching* getErrorMatching(){

        for (int i = 0; i < n_row; i++) {
            for (int j = 0; j < n_col; j++) {
                if (code[i][j] != 0) {
                    error_locations.push_back({i,j});}
            }
        }

        int n_error = error_locations.size();
        std::vector <int> error_label;
        for (int k = 0; k < n_error; ++k) error_label.push_back(k);

        int n_edges = n_error*(n_error-1)/2;
        PerfectMatching *pm = new PerfectMatching(n_error, n_edges);

        struct PerfectMatching::Options options;
        options.verbose = false;
        pm->options = options;

        for (int i = 0; i < n_error; ++i) {
            for (int j = 0; j < i ; ++j) {
                int d;
                d = std::min(abs(error_locations[i][0] - error_locations[j][0]),
                             n_row - abs(error_locations[i][0] - error_locations[j][0])) +
                    std::min(abs(error_locations[i][1] - error_locations[j][1]),
                             n_col - abs(error_locations[i][1] - error_locations[j][1]));
                pm->AddEdge(error_label[i], error_label[j], d);
            }
        }
        pm->Solve();

        return pm;
    }
    std::array<int,2> neighbour(const char direction, int row, int col) {
        std::array<int, 2> nb_pos;
        if (stabiliser_type == X_STB) {
            if (direction == 'N') nb_pos = {2 * row - 1, col};
            else if (direction == 'S') nb_pos = {2 * row + 1, col};
            else if (direction == 'W') nb_pos = {2 * row, col - 1};
            else if (direction == 'E') nb_pos = {2 * row, col};
        }
        if (stabiliser_type == Z_STB) {
            if (direction == 'N') nb_pos = {2 * row, col};
            else if (direction == 'S') nb_pos = {2 * row + 2, col};
            else if (direction == 'W') nb_pos = {2 * row + 1, col};
            else if (direction == 'E') nb_pos = {2 * row + 1, col + 1};
        }
        return nb_pos;
    }

};

class SurfaceCode{
public:
    Data data;
//    std::array<Stabiliser, 2> stabiliser_array;
//    Stabiliser& stabiliserX = stabiliser_array[0];
//    Stabiliser& stabiliserZ = stabiliser_array[1];
    Stabiliser stabiliserX;
    Stabiliser stabiliserZ;

public:
    SurfaceCode(int n_row, int n_col): data(n_row, n_col), stabiliserX(n_row/2, n_col, X_STB),
                                                           stabiliserZ(n_row/2, n_col, Z_STB){}


    //here assume row 0 is X stabiliser, row 1 is Z stabliser, etc.
    void stabiliserUpdate(){
        for (int i = 0; i < data.n_row; i++) {
            for (int j = 0; j < data.n_col; j++) {
                if (data(i, j) != NO_ERROR) {
                    if (data(i, j) != X_ERROR){ //include both the case of Z_ERROR and Y_ERROR
                        if (i%2 == 0) {
                            stabiliserX(i/2,j) ^= 1;
                            stabiliserX(i/2,j+1) ^= 1;
                        }
                        else {
                            stabiliserX((i+1)/2,j) ^= 1;
                            stabiliserX((i-1)/2,j) ^= 1;
                        }
                    }
                    if (data.code[i][j] != Z_ERROR){
                        if (i%2 == 0) {
                            stabiliserZ(i/2,j) ^= 1;
                            stabiliserZ((i/2-1),j) ^= 1;
                        }
                        else {
                            stabiliserZ((i-1)/2, j) ^= 1;
                            stabiliserZ((i-1)/2, j-1) ^= 1;
                        }
                    }
                }
            }
        }
    }

    void stabiliserUpdateSlow(){
        std::array<char, 4> direction_array= {'N', 'S', 'E', 'W'};
        int nX, nZ;
        std::array<int, 2> pos;
        for (int i = 0; i < stabiliserX.n_row; ++i) {
            for (int j = 0; j < stabiliserX.n_col; ++j) {
                nX = 0;
                nZ = 0;
                for (char direction: direction_array){
                    pos = stabiliserX.neighbour(direction, i, j);
                    if (data(pos[0],pos[1]) == Z_ERROR or data(pos[0],pos[1]) == Y_ERROR) nX++;
                    pos = stabiliserZ.neighbour(direction, i, j);
                    if (data(pos[0],pos[1]) == X_ERROR or data(pos[0],pos[1]) == Y_ERROR) nZ++;
                }
                if (nX%2 == 1) stabiliserX(i,j) = 1;
                else stabiliserX(i,j) = 0;
                if (nZ%2 == 1) stabiliserZ(i,j) = 1;
                else stabiliserZ(i,j) = 0;
            }
        }
    }

    //when annihilating errors, we aways go in row direction first, then in col direction.
    void fixError(){
        int n_error = stabiliserX.error_locations.size();
        PerfectMatching* pm = stabiliserX.getErrorMatching();
        std::vector <int> error_label;
        for (int k = 0; k < n_error; ++k) error_label.push_back(k);
        std::set<int>error_corrected = {}; //std::array are implicitly copied.
        for (int error: error_label){
            if (error_corrected.find(error) == error_corrected.end()) { // equivalent to if error is not in error_corrected
                int paired_error = pm->GetMatch(error);
//                std::array<int,2>


//                code[error_locations[error][0]][error_locations[error][1]] ^= 1;
//                code[error_locations[paired_error][0]][error_locations[paired_error][1]] ^= 1;
//                error_corrected.insert(error); //error is already iterated over, hence no need to check.
                error_corrected.insert(paired_error);
            }
        }
    }

//red: 31, grn: 32, yel: 33, blu: 34, mag: 35, cyn: 36, wht: 37
    void printSurfaceCode(){
        for (int i = 0; i < data.n_row; i++) {
            for (int j = 0; j < data.n_col; j++) {
                if (i%2 == 0) {
                    printf("\x1B[31m%2d\x1B[0m ", stabiliserX(i/2, j));
                    printf("%2d ", data(i, j));
                }
                else {
                    printf("%2d ", data(i, j));
                    printf("\x1B[34m%2d\x1B[0m ", stabiliserZ((i-1)/2, j));
                }
            }
            printf("\n");
        }
        printf("\n");
    }


};


int main() {
    SurfaceCode c(8, 4);
    c.data.printCode();
    c.stabiliserX.printCode();
    c.stabiliserZ.printCode();
    c.data.induceError(0.2, X_ERROR);
    c.data.printCode();
    c.data.induceError(0.2, Z_ERROR);
    c.data.printCode();
    c.stabiliserUpdate();
    c.printSurfaceCode();
//    std::array<int,2> pos = c.stabiliserX.neighbour('N', 0, 0);
//    printf("(%d, %d)", pos[0], pos[1]);
    c.stabiliserX.reset();
    c.stabiliserZ.reset();
    c.stabiliserUpdateSlow();
    c.printSurfaceCode();
//    c.stabiliserX.printError();
//    c.stabiliserX.printCode();
//    c.stabiliserZ.printCode();
//    c.printSurfaceCode();


//    c.stabiliserX.fixError();
//    c.stabiliserZ.fixError();
//    c.stabiliserX.printCode();
//    c.stabiliserZ.printCode();
//    c.printSurfaceCode();
}


/*
 * What we want:
 *
 * A CodeArray object:
 * CodeArray.init(height, width)
 * CodeArray.induceErrors(error_percentage, error_type)
 * CodeArray.error_location: an array recording the location of errors
 *
 * A DataArray object: inherent from code array
 *
 * A StabiliserArray object:
 * inherent form code array, with members like StabiliserArray.stabiliser_type, StabiliserArray.relative_location
 * stabiliser.error_location
 * stabiliser.getErrorLocation(): return the location of stabiliser errors.
 * stabiliser.errorPairing(): using min-weight to pair up errors.
 *
 * A SurfaceCode object:
 *
 * SurfaceCode.data: an DataArray object of data qubits
 * SurfaceCode.stabiliser(n): an StabiliserArray object of stabiliser qubits, n is 0 or 1.
 * SurfaceCode.stabiliserUpdate(): which is just stabiliser measurement
 * SurfaceCode.errorFix(): Fix error in SurfaceCode.data using stabiliser.errorPairing()
 * SurfaceCode.
 *
 *
 */