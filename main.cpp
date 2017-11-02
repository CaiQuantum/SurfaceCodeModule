#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <array>

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

    void printCode(){
        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                std::cout << code[i][j]<< " ";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }
};

enum DataError{
    NO_ERROR = 0,
    X_ERROR = 1,
    Z_ERROR = -1,
    Y_ERROR = 2
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

};

//enum StabiliserType{
//    X_MEASUREMENT = 1,
//    Z_MEASUREMENT= -1
//};

//boolean 0 denote no error, 1 denotes error.
class Stabiliser: public Code<int>{
public:
    std::vector<std::array<int,2>> error_locations;
//    StabiliserType stabiliser_type;
    
public:
    Stabiliser(int n_row, int n_col): Code(n_row, n_col){}
//    Stabiliser(int n_row, int n_col, StabiliserType stabiliser_type): Code(n_row, n_col), stabiliser_type(stabiliser_type){}
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
    void getError(){
        for (int i = 0; i < n_row; i++) {
            for (int j = 0; j < n_col; j++) {
                if (code[i][j] != 0) {
                    std::array<int, 2> loc = {i,j};
                    error_locations.push_back(loc);}
            }
        }
    }

//    void printError(){
//        this->getError();
//        for (int i = 0; i < error_locations.size(); ++i) {
//            printf("(%d, %d) ", error_locations[i][0], error_locations[i][1]);
//        }
//    }
    void printError(){
        this->getError();
        for (auto error: error_locations) {
            printf("(%d, %d) ", error[0], error[1]);
//            std::cout<<"("<<error[0]<<", "<<error[1]<<")";
        }
        std::cout<<std::endl;
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
    SurfaceCode(int n_row, int n_col): data(n_row, n_col), stabiliserX(n_row/2, n_col), stabiliserZ(n_row/2, n_col){}


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
    c.stabiliserX.getError();
    c.stabiliserX.printError();
    c.stabiliserX.printCode();
    c.stabiliserZ.printCode();
    c.printSurfaceCode();
//    std::cout<<c(1,1)<<std::endl;
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