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

int getSign(int a){
    return (a>0) - (a<0);
}


class Code{
public:
    std::vector< std::vector<int>> _code;
    int n_row;
    int n_col;
public:
    Code(int n_row, int n_col): n_row(n_row),n_col(n_col){
        _code.resize(n_row);
        for(int i = 0; i < n_row; i++){
            _code[i].resize(n_col);
            for(int j = 0; j < n_col; j++){
                _code[i][j] = (int)0;
            }
        }
    }
//REMEMBER to add destructor.
    int& code(int row, int col){
        row = (row% n_row + n_row)%n_row;
        col = (col% n_col + n_col)%n_col;
        return _code[row][col];
    }

    int& operator()(int row, int col){
        row = (row% n_row + n_row)%n_row;
        col = (col% n_col + n_col)%n_col;
        return _code[row][col];
    }

    std::array<int, 2> idxTransform(int row, int col){
        row = (row% n_row + n_row)%n_row;
        col = (col% n_col + n_col)%n_col;
        return {row, col};
    }

    void printCode(){
        printMatrix(_code);
    }
    void reset(){
        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                _code[i][j] = (int)0;
            }
        }
    }

    int minLength(const int start, const int end, const int n_L){
        int sign = (end > start) - (start > end);
        int d_1 = end - start;
        int d_2 = sign*(abs(d_1)-n_L);
        if (abs(d_1) < abs(d_2)) return d_1;
        else return d_2;
    }

    std::array<int, 2> minPath(std::array<int,2> start, std::array<int,2> end){
        return {minLength(start[0], end[0], n_row),minLength(start[1], end[1], n_col)};
    }

    int distance(std::array<int,2> loc1, std::array<int,2> loc2){
        int d;
        d = std::min(abs(loc1[0] - loc2[0]),
                     n_row - abs(loc1[0] - loc2[0])) +
            std::min(abs(loc1[1] - loc2[1]),
                     n_col - abs(loc1[1] - loc2[1]));
        return d;
    }

};

enum DataError{
    NO_ERROR = 0,
    X_ERROR = 1,
    Z_ERROR = -1,
    Y_ERROR = 2
};

bool isError(int code_value, int error){
    return code_value == error or code_value == 2;
}

int errorComposite(int error0, int error1){
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

enum StabiliserType{
    X_STB = 1,
    Z_STB= -1
};

class Data: public Code{
public:
    StabiliserType stabiliser_type;
    std::array<std::vector<std::array<int,2>>,2> error_locations;
public:
    Data(int n_row, int n_col): Code(n_row, n_col){
        assert(n_row%2 == 0);//for TORIC code, the number of rows and columns of data qubit must be even
        assert(n_col%2 == 0);
    }

    void induceError(double error_prob, int ERROR){
        assert(error_prob <= 1);
        int OTHER_ERROR = (int) -ERROR;

        std::random_device rd;
        std::mt19937 gen(rd());
//        std::mt19937 gen(4);
        std::binomial_distribution<> error_occur(1, error_prob);

        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                if (error_occur(gen)){
                    code(i, j) = errorComposite(code(i, j), ERROR);
                }
            }
        }
    }
    void getErrorLoc(){
        for (int i = 0; i < n_row; i++) {
            for (int j = 0; j < n_col; j++) {
                if (code(i, j) == -1 or code(i, j) == 2) {
                    error_locations[0].push_back({i,j});
                }
                else if (code(i, j) == 1 or code(i, j) == 2){
                    error_locations[1].push_back({i,j});
                }
            }
        }
    }
    void printErrorLoc(){
        getErrorLoc();
        printf("The location of Z errors are \n");
        for (auto loc: error_locations[0]){
            printf("(%d, %d)", loc[0], loc[1]);
        }
        printf("\n");
        printf("The location of X errors are \n");
        for (auto loc: error_locations[1]){
            printf("(%d, %d)", loc[0], loc[1]);
        }
    };

};



//boolean 0 denote no error, 1 denotes error.
class Stabiliser: public Code{
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
                    code(i, j) = not code(i, j);
                }
            }
        }
    }

    void getErrorLoc(){
        for (int i = 0; i < n_row; i++) {
            for (int j = 0; j < n_col; j++) {
                if (code(i, j) != 0) {
                    error_locations.push_back({i,j});}
            }
        }
    }

    PerfectMatching* getErrorMatching(){
        getErrorLoc();
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
                int d = distance(error_locations[i], error_locations[j]);
                pm->AddEdge(error_label[i], error_label[j], d);
            }
        }
        pm->Solve();

        return pm;
    }
};

class SurfaceCode{
public:
    Data data;
    Stabiliser stabiliserX;
    Stabiliser stabiliserZ;
    int n_row;
    int n_col;
private:
    std::array<int,2> error_chain_target;

public:
    SurfaceCode(int n_row, int n_col): n_row(n_row), n_col(2*n_col), data(n_row, n_col),
                                       stabiliserX(n_row/2, n_col, X_STB), stabiliserZ(n_row/2, n_col, Z_STB){}

    int& code(int row, int col){
        row = (row% n_row + n_row)%n_row;
        col = (col% n_col + n_col)%n_col;
        if (row%2 == 0 and col%2 == 0) return stabiliserX.code(row/2,col/2);
        else if (row%2 == 1 and col%2 == 1) return stabiliserZ.code((row-1)/2, (col-1)/2);
        else if (row%2 == 0) return data.code(row,(col-1)/2);
        else if (row%2 == 1) return data.code(row, col/2);
    }

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
                    if (data(i, j)!= Z_ERROR){
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
        std::vector<std::array<int, 2>> pos_array= {{0,1}, {0,(-1)}, {1,0}, {-1,0}};
        int nX, nZ;
        for (int i = 0; i < stabiliserX.n_row; ++i) {
            for (int j = 0; j < stabiliserX.n_col; ++j) {
                nX = 0;
                nZ = 0;
                for (std::array<int,2> pos: pos_array){
                    if (code(2*i+pos[0],2*j+pos[1]) == Z_ERROR or code(2*i+pos[0],2*j+pos[1]) == Y_ERROR) nX++;
                    if (code(2*i+1+pos[0],2*j+1+pos[1]) == X_ERROR or code(2*i+1+pos[0],2*j+1+pos[1]) == Y_ERROR) nZ++;
                }
                if (nX%2 == 1) stabiliserX(i,j) = 1;
                else stabiliserX(i,j) = 0;
                if (nZ%2 == 1) stabiliserZ(i,j) = 1;
                else stabiliserZ(i,j) = 0;
            }
        }
    }

    //when annihilating errors, we aways go in row direction first, then in col direction.
    void fixError(StabiliserType stabiliser_type){
        Stabiliser* stabiliser;
        DataError  ERROR;
        int offset;
        if (stabiliser_type == X_STB){
            stabiliser = &stabiliserX;
            ERROR = Z_ERROR;
            offset = 0;
        }
        else{
            stabiliser = &stabiliserZ;
            ERROR = X_ERROR;
            offset = 1;
        }

        PerfectMatching* pm = stabiliser->getErrorMatching();
        int n_error = stabiliser->error_locations.size();
        std::vector <int> error_label;
        std::array <int, 2> start, end, loc, min_path;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::binomial_distribution<> direction(1, 0.5);


        int  ver_sign, hor_sign, ver_steps, hor_steps;
        for (int k = 0; k < n_error; ++k) error_label.push_back(k);
        std::set<int> error_corrected; //std::array are implicitly copied.
        for (int error: error_label){
            if (error_corrected.find(error) == error_corrected.end()) { // equivalent to if error is not in error_corrected
                int paired_error = pm->GetMatch(error);
                min_path = stabiliser->minPath(stabiliser->error_locations[error],
                                               stabiliser->error_locations[paired_error]);
                ver_sign = getSign(min_path[0]); //+1 if end[0]>start[0], -1 if otherwise.
                hor_sign = getSign(min_path[1]); //+1 if end[0]>start[0], -1 if otherwise.

                start = {2* stabiliser->error_locations[error][0] + offset, 2* stabiliser->error_locations[error][1] + offset};
                end = {2* stabiliser->error_locations[paired_error][0] + offset, 2* stabiliser->error_locations[paired_error][1] + offset};
                loc = start; //implicit copy of start

                //We might want to use step_sign *loc[0] < step_sign * end[0] condition to substitue ver_steps. But the
                //situation is much more complicated when we need to cross the boundary.
                ver_steps = 0;
                hor_steps = 0;
// ///////////////////// Comment out this section if we don't want error correction along random min path.
                while (ver_steps < abs(min_path[0]) and hor_steps < abs(min_path[1])){
                    if (direction(gen)){
                        code(loc[0] + ver_sign, loc[1]) = errorComposite(code(loc[0] + ver_sign, loc[1]), ERROR);
                        loc[0] += 2*ver_sign;
                        ver_steps +=1;
                    }
                    else{
                        code(loc[0], loc[1] + hor_sign) = errorComposite(code(loc[0], loc[1] + hor_sign), ERROR);
                        loc[1] += 2*hor_sign;
                        hor_steps +=1;
                    }
                }
// ///////////////////////////////
                while (ver_steps < abs(min_path[0])){
                    code(loc[0] + ver_sign, loc[1]) = errorComposite(code(loc[0] + ver_sign, loc[1]), ERROR);
                    loc[0] += 2*ver_sign;
                    ver_steps +=1;
                }
                while (hor_steps < abs(min_path[1])){
                    code(loc[0], loc[1] + hor_sign) = errorComposite(code(loc[0], loc[1] + hor_sign), ERROR);
                    loc[1] += 2*hor_sign;
                    hor_steps +=1;
                }
//                error_corrected.insert(error); //error is already iterated over, hence no need to check.
                error_corrected.insert(paired_error);
            }
        }
    }
    /*
     * We mark the endpoint using 30
     * when we visit a point, we add 10 to it.
     * 1 ->11, -1 -> 9, 2 -> 12
     * when starting point and end point are different, we leave the starting point as it is, and a 10 will be added onto
     * it since we've visited it just like any other point on the path.
     *
     * when the starting point and the end point at the same. We set the starting point to be 20. After we first visited it
     * 10 is added to the starting point and it becomes 30, which is the end point.
     *
     * we use error_chain_target to denote the end point instead of using some special number. This is because
     * when we want to use the starting point as endpoint, then the error chain will end at the starting point
     */


//    bool findPath(int row, int col, int error){ // row and col are coor in the surface code grid.
//       /** Accept case - we found the exit **/
//        if (code(row, col) == 30) return true;
// //        else if (code(row, col) == 19) code(row, col) +=1;
//
//        /** Reject case - we hit a wall or our path **/
//        if (code(row, col) != error and code(row, col) != 2 and code(row, col) != 20) return false;
//
//        /** Backtracking Step **/
//
//        // Mark this location as part of our path
//        code(row, col) += 10;
//
//        std::cout<<error;
//
//        std::vector<std::array<int,2>> nb_loc;
//        if (error == -1){
//            if (row%2 == 0) nb_loc = {{row, col +2}, {row -1, col +1},{row+1, col+1}};
//            else nb_loc = {{row+2, col}, {row +1, col -1},{row+1, col+1}};
//        }
//        else if (error == 1){
//            if (row%2 == 1) nb_loc = {{row, col +2}, {row -1, col +1},{row+1, col+1}};
//            else nb_loc = {{row+2, col}, {row +1, col -1},{row+1, col+1}};
//        }
//        //Continue with our path
//        for (std::array<int,2> loc: nb_loc){
//            if (findPath(loc[0], loc[1], error)) return true;
//        }
//
//        /** Deadend - this location can't be part of the solution **/
//
//        // Unmark this location
// //        code(row, col) -= 10;
//
//        // Go back
//        return false;
//    }
//
//
//
//

//    int hasLogicalError(StabiliserType stb){
//        std::vector< std::vector<int>> data_copy = data._code; //implicit copy of vector.
//        int error;
//        int n_logical_error = 0;
//
//        if (stb == X_STB) error = -1;
//        else if (stb == Z_STB) error = 1;
//
//        for (int i = 0; i < n_row; i++) {
//            for (int j = 0; j < n_col; j++) {
//                if (isError(code(i,j), error)) {
//                    code(i,j) = 20;
//                    if(findPath(i,j,error)) n_logical_error+=1;
//                }
//            }
//        }
//        data._code = data_copy;
//        //Identify the errors at the edges first. For logical error to exist, there must be at least one error each
//        // on left and right with similar row index or one error each at top and bottom with similar column index.
//        return n_logical_error;
//    }


// Simply counting the number of error in each row or col is not a good way because this fail when there is a portion of
// path running along the line that we are testing.
// way to counter this: testing along more rows? (NOT conclusive)
// A better way is to track along the path that run along the row?
// The portion of path can be identified by looking at consecutive error along one row.

/* For X_stb, we count along the row 1, for which Z_errors cannot run along the row (unlike row 0
 * for which Z_errors can run along the row going through X_stb.
 * Similarly for Z_stb
 * */
    bool hasLogicalError(StabiliserType stb){
        std::vector< std::vector<int>> data_copy = data._code; //implicit copy of vector.
        int n_v_error = 0;
        int n_h_error = 0;
        if (stb == X_STB){
            for (int j = 0; j < n_col; j +=2) {
                if (isError(code(1,j), -1)) n_v_error++; // check horizontally if any vertical error cuts through
            }
            for (int i = 0; i < n_row; i += 2) {
                if (isError(code(i,1), -1)) n_h_error++;
            }
        }
        else if (stb == Z_STB){
            for (int j = 1; j < n_col; j +=2) {
                if (isError(code(0,j), 1)) n_v_error++;
            }
            for (int i = 1; i < n_row; i += 2) {
                if (isError(code(i,0), 1)) n_h_error++;
            }
        }
        printf("there are %d horizontal logical errors\n", n_h_error%2);
        printf("there are %d vertical logical errors\n", n_v_error%2);
        return ((n_h_error%2) or (n_v_error%2));
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
    c.data.induceError(0.1, X_ERROR);
    c.data.printCode();
    c.data.induceError(0.1, Z_ERROR);
    c.data.printCode();
    c.stabiliserUpdateSlow();
    c.printSurfaceCode();
    c.fixError(X_STB);
    c.stabiliserUpdateSlow();
    c.printSurfaceCode();
    c.fixError(Z_STB);
    c.stabiliserUpdateSlow();
    c.printSurfaceCode();
//    c.data.printCode();
    int X_log_error = c.hasLogicalError(X_STB);
    printf("there are %d X logical errors\n", X_log_error);
    printf("\n");
    int Z_log_error = c.hasLogicalError(Z_STB);
    printf("there are %d Z logical errors\n", Z_log_error);
//    c.printSurfaceCode();
//    int Z_log_error = c.hasLogicalError(Z_STB);
//    printf("there are %d Z logical errors\n", Z_log_error);

//    c.data.printErrorLoc();
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