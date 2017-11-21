#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <set>
#include "PerfectMatching.h"
#include <fstream>
#include "ErrorFunc.h"
#include "HelperFunc.h"

/*To make perfectmatching.h to work, we need to first delete example.cpp in blossom_dir, then we also cannot use the
 * triangle package as suggested due to lack of X11. We use import project in Clion to rewrite Cmake. Remember to exclude
 * the unwanted files.
 *
 * There are two parameters you can tweak in the final version of this program. The first one is the time-spatial relative
 * weight in perfectmatching. i.e. how important is time distance relative to spatial distance when we are calculating
 * the perfect matching distance.
 * */

class Code{
public:
    std::vector<std::vector<int>> _code;
    int n_row;
    int n_col;
public:
    Code(int n_row, int n_col): n_row(n_row),n_col(n_col){
        _code.resize(n_row);
        for(int i = 0; i < n_row; i++){
            _code[i].resize(n_col);
            for(int j = 0; j < n_col; j++){
                _code[i][j] = NO_ERROR;
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
    void resetCode(){
        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                _code[i][j] = NO_ERROR;
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

    std::array<int, 2> minSpatialPath(std::array<int,3> start, std::array<int,3> end){
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
    int spatialDistance(std::array<int,3> loc1, std::array<int,3> loc2){
        int d;
        d = std::min(abs(loc1[0] - loc2[0]),
                     n_row - abs(loc1[0] - loc2[0])) +
            std::min(abs(loc1[1] - loc2[1]),
                     n_col - abs(loc1[1] - loc2[1]));
        return d;
    }

};

enum StabiliserType{
    X_STB,
    Z_STB
};

class Data: public Code{
public:
    std::array<std::vector<std::array<int,2>>,2> error_locations;
public:
    Data(int n_row, int n_col): Code(n_row, n_col){
        assert(n_row%2 == 0);//for TORIC code, the number of rows and columns of data qubit must be even
        assert(n_col%2 == 0);
    }

    void induceError(double error_prob, int ERROR){
        assert(error_prob <= 1);
//        auto OTHER_ERROR = (int) -ERROR;

        std::random_device rd;
        std::mt19937 gen(rd());
//        std::mt19937 gen(4);
        std::binomial_distribution<> error_occur(1, error_prob);

        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                if (error_occur(gen)){
                    _code[i][j] = errorComposite(_code[i][j], ERROR);
                }
            }
        }
    }
    void induceError(double error_prob){
        induceError(error_prob, X_ERROR);
        induceError(error_prob, Z_ERROR);
    }
    void getErrorLoc(){
        for (int i = 0; i < n_row; i++) {
            for (int j = 0; j < n_col; j++) {
                if (isError(_code[i][j], Z_ERROR)) {
                    error_locations[0].push_back({i,j});
                }
                else if (isError(_code[i][j], X_ERROR)){
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
    std::vector<std::array<int,3>> flip_locs;
    std::vector< std::vector<int>> last_code;
    int t = 0;
public:
//    Stabiliser(int n_row, int n_col): Code(n_row, n_col){}
    Stabiliser(int n_row, int n_col, StabiliserType stabiliser_type): Code(n_row, n_col), stabiliser_type(stabiliser_type){
        last_code = _code;
    }
    void induceError(double error_prob){
        assert(error_prob <= 1);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::binomial_distribution<> error_occur(1, error_prob);

        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                if (error_occur(gen)){
                    _code[i][j] = not _code[i][j];
                }
            }
        }
    }

    void addFlipLoc(){
        for (int i = 0; i < n_row; i++) {
            for (int j = 0; j < n_col; j++) {
                if (_code[i][j] != last_code[i][j]) {
                    flip_locs.push_back({i,j,t});}
            }
        }
    }

    PerfectMatching* getErrorMatching(){
//        assert (flip_locs.size() != 0);
        int n_error = flip_locs.size();
        std::vector <int> error_label(n_error);
        for (int k = 0; k < n_error; ++k) error_label[k] = k;

        int n_edges = n_error*(n_error-1)/2;
        auto *pm = new PerfectMatching(n_error, n_edges);

        struct PerfectMatching::Options options;
        options.verbose = false;
        pm->options = options;

        for (int i = 0; i < n_error; ++i) {
            for (int j = 0; j < i ; ++j) {
                //add spatial distance and time distance together as cost.
                int d = spatialDistance(flip_locs[i], flip_locs[j]) + abs(flip_locs[i][2] - flip_locs[j][2]);
                pm->AddEdge(error_label[i], error_label[j], d);
            }
        }
        pm->Solve();

        return pm;
    }
};

class ToricCode{
public:
    Data data;
    Stabiliser stabiliserX;
    Stabiliser stabiliserZ;
    int n_row;
    int n_col;

public:
    ToricCode(int n_row, int n_col): n_row(n_row), n_col(n_col), data(n_row, n_col/2),
                                       stabiliserX(n_row/2, n_col/2, X_STB), stabiliserZ(n_row/2, n_col/2, Z_STB){
        assert(n_row%2==0);
        assert(n_col%2==0);
    }

    int& code(int row, int col){
        row = (row% n_row + n_row)%n_row;
        col = (col% n_col + n_col)%n_col;
        if (row%2 == 0 and col%2 == 0) return stabiliserX._code[row/2][col/2];
        else if (row%2 == 1 and col%2 == 1) return stabiliserZ._code[(row-1)/2][(col-1)/2];
        else if (row%2 == 0) return data._code[row][(col-1)/2];
        else return data._code[row][col/2]; //equivalent to else if (row%2 == 1)
    }
    //red: 31, grn: 32, yel: 33, blu: 34, mag: 35, cyn: 36, wht: 37
    void printCode(){
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

    void stabiliserUpdate(){
        std::vector<std::array<int, 2>> pos_array= {{0,1}, {-1,0}, {1,0}, {0,(-1)}};
        int nX, nZ;
        for (int i = 0; i < stabiliserX.n_row; ++i) {
            for (int j = 0; j < stabiliserX.n_col; ++j) {
                nX = 0;
                for (std::array<int,2> pos: pos_array){
                    if (isError(code(2*i+pos[0],2*j+pos[1]), Z_ERROR)) nX++;
                }
                stabiliserX(i,j) = nX%2;
            }
        }
        for (int i = 0; i < stabiliserZ.n_row; ++i) {
            for (int j = 0; j < stabiliserZ.n_col; ++j) {
                nZ = 0;
                for (std::array<int,2> pos: pos_array){
                    if (isError(code(2*i+1+pos[0],2*j+1+pos[1]), X_ERROR)) nZ++;
                }
                stabiliserZ(i,j) = nZ%2;
            }
        }
    }
    void stabiliserUpdate(double error_prob){

        std::array<double,512> prob_array;
        std::array<std::array<int, 5>,512> error_table;

        std::ifstream inFile;

        //If file can't be read, maybe the filename buffer is NOT LONG ENOUGH!!!
        char filename [100];
//        sprintf (filename, "../ParityCheckErrorTable/Data/exact_error_table%.4f.txt", error_prob);
        sprintf (filename, "../ParityCheckErrorTable/Data/new_exact_error_table%.4f.txt", error_prob);
        inFile.open(filename);
        if (! inFile) {
            std::cerr << "unable to open file for reading" << std::endl;
        }
        for (int i = 0; i < 512; i++) {
            inFile >> prob_array[i];
            for (int j = 0; j < 5; j++) {
                inFile >> error_table[i][j];
            }
        }
        inFile.close();

        std::random_device rd;
        std::mt19937 gen(rd());
        std::discrete_distribution<> error_distr(prob_array.begin(), prob_array.end());
        std::vector<std::array<int, 2>> pos_array = {{0,1}, {0,(-1)}, {1,0}, {-1,0}};
        std::array<int, 5> X_stb_err;
        int nX, nZ;
        int ancilla;
        for (int i = 0; i < stabiliserX.n_row; ++i) {
            for (int j = 0; j < stabiliserX.n_col; ++j) {
                ancilla= 0;
//                nX=0;
                //randomly select one row in the error_table
                X_stb_err = error_table[error_distr(gen)];
                int k = 1;
                for (std::array<int,2> pos: pos_array){
//                    //calculate real parity
//                    if (isError(code(2*i+pos[0],2*j+pos[1]), Z_ERROR)) nX++;

                    '''add in hadamard!!!!!'''

                    std::tie(code(2*i+pos[0],2*j+pos[1]), ancilla) = pass_through_cnot(code(2*i+pos[0],2*j+pos[1]),ancilla);

                    //add errors to data qubit
                    code(2*i+pos[0],2*j+pos[1]) = errorComposite(code(2*i+pos[0],2*j+pos[1]), X_stb_err[k]);
                    k++;
                }
                //add readout errors
                if (ancilla == )
                stabiliserX(i,j) = (nX+X_stb_err[0])%2;
            }
        }
        std::array<int, 5> Z_stb_err;
        for (int i = 0; i < stabiliserZ.n_row; ++i) {
            for (int j = 0; j < stabiliserZ.n_col; ++j) {
                nZ = 0;
                Z_stb_err = error_table[error_distr(gen)];
                int k = 1;
                for (std::array<int,2> pos: pos_array){
                    if (isError(code(2*i+1+pos[0],2*j+1+pos[1]), X_ERROR)) nZ++;
                    code(2*i+1+pos[0],2*j+1+pos[1]) = errorComposite(code(2*i+1+pos[0],2*j+1+pos[1]), Z_stb_err[k]);
                    k++;
                }
                stabiliserZ(i,j) = (nZ+Z_stb_err[0])%2;
            }
        }

    }

//    void stabiliserUpdate(double error_prob){
//
//        std::array<double,512> prob_array;
//        std::array<std::array<int, 5>,512> error_table;
//        std::vector<std::array<int, 2>> pos_array =  {{0,1}, {0,(-1)}, {1,0}, {-1,0}};
//
//        std::ifstream inFile;
//
//        //If file can't be read, maybe the filename buffer is NOT LONG ENOUGH!!!
//        char filename [100];
////        sprintf (filename, "../ParityCheckErrorTable/Data/exact_error_table%.4f.txt", error_prob);
//        sprintf (filename, "../ParityCheckErrorTable/Data/fowler_error_table_X%.4f.txt", error_prob);
//        inFile.open(filename);
//        if (! inFile) {
//            std::cerr << "unable to open file for reading" << std::endl;
//        }
//        for (int i = 0; i < 512; i++) {
//            inFile >> prob_array[i];
//            for (int j = 0; j < 5; j++) {
//                inFile >> error_table[i][j];
//            }
//        }
//        inFile.close();
//
//        std::random_device rd;
//        std::mt19937 gen(rd());
//        std::discrete_distribution<> error_distr(prob_array.begin(), prob_array.end());
//        std::array<int, 5> X_stb_err;
//        int nX, nZ;
//        for (int i = 0; i < stabiliserX.n_row; ++i) {
//            for (int j = 0; j < stabiliserX.n_col; ++j) {
//                nX = 0;
//                //randomly select one row in the error_table
//                X_stb_err = error_table[error_distr(gen)];
//                int k = 1;
//                for (std::array<int,2> pos: pos_array){
//                    //calculate real parity
//                    if (isError(code(2*i+pos[0],2*j+pos[1]), Z_ERROR)) nX++;
//                    //add errors to data qubit
//                    code(2*i+pos[0],2*j+pos[1]) = errorComposite(code(2*i+pos[0],2*j+pos[1]), X_stb_err[k]);
//                    k++;
//                }
//                //add readout errors
//                stabiliserX(i,j) = (nX+X_stb_err[0])%2;
//            }
//        }
//
//        std::ifstream inFile2;
//        //If file can't be read, maybe the filename buffer is NOT LONG ENOUGH!!!
//        char filename2 [100];
////        sprintf (filename2, "../ParityCheckErrorTable/Data/exact_error_table%.4f.txt", error_prob);
//        sprintf (filename2, "../ParityCheckErrorTable/Data/fowler_error_table_Z%.4f.txt", error_prob);
//        inFile2.open(filename2);
//        if (! inFile2) {
//            std::cerr << "unable to open file for reading" << std::endl;
//        }
//        for (int i = 0; i < 512; i++) {
//            inFile2 >> prob_array[i];
//            for (int j = 0; j < 5; j++) {
//                inFile2 >> error_table[i][j];
//            }
//        }
//        inFile2.close();
//        std::discrete_distribution<> error_distr2(prob_array.begin(), prob_array.end());
//        std::array<int, 5> Z_stb_err;
//        for (int i = 0; i < stabiliserZ.n_row; ++i) {
//            for (int j = 0; j < stabiliserZ.n_col; ++j) {
//                nZ = 0;
//                Z_stb_err = error_table[error_distr2(gen)];
//                int k = 1;
//                for (std::array<int,2> pos: pos_array){
//                    if (isError(code(2*i+1+pos[0],2*j+1+pos[1]), X_ERROR)) nZ++;
//                    code(2*i+1+pos[0],2*j+1+pos[1]) = errorComposite(code(2*i+1+pos[0],2*j+1+pos[1]), Z_stb_err[k]);
//                    k++;
//                }
//                stabiliserZ(i,j) = (nZ+Z_stb_err[0])%2;
//            }
//        }
//    }



    //when mode = 1, we use the full error model for the parity check circuit.
    void timeStep(double error_prob, int mode = 1){
        if (mode == 1){
            stabiliserUpdate(error_prob);
        }
        else{
            data.induceError(error_prob);
            stabiliserUpdate();
            stabiliserX.induceError(error_prob);
            stabiliserZ.induceError(error_prob);
        }
        stabiliserX.t +=1;
        stabiliserX.addFlipLoc();
        stabiliserX.last_code = stabiliserX._code;
        stabiliserZ.t +=1;
        stabiliserZ.addFlipLoc();
        stabiliserZ.last_code = stabiliserZ._code;
    }

    void lastStep(double error_prob, int mode = 1){
        if (mode ==1){
            stabiliserUpdate(error_prob);
        }
        else{
            data.induceError(error_prob);
        }
        stabiliserUpdate();
        stabiliserX.t +=1;
        stabiliserX.addFlipLoc();
        stabiliserX.last_code = stabiliserX._code;
        stabiliserZ.t +=1;
        stabiliserZ.addFlipLoc();
        stabiliserZ.last_code = stabiliserZ._code;
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
        int n_error = stabiliser->flip_locs.size();
        std::vector <int> error_label;
        std::array <int, 2> start, loc, min_path;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::binomial_distribution<> direction(1, 0.5);


        int  ver_sign, hor_sign, ver_steps, hor_steps;
        for (int k = 0; k < n_error; ++k) error_label.push_back(k);
        std::set<int> error_corrected; //std::array are implicitly copied.
        for (int error: error_label){
            if (error_corrected.find(error) == error_corrected.end()) { // equivalent to if error is not in error_corrected
                int paired_error = pm->GetMatch(error);
                min_path = stabiliser->minSpatialPath(stabiliser->flip_locs[error],
                                               stabiliser->flip_locs[paired_error]);
                ver_sign = getSign(min_path[0]); //+1 if end[0]>start[0], -1 if otherwise.
                hor_sign = getSign(min_path[1]); //+1 if end[0]>start[0], -1 if otherwise.

                start = {2* stabiliser->flip_locs[error][0] + offset, 2* stabiliser->flip_locs[error][1] + offset};
                loc = start; //implicit copy of start

                //We might want to use step_sign *loc[0] < step_sign * end[0] condition to substitue ver_steps. But the
                //situation is much more complicated when we need to cross the boundary.
                ver_steps = 0;
                hor_steps = 0;
                // ////Comment out this section if we don't want error correction along random min path.
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
                error_corrected.insert(paired_error); //we did not add error because error is already iterated over.
            }
        }
        delete pm;
    }

    void fixError(){
        fixError(X_STB);
        fixError(Z_STB);
    }

    /* For X_stb, we count along the row 1, for which Z_errors cannot run along the row (unlike row 0
     * for which Z_errors can run along the row going through X_stb.
     * Similarly for Z_stb
     * */

    bool hasLogicalError(StabiliserType stb){
        int n_v_error = 0;
        int n_h_error = 0;
        if (stb == X_STB){
            for (int j = 0; j < n_col; j += 2) {
                if (isError(code(1,j), Z_ERROR)) n_v_error++; // check horizontally if any vertical error cuts through
            }
            for (int i = 0; i < n_row; i += 2) {
                if (isError(code(i,1), Z_ERROR)) n_h_error++;
            }
        }
        else if (stb == Z_STB){
            for (int j = 1; j < n_col; j += 2) {
                if (isError(code(0,j), X_ERROR)) n_v_error++;
            }
            for (int i = 1; i < n_row; i += 2) {
                if (isError(code(i,0), X_ERROR)) n_h_error++;
            }
        }
        return ((n_h_error%2) or (n_v_error%2));
    }

    bool hasLogicalError(){
        return hasLogicalError(X_STB) or hasLogicalError(Z_STB);
    }
};

//In this function we assume the number of time steps is the same as length L.
double averageLogicalError(int L, double data_error_rate, int n_runs, int error_mode){
    assert(L%2 == 0);
    int logical_errors_counter = 0;
    for (int i = 0; i < n_runs; ++i) {
        ToricCode c(L, L);
        for (int t = 0; t < L/2-1; ++t) {
            c.timeStep(data_error_rate, error_mode);
        }
        c.lastStep(data_error_rate, error_mode);
        c.fixError();
        logical_errors_counter += c.hasLogicalError();
    }
    return (double)logical_errors_counter/(double)n_runs;
}

void errorDataOutput(int n_runs, int error_mode, int job_array_id = 1){
    std::vector<double> data_error_rate_array;
    for (double data_error_rate = 0.001; data_error_rate <0.011; data_error_rate += 0.001) {
        data_error_rate_array.push_back(data_error_rate);
    }
    std::vector<int> code_size_array;
    for (int code_size = 12; code_size <= 24; code_size += 4) {
        code_size_array.push_back(code_size);
    }
    std::ofstream file;
    char filename [100];
    if(error_mode == 1){
        sprintf (filename, "../data_files/FullCircuitErrorCumulativeErrorData%d.txt", job_array_id);
    }
    else{
        sprintf (filename, "../data_files/CumulativeErrorData%d.txt", job_array_id);
    }
    double avg_log_error;
    for (double data_error_rate : data_error_rate_array) {
        for (int code_size: code_size_array){
//            printf("calculating avg log errors for code size %d with data error rate %.3f\n", code_size, data_error_rate);
//            fflush(stdout);
            avg_log_error = averageLogicalError(code_size, data_error_rate, n_runs, error_mode);
            std::cout<<data_error_rate<<","<<code_size<<","<<avg_log_error<<","<<n_runs<<std::endl;
            file.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);
            file<<data_error_rate<<","<<code_size<<","<<avg_log_error<<","<<n_runs<<std::endl;
            // code_size is the size of stabiliser grid, the size of the
            file.close();
        }
    }
}

int main() {
//    load_cerr_table();
    int cbit, nbit;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::tie(cbit, nbit) = pass_through_cnot(i,j);
            printf("%d %d -> %d %d \n", i,j,cbit,nbit);
        }
    }



//    errorDataOutput(500, 1, 100);
}





//int main(int argc, char *argv[]) {
//    int job_array_id = std::atoi(argv[1]);
//    errorDataOutput(1000, 1, job_array_id);
//    return 0;
//}
