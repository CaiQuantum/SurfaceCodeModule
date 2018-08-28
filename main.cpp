#include <utility>

#include <Python.h>
#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <set>
#include <PerfectMatching.h>
#include <fstream>
#include "ErrorFunc.h"
#include "HelperFunc.h"
#include <tuple>
#include <stdexcept>

/*To make perfectmatching.h to work, we need to first delete example.cpp in blossom_dir, then we also cannot use the
 * triangle package as suggested due to lack of X11. We use import project in Clion to rewrite Cmake. Remember to exclude
 * the unwanted files.
 *
 * There are two parameters you can tweak in the final version of this program. The first one is the time-spatial relative
 * weight in perfectmatching. i.e. how important is time distance relative to spatial distance when we are calculating
 * the perfect matching distance.
 * */

/* To do list:
 * Think about how to merge the Hadamard over all data qubits in sym stb update into error tables.
 *
 * Implementation of planar code:
 *
 * Add extra two row and two column to the toric code.
 *
 * */



std::random_device rd;
std::mt19937 rand_gen(rd());
const int TD_DISTANCE_RATIO = 1;

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
        // Here we assume we will only call row that are in  -n_row < row < n_row
        // Otherwise we need to use row = (row% n_row + n_row)%n_row;
        row = (row + n_row)%n_row;
        col = (col + n_col)%n_col;
        return _code[row][col];
    }

    int& operator()(int row, int col){
        row = (row + n_row)%n_row;
        col = (col + n_col)%n_col;
        return _code[row][col];
    }

    void printCode(){
        printMatrix(_code);
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

        std::binomial_distribution<> error_occur(1, error_prob);

        for(int i = 0; i < n_row; i++){
            for(int j = 0; j < n_col; j++){
                if (error_occur(rand_gen)){
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
    std::vector<std::vector<int>> last_code;
    std::discrete_distribution<> error_distr;
    std::array<int,2> loc_in_toric;
    int planar;
    int t = 0;
    int n_error=0; //number of errors in stb excluding boundary stb errors

public:
    Stabiliser(int n_row, int n_col, StabiliserType stabiliser_type, std::array<int,2> loc_in_toric, int planar=0):
            Code(n_row, n_col), stabiliser_type(stabiliser_type), loc_in_toric(loc_in_toric), planar(planar){
        std::array<double,512> error_prob{};
        error_prob[0] = 1;
        std::discrete_distribution<> error_distr(error_prob.begin(), error_prob.end());
        last_code = _code;
    }
    void induceError(double error_prob){
        assert(error_prob <= 1);
        std::binomial_distribution<> error_occur(1, error_prob);

        for(int i = 0; i < n_row - loc_in_toric[0] * planar; i++){
            for(int j = 0; j < n_col - loc_in_toric[1] * planar; j++){
                if (error_occur(rand_gen)){
                    _code[i][j] = not _code[i][j];
                }
            }
        }
    }

    void addFlipLoc(){
        for (int i = 0; i < n_row - loc_in_toric[0] * planar; i++) {
            for (int j = 0; j < n_col - loc_in_toric[1] * planar; j++) {
                if (_code[i][j] != last_code[i][j]) {
                    flip_locs.push_back({i,j,t});
                }
            }
        }
    }

    int minLinePath(const int start, const int end, const int n_L){
        int d_1 = end - start;
        int sign_d1 = (end > start) - (start > end);
        int d_2 = d_1 - sign_d1 * n_L;
        if (abs(d_1) < abs(d_2)) return d_1;
        else return d_2;
    }

    inline std::array<int, 2> minSpatialPath(std::array<int,3> start, std::array<int,3> end){
        return {minLinePath(start[0], end[0], n_row), minLinePath(start[1], end[1], n_col)};
    }

//    inline int spatialDistance(std::array<int,3> loc1, std::array<int,3> loc2){
//        int d0 = abs(loc1[0] - loc2[0]);
//        int d1 = abs(loc1[1] - loc2[1]);
//        int d = std::min(d0, n_row - d0) + std::min(d1, n_col - d1);
//
//        return d;
//    }

    inline int distance(std::array<int,3> loc1, std::array<int,3> loc2){
        // we create d here to reduce the need to access element of loc, hence increasing ths speed
        int d0 = abs(loc1[0] - loc2[0]);
        int d1 = abs(loc1[1] - loc2[1]);
        int d = std::min(d0, n_row - d0) + std::min(d1, n_col - d1) + TD_DISTANCE_RATIO * abs(loc1[2] - loc2[2]);
        return d;
    }

    std::array<int, 3> getBoundaryLoc(std::array<int,3> loc){
        std::array<int, 3> boundary_loc = loc;

        if (loc_in_toric[0] == 1){
            boundary_loc[0] = n_row - 1;
        }
        else{
            boundary_loc[1] = n_col - 1;
        }

        return boundary_loc;
    }


    inline int distanceToBoundary(std::array<int,3> loc){
        if (loc_in_toric[0] == 1){
            // loc[0] + 1 is actually, loc[0] - (-1). -1 and n_row -1 are both index of the boundary.
            // we create d here to reduce the need to access element of loc, hence increasing ths speed
            int d = loc[0];
            return std::min(d + 1, n_row - 1 - d);
        }
        else{
            int d = loc[1];
            return std::min(d + 1, n_col - 1 - d);
        }

//        return spatialDistance(loc, getBoundaryLoc(loc));
    }

    PerfectMatching* getErrorMatching(){
//        assert (flip_locs.size() != 0);
        n_error = flip_locs.size();

        std::vector<std::array<int,3>> edges;
        if (planar){
            // reserving space for all the edges related to mirror errors.
            edges.reserve(n_error + n_error* (n_error-1)/2);
            // Adding bouandary mirror errors
            for (int i = 0; i < n_error; ++i) {
                flip_locs.push_back(getBoundaryLoc(flip_locs[i]));
                edges.push_back({i, (i + n_error), distanceToBoundary(flip_locs[i])});
            }

            // Adding edges between mirror errors
            for (int i = n_error; i < 2 * n_error; ++i) {
                for (int j = n_error; j < i; ++j) {
                    edges.push_back({i, j , 0});
                }
            }

            for (int i = 0; i < n_error; ++i) {
                for (int j = 0; j < i ; ++j) {
                    //add spatial distance and time distance together as cost.
                    int d_link = distance(flip_locs[i], flip_locs[j]);

                    int d_bound = edges[i][2] + edges[j][2];

                    if (d_link < d_bound){
                        edges.push_back({i, j , d_link});
                    }
                }
            }

        }
        else{
            edges.reserve(n_error* (n_error-1)/2);
            for (int i = 0; i < n_error; ++i) {
                for (int j = 0; j < i ; ++j) {
                    //add spatial distance and time distance together as cost.
                    int d = distance(flip_locs[i], flip_locs[j]);
                    edges.push_back({i, j ,d});
                }
            }
        }

        int n_nodes = flip_locs.size();
        int n_edges = edges.size();

        auto *pm = new PerfectMatching(n_nodes, n_edges);
        struct PerfectMatching::Options options;
        options.verbose = false;
        pm->options = options;
        for (std::array<int,3> e: edges){
            pm->AddEdge(e[0], e[1], e[2]);
        }

        pm->Solve();

        return pm;
    }
};

class SurfaceCode{
public:
    int n_row;
    int n_col;
    Data primal_data;
    Data dual_data;
    Stabiliser stabiliserX;
    Stabiliser stabiliserZ;
    int planar;

public:
    SurfaceCode(int n_row, int n_col, std::discrete_distribution<> X_error_distr,
                std::discrete_distribution<> Z_error_distr, int planar=0):
            n_row(n_row + planar), n_col(n_col+planar),
            primal_data((n_row+planar)/2, (n_col+planar)/2),
            dual_data((n_row+planar)/2, (n_col+planar)/2),
            stabiliserX((n_row+planar)/2, (n_col+planar)/2, X_STB, {1, 0}, planar),
            stabiliserZ((n_row+planar)/2, (n_col+planar)/2, Z_STB, {0, 1}, planar),
            planar(planar){
        stabiliserX.error_distr = std::move(X_error_distr);
        stabiliserZ.error_distr = std::move(Z_error_distr);
        assert((n_row+planar)%2==0);
        assert((n_col+planar)%2==0);
    }

    int& code(int row, int col){
        /* The top right corner of out code looks like
         *
         *  primal_d      Z      primal_d      Z
         *      X      dual_d        X      dual_d
         *  primal_d      Z      primal_d      Z
         *      X      dual_d        X      dual_d
         * */

        row = (row + n_row)%n_row;
        col = (col + n_col)%n_col;
        if (row%2 == 0 and col%2 == 0) return primal_data._code[row/2][col/2];
        else if (row%2 == 1 and col%2 == 1) return dual_data._code[(row -1)/2][(col-1)/2];
        else if (row%2 == 0) return stabiliserZ._code[row/2][(col-1)/2];
        else return stabiliserX._code[(row-1)/2][col/2];
    }
    //red: 31, grn: 32, yel: 33, blu: 34, mag: 35, cyn: 36, wht: 37
    void printCode(){
        for (int i = 0; i < primal_data.n_row * 2; i++) {
            for (int j = 0; j < primal_data.n_col; j++) {
                if (i%2 == 0) {
                    printf("%2d ", primal_data(i/2, j));
                    printf("\x1B[34m%2d\x1B[0m ", stabiliserZ(i/2, j));
                }
                else {
                    printf("\x1B[31m%2d\x1B[0m ", stabiliserX((i-1)/2, j));
                    printf("%2d ", dual_data((i-1)/2, j));
                }
            }
            printf("\n");
        }
        printf("\n");
    }

    void perfectStabiliserUpdate(){
        // Cannot change the order, not even to {{1, 0}, {0,(-1)}, {0, 1}, {-1, 0}}
        std::array<std::array<int, 2>, 4> pos_array { {{0,1}, {0,(-1)}, {1, 0}, {-1, 0}} };
        for (int stb = 0; stb < 2; ++stb)  {
            int tracked_error;
            Stabiliser *stabiliser;

            if (stb == Z_STB) {
                stabiliser = &stabiliserZ;
                tracked_error = X_ERROR;
            }
            else {
                stabiliser = &stabiliserX;
                tracked_error = Z_ERROR;
            }

            std::array<int, 2> pos_offset = stabiliser->loc_in_toric;
            for (int i = 0; i < (stabiliser->n_row - pos_offset[0] * planar); ++i) {
                for (int j = 0; j < (stabiliser->n_col - pos_offset[1] * planar); ++j) {
                    int nP = 0;
                    for (std::array<int, 2> pos: pos_array) {
                        int &data_qubit = code( 2*i+pos_offset[0]+pos[0], 2*j+pos_offset[1]+pos[1]);
                        //calculate real parity
                        if (isError(data_qubit, tracked_error)) nP++;
                    }
                    stabiliser->code(i, j) = nP % 2;
                }
            }
        }
    }

    void reset_boundary(){
        assert(planar);
        for (int i = 0; i < dual_data.n_row; i++){
            dual_data(i, dual_data.n_col-1) = NO_ERROR;
        }
        for (int j = 0; j < dual_data.n_col; j++){
            dual_data(dual_data.n_row - 1, j) = NO_ERROR;
        }
    }

    /*
     * sym = 0: two hadamard in X_stb ancilla, but no Hadamard in Z_stb ancilla
     * sym = 1: hadamard applied to whole qubit array (NOT within parity measurement circuit) before X/Z measurements
     */
    void stabiliserUpdate(){

        //First row is order of X_stb update, second row is Z_stb.
        //About the number of braces see https://stackoverflow.com/questions/43628497/3d-stdarray-in-c
        std::array<std::array<std::array<int, 2>, 4>, 2> pos_array{{
                                                                           {{{0, 1}, {0, (-1)}, {1, 0}, {-1, 0}}},
                                                                           {{{0, 1}, {0, (-1)}, {1, 0}, {-1, 0}}}
                                                                   }};

        for (int stb = 0; stb < 2; ++stb)  {
            int tracked_error;
            Stabiliser *stabiliser;

            if (stb == Z_STB) {
                stabiliser = &stabiliserZ;
                tracked_error = X_ERROR;
            }
            else {
                stabiliser = &stabiliserX;
                tracked_error = Z_ERROR;
            }

            std::array<int, 2> pos_offset = stabiliser->loc_in_toric;
            for (int i = 0; i < (stabiliser->n_row - pos_offset[0] * planar); ++i) {
                for (int j = 0; j < (stabiliser->n_col - pos_offset[1] * planar); ++j) {
                    int nP = 0;
                    //randomly select one row in the error_table
                    std::array<int, 5> stb_err = error_table[stabiliser->error_distr(rand_gen)];
                    int k = 1;
                    for (std::array<int, 2> pos: pos_array[stb]) {
                        int d_i = 2*i+pos_offset[0]+pos[0];
                        int d_j = 2*j+pos_offset[1]+pos[1];
                        if (planar and (d_i == -1 or d_j == -1)){ continue;}
                        int &data_qubit = code(d_i, d_j);
                        //calculate real parity
                        if (isError(data_qubit, tracked_error)) nP++;
                        //add errors to data qubit
                        data_qubit = errorComposite(data_qubit, stb_err[k]);
                        k++;
                    }
                    //add readout errors
                    stabiliser->code(i, j) = (nP + stb_err[0]) % 2;
                }
            }
        }
        if (planar){
            reset_boundary();
        }
    }

    void syncStabiliserUpdate(){
        //First row is order of X_stb update, second row is Z_stb.
        //About the number of braces see https://stackoverflow.com/questions/43628497/3d-stdarray-in-c
        std::array<std::array<std::array<int, 2>, 4>, 2> pos_array{{
                                                                           {{{0, -1}, {-1, 0}, {1, 0}, {0, 1}}},
                                                                           {{{0, -1}, {-1, 0}, {1, 0}, {0, 1}}}
                                                                   }};
        std::array<std::vector<std::vector<int>>, 2> stb_error_array;
        std::array<Stabiliser*, 2>stabiliser {{&stabiliserX, &stabiliserZ}};
        std::array<DataError, 2>tracked_error{Z_ERROR, X_ERROR};

        for (int stb = 0; stb < 2; ++stb) {
            //Initialisation
            std::array<int, 2> pos_offset = stabiliser[stb]->loc_in_toric;
            stb_error_array[stb].resize(stabiliser[stb]->n_row - pos_offset[0] * planar);
            for (int i = 0; i < (stabiliser[stb]->n_row - pos_offset[0] * planar); ++i) {
                stb_error_array[stb][i].resize(stabiliser[stb]->n_col - pos_offset[1] * planar);
                for (int j = 0; j < (stabiliser[stb]->n_col - pos_offset[1] * planar); ++j) {
                    stabiliser[stb]->code(i, j) = 0;
                    //randomly select one row in the error_table
                    stb_error_array[stb][i][j] = stabiliser[stb]->error_distr(rand_gen);
                }
            }
        }

        for (int qb = 0; qb < 4; ++qb) {
            for (int stb = 0; stb < 2; ++stb) {
                std::array<int, 2> pos_offset = stabiliser[stb]->loc_in_toric;
                int starting_i = 0;
                int starting_j = 0;
                if (planar and pos_offset[0]==0 and pos_array[stb][qb][0] == -1){starting_i = 1;}
                if (planar and pos_offset[1]==0 and pos_array[stb][qb][1] == -1){starting_j = 1;}
                for (int i = starting_i; i < (stabiliser[stb]->n_row - pos_offset[0] * planar); ++i) {
                    for (int j = starting_j; j < (stabiliser[stb]->n_col - pos_offset[1] * planar); ++j) {
                        int &data_qubit = code(2 * i + pos_offset[0] + pos_array[stb][qb][0],
                                               2 * j + pos_offset[1] + pos_array[stb][qb][1]);
                        //calculate real parity
                        if (isError(data_qubit, tracked_error[stb])) stabiliser[stb]->code(i, j)++;
                        //add errors to data qubit
                        data_qubit = errorComposite(data_qubit, error_table[stb_error_array[stb][i][j]][qb + 1]);
                    }
                }
            }
        }

        for (int stb = 0; stb < 2; ++stb)  {
            //Initialisation
            std::array<int, 2> pos_offset = stabiliser[stb]->loc_in_toric;
            for (int i = 0; i < (stabiliser[stb]->n_row - pos_offset[0] * planar); ++i) {
                for (int j = 0; j < (stabiliser[stb]->n_col - pos_offset[1] * planar); ++j) {
                    stabiliser[stb]->code(i, j) = (stabiliser[stb]->code(i, j) + error_table[stb_error_array[stb][i][j]][0]) % 2;
                }
            }
        }
        if (planar){
            reset_boundary();
        }
    }


    // we use error table for asymmetric parity check circuit (in Fowler's paper).
    void timeStep(bool is_last_step = true){

        //In the last step, the data qubit are measured, hence, we can know the exact parity of the data.
        if (is_last_step) perfectStabiliserUpdate();
        else stabiliserUpdate();

        stabiliserX.t +=1;
        stabiliserX.addFlipLoc();
        stabiliserX.last_code = stabiliserX._code;
        stabiliserZ.t +=1;
        stabiliserZ.addFlipLoc();
        stabiliserZ.last_code = stabiliserZ._code;
    }

    // we induce random errors to qubits (ignore the circuit structure)
    void randomErrorTimeStep(double error_prob, bool is_last_step = false){

        primal_data.induceError(error_prob);
        dual_data.induceError(error_prob);
        perfectStabiliserUpdate();
        stabiliserX.induceError(error_prob);
        stabiliserZ.induceError(error_prob);

        //In the last step, the data qubit are measured, hence, we can know the exact parity of the data.
        if (is_last_step) perfectStabiliserUpdate();

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
        if (stabiliser_type == X_STB){
            stabiliser = &stabiliserX;
            ERROR = Z_ERROR;
        }
        else{
            stabiliser = &stabiliserZ;
            ERROR = X_ERROR;
        }
        std::array<int, 2> offset = stabiliser->loc_in_toric;
        PerfectMatching* pm = stabiliser->getErrorMatching();
//        int n_error = stabiliser->n_error;
//        std::array <int, 2> start, loc, min_path;

//        std::binomial_distribution<> direction(1, 0.5);


//        int  ver_sign, hor_sign, ver_steps, hor_steps;
        std::set<int> error_corrected; //std::array are implicitly copied.
        for (int error = 0; error < stabiliser->n_error; ++error){
            if (error_corrected.find(error) == error_corrected.end()) { // equivalent to if error is not in error_corrected
                int paired_error = pm->GetMatch(error);
                std::array <int, 2> min_path = stabiliser->minSpatialPath(stabiliser->flip_locs[error],
                                                                          stabiliser->flip_locs[paired_error]);
                int ver_sign = getSign(min_path[0]); //+1 if end[0]>start[0], -1 if otherwise.
                int hor_sign = getSign(min_path[1]); //+1 if end[0]>start[0], -1 if otherwise.

                std::array <int, 2> start = {2 * stabiliser->flip_locs[error][0] + offset[0], 2 * stabiliser->flip_locs[error][1] + offset[1]};
                std::array <int, 2> loc = start; //implicit copy of start

                //We might want to use step_sign *loc[0] < step_sign * end[0] condition to substitue ver_steps. But the
                //situation is much more complicated when we need to cross the boundary.
                int ver_steps = 0;
                int hor_steps = 0;
                // ////Comment out this section if we don't want error correction along random min path.
//                while (ver_steps < abs(min_path[0]) and hor_steps < abs(min_path[1])){
//                    if (direction(rand_gen)){
//                        code(loc[0] + ver_sign, loc[1]) = errorComposite(code(loc[0] + ver_sign, loc[1]), ERROR);
//                        loc[0] += 2*ver_sign;
//                        ver_steps +=1;
//                    }
//                    else{
//                        code(loc[0], loc[1] + hor_sign) = errorComposite(code(loc[0], loc[1] + hor_sign), ERROR);
//                        loc[1] += 2*hor_sign;
//                        hor_steps +=1;
//                    }
//                }
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

    double hasLogicalError(){
        bool primal_X_error = false;
        bool primal_Z_error = false;
        // Checking X error of primal lattice (cutting through Z_stb):
        for (int i = 0; i < primal_data.n_row; i ++) {
            if (isError(primal_data(i, 0), X_ERROR)) primal_X_error = not primal_X_error;
        }
        // Checking Z error of primal lattice (cutting through X_stb):
        for (int j = 0; j < primal_data.n_col; j ++) {
            if (isError(primal_data(0, j), Z_ERROR)) primal_Z_error = not primal_Z_error;
        }
        auto primal_error = (double)(primal_X_error or primal_Z_error);

        if (not planar) {
            bool dual_X_error = false;

            bool dual_Z_error = false;
            // Checking X error of primal lattice (cutting through Z_stb):
            for (int j = 0; j < dual_data.n_col; j ++) {
                if (isError(dual_data(0, j), X_ERROR)) dual_X_error = not dual_X_error;
            }
            // Checking Z error of primal lattice (cutting through X_stb):
            for (int i = 0; i < dual_data.n_row; i++) {
                if (isError(dual_data(i, 0), Z_ERROR)) dual_Z_error = not dual_Z_error;
            }
            auto dual_error = (double)(dual_X_error or dual_Z_error);
            return (primal_error + dual_error) / 2;
        }
        else{
            return primal_error;
        }
    }
};


//In this function we assume the number of time steps is the same as length L.
double averageLogicalError(int L, int n_runs, const char* main_file_name, int planar, double time_step_ratio = 0.5){
//    assert(L%2 == 0);
    double logical_errors_counter = 0;
    std::array<std::array<double,512>, 2> error_prob{};

    for (int stb = 0; stb < 2; ++stb) {
        std::ifstream inFile;
        char error_table_filename [100];
        if (stb == X_STB) {
            sprintf(error_table_filename, "%s_X.txt", main_file_name);
        }
        else{
            sprintf(error_table_filename, "%s_Z.txt", main_file_name);
        }
        inFile.open(error_table_filename);
        if (!inFile) {
            char error_message[150];
            sprintf(error_message, "Could not open %s", error_table_filename);
            throw std::runtime_error(error_message);
        }
        for (int i = 0; i < 512; i++) {
            inFile >> error_prob[stb][i];
            inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        inFile.close();
    }
    std::discrete_distribution<> X_error_distr(error_prob[0].begin(), error_prob[0].end());
    std::discrete_distribution<> Z_error_distr(error_prob[1].begin(), error_prob[1].end());

    for (int i = 0; i < n_runs; ++i) {
        SurfaceCode surface_code(L, L, X_error_distr, Z_error_distr, planar);
//        surface_code.printCode();
        for (int t = 0; t < (L + planar) * time_step_ratio; ++t) {
            surface_code.timeStep(false);
        }
        surface_code.timeStep(true); // the last argument true means it is the last step
        surface_code.fixError();
        logical_errors_counter += surface_code.hasLogicalError();
    }
    return logical_errors_counter/(double)n_runs;
}


static PyObject * _averageLogicalError(PyObject *self, PyObject *args) {
    int L;
    int n_runs;
    const char* main_file_name;
    double res;
    int planar;
    double time_step_ratio;
    if (!PyArg_ParseTuple(args, "iisid", &L, &n_runs, &main_file_name, &planar, &time_step_ratio))
        return nullptr;
    res = averageLogicalError(L, n_runs, main_file_name, planar, time_step_ratio);
    return PyFloat_FromDouble(res);
}


double averagePhysicalError(int L, int n_runs, const char* main_file_name, int planar, double time_step_ratio = 0.5){
//    assert(L%2 == 0);
    int physical_errors_counter = 0;
    std::array<std::array<double,512>, 2> error_prob{};

    for (int stb = 0; stb < 2; ++stb) {
        std::ifstream inFile;
        char error_table_filename [100];
        if (stb == X_STB) {
            sprintf(error_table_filename, "%s_X.txt", main_file_name);
        }
        else{
            sprintf(error_table_filename, "%s_Z.txt", main_file_name);
        }
        inFile.open(error_table_filename);
        if (!inFile) {
            char error_message[150];
            sprintf(error_message, "Could not open %s", error_table_filename);
            throw std::runtime_error(error_message);
        }
        for (int i = 0; i < 512; i++) {
            inFile >> error_prob[stb][i];
            inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        inFile.close();
    }
    std::discrete_distribution<> X_error_distr(error_prob[0].begin(), error_prob[0].end());
    std::discrete_distribution<> Z_error_distr(error_prob[1].begin(), error_prob[1].end());

    for (int i = 0; i < n_runs; ++i) {
        SurfaceCode surface_code(L, L, X_error_distr, Z_error_distr, planar);
//        surface_code.printCode();
        for (int t = 0; t < (L + planar) * time_step_ratio; ++t) {
            surface_code.timeStep(false);
        }
        surface_code.timeStep(true); // the last argument true means it is the last step
        physical_errors_counter += surface_code.stabiliserX.flip_locs.size();
        physical_errors_counter += surface_code.stabiliserZ.flip_locs.size();
    }
    return (double)physical_errors_counter/(double)n_runs;
}


static PyObject * _averagePhysicalError(PyObject *self, PyObject *args) {
    int L;
    int n_runs;
    const char* main_file_name;
    double res;
    int planar;
    double time_step_ratio;
    if (!PyArg_ParseTuple(args, "iisid", &L, &n_runs, &main_file_name, &planar, &time_step_ratio))
        return nullptr;
    res = averagePhysicalError(L, n_runs, main_file_name, planar, time_step_ratio);
    return PyFloat_FromDouble(res);
}

static PyMethodDef SurfaceCodeMethods[] = {
        {
                "average_logical_error",
                _averageLogicalError,
                METH_VARARGS,
                ""
        },
        {
                "average_physical_error",
                _averagePhysicalError,
                METH_VARARGS,
                ""
        },
        {nullptr, nullptr, 0, nullptr}
};


static struct PyModuleDef SurfaceCodeModule = {
        PyModuleDef_HEAD_INIT,
        "SurfaceCode",
        nullptr,
        -1,
        SurfaceCodeMethods
};


PyMODINIT_FUNC PyInit_SurfaceCode(void)
{
    PyObject *m;
    m = PyModule_Create(&SurfaceCodeModule);
    if (m == nullptr)
        return nullptr;
    printf("init SurfaceCode module\n");
    return m;
}


//int main() {
//    const char* main_file_name = "../../../ElectronShuttling/Data/ErrorTable/error_table_0.002_0.000_0.005";
//    double error = averageLogicalError(11, 10000, main_file_name, true, 4);
//    std::cout<< error << std::endl;
//}

