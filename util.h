#ifndef MUTIL_H
#define MUTIL_H

#include <cmath>
#include <string>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

struct data {
    data(double num, int ind1, int ind2) {
        number = num;
        index1 = ind1;
        index2 = ind2;
    }
    double number;
    int index1;
    int index2;
};

struct by_number {
    bool operator()(data const &left, data const &right) {
        return abs(left.number) > abs(right.number);
    }
};

int safe_mmap_read_only(string filename, char **ptr_mmaped_file, size_t *ptr_file_size);

/*
 convert int to string
 */
string convertInt(int number);

/*
 find n factorial
*/
long int fact(int n) ;

/*
 find minimum between a and b
 */
double min(double a, double b) ;

/*
 find combinration n choose r
 */
long int nCr(int n, int r) ;

void printVec(vector<int> v);
void printCharVec(vector<char> v);


/*
 * import data from snp map file
 */
void importSnpMap(string snpMapFile, int numCols, vector<string> * firstCol, vector<vector<int>> * remainingCols);


/*
 import data from file
 */
void importData(string fileName, vector<double> *& vector);

/*
 import second column of data from file
 */
void importDataSecondColumn(string fileName, vector<double>& vector);

/*
 import first column of data from file
 */
void importDataFirstColumn(string fileName, vector<string>& list, int ignore=0);

/*
 export data type char from file
 */
void exportVector2File(string fileName, char * data, int size);

/*
 export data type double from file
 */
void exportVector2File(string fileName, double * data, int size);

/*
 export data type int from file
 */
void export2File(string fileName, int data);

/*
 export data type double from file
 */
void export2File(string fileName, double data);

/*
 make a matrix semi-positive definite, ie. full rank
 */
void makeSigmaPositiveSemiDefinite(mat * sigma, int size);

/*
 handles low rank case, performs eigenvalue decomposition on the matrix sigma = LDL(trans)
 On output the diagonal of the input matrix sigma stores the diagonal elements of D, 
 and the lower triangular portion of A contains the matrix L. 
 Since L has ones on its diagonal these do not need to be explicitely stored. 
 The upper triangular portion of A is unmodified.
*/
mat* eigen_decomp(mat* sigma, int size);
#endif
