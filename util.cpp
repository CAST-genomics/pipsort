#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <armadillo>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>


using namespace std;
using namespace arma;


int safe_mmap_read_only(string filename, char **ptr_mmaped_file, size_t *ptr_file_size) {
  int fd;
  struct stat filestat;
  size_t len;
  fd = open(filename.c_str(), O_RDONLY);
  if ( fd < 0 ) {
    printf("Could not open %s\n", filename.c_str());
    return -1;
  }
  int status = fstat(fd, &filestat);
  if ( status < 0 ) {
    return -1;
  }
  len = filestat.st_size;
  if ( len == 0 ) {
    *ptr_mmaped_file = NULL;
    *ptr_file_size = 0;
  } else {
    *ptr_mmaped_file = (char *)mmap(NULL, (size_t) len, PROT_READ, MAP_SHARED, fd, 0);
  }
  close(fd);
  *ptr_file_size = (size_t) len;
  return 0;
}

long int fact(int n) {
    if(n==0)
        return 1;
    return n* fact(n-1);
}

long int nCr(int n, int r) {
    long int result = 1;
    for(int i = n; i > n-r; i--)
        result *= i;
    return result/fact(r);
}

double min(double a, double b) {
    if(a>b)
        return b;
    else
        return a;
}

void printVec(vector<int> v) {
   for ( int i = 0; i < v.size(); i++ ) {
     printf("%d", v[i]);
   }
   printf("\n");
}
void printCharVec(vector<char> v) {
   for ( int i = 0; i < v.size(); i++ ) {
     printf("%c", v[i]);
   }
   printf("\n");
}



void importData(string fileName, vector<double> *& vector) {
    ifstream file(fileName.c_str(), ifstream::in);
    if (!file) {
        cout << "Unable to open file; This is why";
        exit(1); // terminate with error
    }

    double data;
    while(file >> data){ vector->push_back(data); }
    file.close();
}


void importSnpMap(string snpMapFile, int numCols, vector<string> * firstCol, vector<vector<int>> * remainingCols) {
    fstream f;
    f.open(snpMapFile, ios::in);
    if ( !f.is_open()) {
       cout << "Could not open file\n";
       exit(1);
    }

    vector<string> row;
    string line, word;
  
    while (getline(f, line)) {
  
        row.clear();
  
        stringstream s(line);
  
	for ( int i = 0; i < numCols; i++ ) {
            getline(s, word, ',');
	    if ( i == 0 ) {
                firstCol->push_back(word); //rsid
	    } else {
                (*remainingCols)[i-1].push_back(stoi(word)); //idx in study
	    }
	}
    }  
    f.close();
}


/*
 The column index starts by 1 in this implemenation
 */
void importDataSecondColumn(string fileName, vector<double>& vector) {
    string line = "";
    string dataS = "";
    double data = 0.0;
    ifstream fin(fileName.c_str(), std::ifstream::in);
    while( getline(fin, line) ){
        istringstream iss(line);
        iss >> dataS;
        iss >> data;
        vector.push_back((double)data);
    }
    fin.close();
}

void importDataFirstColumn(string fileName, vector<string>& list, int ignore=0) {
    string data = "";
    string line = "";
    ifstream fin(fileName.c_str(), std::ifstream::in);
    for(int i = 0; i < ignore; i++)
        getline(fin, line);

    while( getline(fin, line) ){
        istringstream iss(line);
        iss >> data;
        list.push_back(data);
    }
    fin.close();
}

string convertInt(int number) {
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

void exportVector2File(string fileName, char * data, int size) {
    ofstream outfile(fileName.c_str(), ios::out | ios::app);
    for (int i = 0; i < size; i++)
        outfile << data[i] << " ";
    //outfile << endl;
    outfile.close();
}

void exportVector2File(string fileName, double * data, int size) {
    ofstream outfile(fileName.c_str(), ios::out | ios::app);
    for (int i = 0; i < size; i++)
        outfile << data[i] << " ";
    //outfile << endl;
    outfile.close();
}

void export2File(string fileName, double data) {
    ofstream outfile(fileName.c_str(), ios::out | ios::app);
    outfile << data << endl;
    outfile.close();
}

void export2File(string fileName, int data) {
    ofstream outfile(fileName.c_str(), ios::out | ios::app);
    outfile << data << endl;
    outfile.close();
}

void makeSigmaPositiveSemiDefinite(mat* sigma, int size) {
    int gsl_tmp = 0;
    double matDet  = 0;
    double addDiag = 0;
    bool positive = false;

    //gsl_set_error_handler_off();
    gsl_matrix * tmpResultMatrix = gsl_matrix_calloc (size, size);
    gsl_permutation *p = gsl_permutation_alloc(size);
    do{
        for(int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if(i==j)
                    gsl_matrix_set(tmpResultMatrix,i,j,(*sigma)(i, j)+addDiag);
                else
                    gsl_matrix_set(tmpResultMatrix,i,j,(*sigma)(i, j));
            }
        }

        gsl_linalg_LU_decomp(tmpResultMatrix, p, &gsl_tmp);
        matDet = gsl_linalg_LU_det(tmpResultMatrix,gsl_tmp);
        if(matDet > 0)
            positive = true;
        else {
            addDiag += 0.01;
        }
    } while(!positive);

    for(int i = 0; i < size; i++){
        (*sigma)(i,i) = (*sigma)(i,i) + addDiag;
    }
}

mat* eigen_decomp(mat* sigma, int size){
    gsl_matrix *tmpSigma = gsl_matrix_calloc(size, size);
    gsl_matrix *eigenvec = gsl_matrix_calloc(size, size);
    gsl_vector *eigenval = gsl_vector_calloc(size);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(size);
    //gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(size);

    for(int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            gsl_matrix_set(tmpSigma,i,j,(*sigma)(i, j));
        }
    }

    //decompose to sigma = QDQ(trans), Q is orthogonal matrix of eigenvectors, D is diagonal matrix of eigenvalues
    gsl_eigen_symmv(tmpSigma, eigenval, eigenvec, w);

    //store all the eigenvalues
    mat* vals = new mat(size, size,fill::zeros);
    for(int i = 0; i < size; i++){
        (*vals)(i,i) = gsl_vector_get(eigenval,i);
    }

    //store all the eigenvectors in sigma
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            (*sigma)(i,j) = gsl_matrix_get(eigenvec,i,j);
        }
    }

    gsl_matrix_free(tmpSigma);
    gsl_matrix_free(eigenvec);
    gsl_vector_free(eigenval); 
    gsl_eigen_symmv_free(w);
    
    return vals;
}
