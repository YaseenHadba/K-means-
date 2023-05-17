#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <Python.h>

void run(int k , int max_iter, double epsilon,int d,int num_vectors, double **vectors, double **init_centroids);
double static **centroids; /*defined as a global variable in order to be able to free the space it occupied */

static PyObject *fit(int k , int max_iter, double epsilon,int d, PyObject *vectors, PyObject *init_centroids){
    double **vMatrix, **starting;
    PyObject  *py_vector, *py_centroid, *py_result, *tmp;
    Py_ssize_t i, j;
    int N;
    N = (int) PyObject_Length(vectors);
    vMatrix = (double**)malloc( N*sizeof(double*));
    if(vMatrix == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for (i=0;i<N;i++){
        vMatrix[i] = (double*) malloc(d*sizeof(double ));
        if(vMatrix[i] == NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    starting = (double**) malloc(k*sizeof (double*));
    if(starting == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for (i=0;i<k;i++){
        starting[i]=(double*)malloc(d*sizeof (double));
        if(starting[i] == NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0;i<N;i++){
        py_vector = PyList_GetItem(vectors, i);
        for (j=0;j<d;j++){
           vMatrix[i][j] = PyFloat_AsDouble(PyList_GetItem(py_vector, j));
        }
    }
    for (i=0;i<k;i++){
        py_centroid = PyList_GetItem(init_centroids, i);
        for (j=0;j<d;j++){
            starting[i][j] = PyFloat_AsDouble(PyList_GetItem(py_centroid, j));
        }
    }
    run(k, max_iter, epsilon, d, N, vMatrix, starting);
    py_result = PyList_New(k);
    for(i=0;i<k;i++)
    {
        tmp = PyList_New(d);
        for(j=0;j<d;j++)
        {
            PyList_SetItem(tmp,j,PyFloat_FromDouble(centroids[i][j]));
        }
        PyList_SetItem(py_result,i,tmp);
    }
    for(i=0;i<k;i++){free(centroids[i]);}
    free(centroids);
    return py_result;
}


static PyObject *fit_api(PyObject *self,PyObject *args){
    int k, max_iter,d;
    double epsilon;
    PyObject *vectors;
    PyObject *init_centroids;
    if(!PyArg_ParseTuple(args,"iidiOO",&k,&max_iter,&epsilon,&d,&vectors,&init_centroids)){return NULL;}
    return Py_BuildValue("O", fit(k,max_iter,epsilon,d,vectors,init_centroids));
}


static PyMethodDef fit_apiMethods[]={
        {"fit",(PyCFunction)fit_api, METH_VARARGS,
                        PyDoc_STR("kmeans++ hw2")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        fit_apiMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}

double* divide(double* vector,int num,int vDim){
    double *result;
    int i;
    result = (double*) malloc(vDim*sizeof(double ));
    if(result == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for (i = 0; i < vDim; i++) {
        result[i]=  vector[i]/num;
    }
    return result;
}

double help(double* vector1, double* vector2, int vDim){ /* help func for calculating the norm*/ 
    double* arr;
    int i;
    double sum;
    arr = (double*) malloc(vDim*sizeof (double));
    if(arr ==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0;i<vDim;i++){
        arr[i]=vector1[i]-vector2[i];
    }
    sum=0;
    for(i=0;i<vDim;i++){
        sum = sum + arr[i]*arr[i];
    }
    free(arr);
    return sum;
}

int calc_norm(double** centroids, double* point,int k, int vDim){
    double closest;
    double x;
    int i;
    int closest_norm;
    closest_norm =0;
    closest = help(centroids[0],point,vDim);
    for(i=1;i<k;i++){
        x = help(centroids[i],point,vDim);
        if(x < closest){
            closest_norm=i;
            closest =x;
        }
    }
    return closest_norm;
}

void run(int k , int max_iter, double epsilon,int d,int num_vectors, double **vectors, double **init_centroids) {
    int vDim;
    int i;
    int j;
    int index;
    int q;
    int flag;
    int N;
    double **matrix;
    double **copy;
    int *clusters_size;
    flag=-1;
    N=num_vectors;
    vDim = d;
    q=0;
    centroids = (double **) malloc(k*sizeof(double*));
    for (i=0;i<k;i++){
        centroids[i] = (double*) malloc(d*sizeof (double));
    }
    for(i=0;i<k;i++){
        for (j=0;j<d;j++) {
            centroids[i][j] = init_centroids[i][j];
        }
    }
    for (i =0 ;i<k;i++){
        free(init_centroids[i]);
    }
    free(init_centroids);
    matrix = (double **) malloc(k * sizeof(double *));
    if (matrix == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    for (i = 0; i < k; i++) {
        matrix[i] = (double *) malloc(vDim * sizeof(double));
        if (matrix[i] == NULL) {
            printf("An Error Has Occurred");
            exit(1);
        }
    }  
    clusters_size = (int *) malloc(k * sizeof(int));
    if (clusters_size == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }  
    copy = (double **) malloc(k * sizeof(double *));
    if (copy == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    for (j = 0; j < k; j++) {
        copy[j] = (double *) malloc(vDim * sizeof(double));
    }
    while (q < max_iter && flag==-1) {    
        for (i = 0; i < k; i++) {
            clusters_size[i] = 0;
            for (j = 0; j < vDim; j++) {
                matrix[i][j] = 0;
            }
        }
        
        for (j = 0; j < N; j++) {
            index = calc_norm(centroids, vectors[j], k,vDim);
            clusters_size[index] += 1;
            for (i = 0; i < vDim; i++) {
                matrix[index][i] += vectors[j][i];
            }
        }
        for (i = 0; i < k; i++) {
            for (j = 0; j < vDim; j++) {
                copy[i][j] = centroids[i][j];
            }
        }      
        for (i = 0; i < k; i++) {
            free(centroids[i]);
            centroids[i] = divide(matrix[i], clusters_size[i], vDim);
        }     
        for (i = 0; i < k; i++) {
            if (help(centroids[i], copy[i], vDim) > (epsilon * epsilon)) {
                flag = -1;
                break;
            }
            flag = 1;
        }
        q += 1;
    }   
    for (i = 0; i < N; ++i) {free(vectors[i]);}
    free(vectors);
    for (i = 0; i < k; i++) {
        free(copy[i]);
        free(matrix[i]);
    }
    free(matrix);
    free(clusters_size);
    free(copy);
}

