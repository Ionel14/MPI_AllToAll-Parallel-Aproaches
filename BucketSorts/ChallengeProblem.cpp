#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>

int MPI_Bucket_sort(int n, double* v, double max, int root, double &communicationTime, MPI_Comm comm);

void GetBuckets(int local_n, int size, double max, double** buckets_array, double* local_array, int* counts);

void matrixToArray(double** matrix, int rows, int* cols, double* array);
double* merge_array(int n, double* a, int m, double* b);
void merge_sort(int n, double* a);
void swap(double* a, double* b);

int main(int argc, char** argv)
{
    int rank, size, n = 5000000;
    double* v, max = 30;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    v = (double*)calloc(n, sizeof(double));

    if (rank == 0)
    {
        for (int i = 0; i < n; i++)
            v[i] = (rand() % (1000 * (int)max)) / 1000.;
        /*for (int i = 0; i < n; i++)
                printf("%lf ", v[i]);
        puts("");*/
    }


    double t = MPI_Wtime(), allTime, auxTime, communicationTime = 0;

    MPI_Bucket_sort(n, v, max, 0, communicationTime, MPI_COMM_WORLD);

    t = MPI_Wtime() - t;
    MPI_Reduce(&t, &allTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&communicationTime, &auxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("AllTime: %lf\n", allTime);
        printf("ComputationTime: %lf\n", allTime - auxTime);
        printf("CommunicationTime: %lf\n", auxTime);
        /*for (int i = 0; i < n; i++)
                printf("%lf ", v[i]);*/
    }
    free(v);
    MPI_Finalize();
    return 0;
}

int MPI_Bucket_sort(int n, double* v, double max, int root, double &communicationTime, MPI_Comm comm)
{
    //declare vars
    int rank, size, local_n, *counts, *newCounts, *displs, *displs_send;
    double *local_array, *big_bucket, **buckets_array, time;
    int bigCount = 0;
    communicationTime = MPI_Wtime();
    //get rank and size
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    communicationTime = MPI_Wtime() - communicationTime;

    local_n = n / size;

    //allocate local_array
    buckets_array = (double**)calloc(size, sizeof(double*));
    big_bucket = (double*)calloc(n, sizeof(double));
    local_array = (double*)calloc(local_n, sizeof(double));
    counts = (int*)calloc(size, sizeof(int));
    newCounts = (int*)calloc(size, sizeof(int));
    displs = (int*)calloc(size, sizeof(int));
    displs_send = (int*)calloc(size, sizeof(int));
    
    //scatter array to local_array
    time = MPI_Wtime();
    int error = MPI_Scatter(v, local_n, MPI_DOUBLE, local_array, local_n, MPI_DOUBLE, root, comm);
    time = MPI_Wtime() - time;
    communicationTime += time;

    if (error != MPI_SUCCESS)
    {
        free(local_array);
        return error;
    }

    //Separate local array on buckets
    GetBuckets(local_n, size, max, buckets_array, local_array, counts);

    //Get counts of each processor about the buckets
    time = MPI_Wtime();
    MPI_Alltoall(counts, 1, MPI_INT, newCounts, 1, MPI_INT, comm);
    time = MPI_Wtime() - time;
    communicationTime += time;

    //Calculate the nr of elements that I have to get, local displacements and the displacements of the bucket I get
    bigCount = newCounts[0];
    for (int i = 1; i < size; i++)
    {
        displs_send[i] = displs_send[i - 1] + counts[i - 1];
        displs[i] = displs[i - 1] + newCounts[i - 1];
        bigCount += newCounts[i];
    }

    matrixToArray(buckets_array, size, counts, local_array);

    //Build the big bucket vith values gathered from each processor
    time = MPI_Wtime();
    MPI_Alltoallv(local_array, counts, displs_send, MPI_DOUBLE, big_bucket, newCounts, displs, MPI_DOUBLE, comm);
    time = MPI_Wtime() - time;
    communicationTime += time;

    //In case you want to see each big_bucket uncomment this
    //printf("\n\nrank: %d\n", rank);
    //for (int i = 0; i < bigCount; i++)
    //{
    //    printf("buckets[%d] = %lf\n", i, big_bucket[i]);

    //}

    merge_sort(bigCount, big_bucket);

    //Get counts and displs 
    time = MPI_Wtime();
    MPI_Gather(&bigCount, 1, MPI_INT, counts, 1, MPI_INT, root, comm);
    time = MPI_Wtime() - time;
    communicationTime += time;

    displs[0] = 0;
    for (int i = 1; i < size; i++)
        displs[i] = displs[i - 1] + counts[i - 1];

    //gather on root
    time = MPI_Wtime();
    MPI_Gatherv(big_bucket, bigCount, MPI_DOUBLE, v, counts, displs, MPI_DOUBLE, root, comm);
    time = MPI_Wtime() - time;
    communicationTime += time;

    free(local_array);
    free(counts);
    free(newCounts);
    free(displs);
    free(displs_send);
    free(big_bucket);
    return MPI_SUCCESS;
}

//function to get the array from a matrix and free the memory of the matrix
void matrixToArray(double** matrix, int rows, int* cols, double* array) {
    int index = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols[i]; j++) {
            array[index++] = matrix[i][j];
        }
        free(matrix[i]);
    }
    free(matrix);
}

// function to separate the local array on buckets
void GetBuckets(int local_n, int size, double max, double** buckets_array, double* local_array, int *counts)
{
    for (int i = 0; i < size; i++)
    {
        buckets_array[i] = (double*)calloc(local_n, sizeof(double));
    }

    double maxBucket0 = max / size;
    for (int i = 0; i < local_n; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (local_array[i] >= maxBucket0 * j && local_array[i] < maxBucket0 * (j + 1))
            {
                buckets_array[j][counts[j]++] = local_array[i];
                break;
            }
        }
    }
}

double* merge_array(int n, double* a, int m, double* b) {

    int i, j, k;
    double* c = (double*)calloc(n + m, sizeof(double));

    for (i = j = k = 0; (i < n) && (j < m);)

        if (a[i] <= b[j])c[k++] = a[i++];
        else c[k++] = b[j++];

    if (i == n)for (; j < m;)c[k++] = b[j++];
    else for (; i < n;)c[k++] = a[i++];

    return c;
}

// function to merge sort the array a with n elements
void merge_sort(int n, double* a) {

    double* c;
    int i;

    if (n <= 1) return;

    if (n == 2) {

        if (a[0] > a[1])swap(&a[0], &a[1]);
        return;
    }

    merge_sort(n / 2, a); merge_sort(n - n / 2, a + n / 2);

    c = merge_array(n / 2, a, n - n / 2, a + n / 2);

    for (i = 0; i < n; i++)a[i] = c[i];

    return;
}

// swap two doubles
void swap(double* a, double* b) {

    double temp;

    temp = *a; *a = *b; *b = temp;

}