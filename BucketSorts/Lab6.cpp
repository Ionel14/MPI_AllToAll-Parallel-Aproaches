//#include <mpi.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <malloc.h>
//#include <time.h>
//
//int MPI_Bucket_sort(int n, double* v, double max, int root, double& communicationTime, MPI_Comm comm);
//
//double* merge_array(int n, double* a, int m, double* b);
//void     merge_sort(int n, double* a);
//void     swap(double* a, double* b);
//
//int main(int argc, char** argv)
//{
//    int rank, size, n = 10000000;
//    double* v, max = 30;
//
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);
//
//    v = (double*)calloc(n, sizeof(double));
//
//    if (rank == 0)
//    {
//        for (int i = 0; i < n; i++)
//            v[i] = (rand() % (1000 * (int)max)) / 1000.;
//        //for (int i = 0; i < n; i++)
//        //        printf("%lf ", v[i]);
//        puts("");
//    }
//
//
//    double t = MPI_Wtime(), allTime, communicationTime = 0, auxTime;
//
//    MPI_Bucket_sort(n, v, max, 0, communicationTime, MPI_COMM_WORLD);
//
//    t = MPI_Wtime() - t;
//    MPI_Reduce(&t, &allTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&communicationTime, &auxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
//
//    if (rank == 0)
//    {
//        printf("AllTime: %lf\n", allTime);
//        printf("ComputationTime: %lf\n", allTime - auxTime);
//        printf("CommunicationTime: %lf\n", auxTime);
//        //for (int i = 0; i < n; i++)
//        //        printf("%lf ", v[i]);
//    }
//    free(v);
//    MPI_Finalize();
//    return 0;
//}
//
//int MPI_Bucket_sort(int n, double* v, double max, int root, double &communicationTime, MPI_Comm comm)
//{
//    //declarare varibile
//    int rank, size, local_n = 0, * counts, * displs;
//    double* local_v, time;
//
//    time = MPI_Wtime();
//    //determinare rank si size
//    MPI_Comm_rank(comm, &rank);
//    MPI_Comm_size(comm, &size);
//    time = MPI_Wtime() - time;
//    communicationTime += time;
//
//    //alocare memorie
//    local_v = (double*)calloc(n, sizeof(double));
//    counts = (int*)calloc(size, sizeof(int));
//    displs = (int*)calloc(size, sizeof(int));
//    
//    //Bcast
//    time = MPI_Wtime();
//    MPI_Bcast(v, n, MPI_DOUBLE, root, comm);
//    time = MPI_Wtime() - time;
//    communicationTime += time;
//
//    //generare vector local
//    for (int i = 0; i < n; i++)
//        if (v[i] >= rank * max / size && v[i] < (rank + 1) * max / size)
//            local_v[local_n++] = v[i];
//
//    //sortare vector local
//    merge_sort(local_n, local_v);
//
//    //determinare counts si displs 
//    time = MPI_Wtime();
//    MPI_Gather(&local_n, 1, MPI_INT, counts, 1, MPI_INT, root, comm);
//    time = MPI_Wtime() - time;
//    communicationTime += time;
//
//    for (int i = 1; i < size; i++)
//        displs[i] = displs[i - 1] + counts[i - 1];
//    
//    //gather pe root
//    time = MPI_Wtime();
//    MPI_Gatherv(local_v, local_n, MPI_DOUBLE, v, counts, displs, MPI_DOUBLE, root, comm);
//    time = MPI_Wtime() - time;
//    communicationTime += time;
//    
//    free(displs);
//    free(counts);
//    free(local_v);
//    return MPI_SUCCESS;
//}
//
//double* merge_array(int n, double* a, int m, double* b) {
//
//    int i, j, k;
//    double* c = (double*)calloc(n + m, sizeof(double));
//
//    for (i = j = k = 0; (i < n) && (j < m);)
//
//        if (a[i] <= b[j])c[k++] = a[i++];
//        else c[k++] = b[j++];
//
//    if (i == n)for (; j < m;)c[k++] = b[j++];
//    else for (; i < n;)c[k++] = a[i++];
//
//    return c;
//}
//
//// function to merge sort the array a with n elements
//
//void merge_sort(int n, double* a) {
//
//    double* c;
//    int i;
//
//    if (n <= 1) return;
//
//    if (n == 2) {
//
//        if (a[0] > a[1])swap(&a[0], &a[1]);
//        return;
//    }
//
//
//
//    merge_sort(n / 2, a); merge_sort(n - n / 2, a + n / 2);
//
//    c = merge_array(n / 2, a, n - n / 2, a + n / 2);
//
//    for (i = 0; i < n; i++)a[i] = c[i];
//
//    return;
//}
//
//
//// swap two doubles
//void swap(double* a, double* b) {
//
//    double temp;
//
//    temp = *a; *a = *b; *b = temp;
//
//}
