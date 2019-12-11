#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"
// Define output file name
#define OUTPUT_FILE "stencil.pgm"
#define MASTER 0 //MPI
int calc_ncols_from_rank(int rank, int size,int NCOLS);//MPI
void stencil(const int nx, const int ny, float * restrict image, float * restrict tmp_image);
void init_image(const int nx, const int ny, const int width, const int height,
                float* image, float* tmp_image);
void output_image(const char* file_name, const int nx, const int ny,
                  const int width, const int height, float* image);
double wtime(void);
int main(int argc, char* argv[])
{ 
  /**********MPI Parameters****************************************/
  int rank;              /* the rank of this process */
  int left;              /* the rank of the process to the left */
  int right;             /* the rank of the process to the right */
  int size;              /* number of processes in the communicator */
  int tag = 0;           /* scope for adding extra information to a message */
  MPI_Status status;     /* struct used by MPI_Recv */
  int local_nrows;       /* number of rows apportioned to this rank */
  int local_ncols;       /* number of columns apportioned to this rank */
  int local_ncole_cal;
  float *local_image;    /* local image grid */
  float *tmp_local_image;  /*tmp local image grid*/
  /****************************************************************/
  //MPI_Init returns once it has started up processor
  //Get size of cohort and rank for this process
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  // printf(" ***The total size is  %d !!! ***** !\n",size);
  // Check usage
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  // Initiliase problem dimensions from command line arguments
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  // we pad the outer edge of the image to avoid out of range address issues in
  // stencil
  int width = nx + 2;
  int height = ny + 2;
  /** determine process ranks to the left and right of rank
  ** respecting periodic boundary conditions*/
  left = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
  right = (rank + 1) % size;

  /** determine local grid size
  ** each rank gets all the rows, but a subset of the number of columns*/
  local_nrows = height;
  local_ncols = calc_ncols_from_rank(rank, size, width);
  // printf(" The rank  %d 's block number: %d!\n",rank,local_ncols);
  if (local_ncols < 1) {
    fprintf(stderr,"Error: too many processes:- local_ncols < 1\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // Allocate the image
  float* image = malloc(sizeof(float) * width * height);
  float* tmp_image = malloc(sizeof(float) * width * height);
  float* final_image = malloc(sizeof(float) * width * height);
  // Set the input image
  init_image(nx, ny, width, height, image, tmp_image);
  // If only one processor , perform as serial code
  //
  if (size == 1){
    double tic = wtime();
    for (int t = 0; t < niters; ++t) {
      stencil(width, height, image, tmp_image);
      stencil(width, height, tmp_image, image);
    }
    double toc = wtime();
    // Output
    printf("------------------------------------\n");
    printf(" runtime: %lf s\n", toc - tic);
    printf("------------------------------------\n");

    output_image(OUTPUT_FILE, nx, ny, width, height, image);

    free(image);
    free(tmp_image);
    free(final_image);

    return EXIT_SUCCESS;
  }
  // If hava more than one processor, perform by using mpi
  // Allocate the local image and local temp image
  /*Need to add 2 for halo exchanges*/
  if (rank == MASTER)  //Add only one for rank 0 and rank size-1
  {
  local_image = (float*) malloc (sizeof(float*) * local_nrows * (local_ncols + 1));
  tmp_local_image = (float*) malloc (sizeof(float*) * local_nrows * (local_ncols + 1));
  local_ncole_cal = local_ncols + 1;
  #ifdef GATHER
    for(int i = 0;i<(local_ncols+1)*local_nrows;i++)
      local_image[i] = 0.0f;
  #else
  #endif
  // printf(" The rank split %d successful!\n",rank);
  }else if (rank == (size -1))
  {
  local_image = (float*) malloc (sizeof(float*) * local_nrows * (local_ncols + 1));
  tmp_local_image = (float*) malloc (sizeof(float*) * local_nrows * (local_ncols + 1));
  local_ncole_cal = local_ncols + 1;
  #ifdef GATHER
    for(int i = 0;i<(local_ncols+1)*local_nrows;i++)
      local_image[i] = 0.0f;
  #else
  #endif
  // printf(" The rank split %d successful!\n",rank);
  }else
  {
  local_image = (float*) malloc (sizeof(float*) * local_nrows * (local_ncols + 2));
  tmp_local_image = (float*) malloc (sizeof(float*) * local_nrows * (local_ncols + 2));
  local_ncole_cal = local_ncols + 2;
  #ifdef GATHER
    for(int i = 0;i<(local_ncols+2)*local_nrows;i++)
      local_image[i] = 0.0f;
  #else
  #endif
  // printf(" The rank split %d successful!\n",rank);
  }
  /****************************** initialize **********************************/
  #ifdef GATHER
    double ini_tic = wtime();
    int *displs; //Use to calculate the offset of the scatter and gather
    int *sendcounts; //Use to calculate the number of the scatter and gather
    displs = (int *) malloc (sizeof (int) * size);
    sendcounts = (int *) malloc (sizeof (int) * size);
    for (int i = 0; i <size; i++){
    int cal_tmp = calc_ncols_from_rank(i, size, width);
    sendcounts[i] = cal_tmp * height; 
    displs[i] = i * local_ncols * height;
    if (i == 0){
      sendcounts[i] -= 1;
      displs[i] += height;
    }
    if (i == size -1)
      sendcounts[i] -= 1;
    }
    MPI_Scatterv(image,sendcounts,displs,MPI_FLOAT,&local_image[height],local_ncols * height,MPI_FLOAT,
                  MASTER, MPI_COMM_WORLD);
    double ini_toc = wtime();
    // printf(" Using Gather function!\n");
  #else
    // printf(" Using manual gather function!\n");
  //initialize the local some boundary for grid
  double ini_tic = wtime();
  if (rank == MASTER)
  {
  // the first block just add one new halo
    for (int i = 0; i < local_ncols ; ++i) { // total number of the the cols is the cale(rank,size)
      for (int j = 0; j < local_nrows; ++j){ //total number of the rows is height
        local_image [j + i * local_nrows] = image [j + ((rank*(width/size))+i)*local_nrows];
        tmp_local_image [j + i * local_nrows] = image [j + ((rank*(width/size))+i)*local_nrows];
        
      }
    }
  //  printf(" ***The rank %d initialize success ***** !\n",rank);
  }else if (rank == (size -1) )
  {
  // the last block just add one new halo
    // the first block just one halo
    for (int i = 1; i < local_ncols + 1; ++i) { // total number of the the cols is the cale(rank,size)
      for (int j = 0; j < local_nrows; ++j) {  //total number of the rows is height
        local_image [j + i * height] = image [j + ((rank*(width/size)-1)+i)*height];
        tmp_local_image [j + i * height] = image [j + ((rank*(width/size)-1)+i)*height];
      }
    
  }
  // printf(" ***The rank %d initialize success ***** !\n",rank);
  }else 
  {
  // inter blocks
  for (int i = 1; i < local_ncols + 1; ++i) { // total number of the the cols is the cale(rank,size)
    for (int j = 0; j < local_nrows; ++j) { //total number of the rows is height
      local_image [j + i * height] = image [j + ((rank*(width/size)-1)+i)*height];
      tmp_local_image [j + i * height] = image [j + ((rank*(width/size)-1)+i)*height];
    }

  }
  // printf(" ***The rank %d initialize success ***** !\n",rank);
  }
  double ini_toc = wtime();
  #endif
  //Evaluate the different types of initialize
    if (rank == MASTER)
  {
    // Output
    printf("------------------------------------\n");
    printf(" runtime of the initilize: %lf s\n", ini_toc-ini_tic);
    printf("------------------------------------\n");
    printf("With %d processes\n", size);
  }
  /****************************** Asm the halo before the stencil **********************************/
   ///SEND and Receive
   //send left, receive right
  // printf(" The initial Halo begin");
  if (rank != 0)
  {
    // printf(" ***The rank %d send left successful ***** !\n",rank);
    MPI_Ssend(&local_image[1*height], height, MPI_FLOAT, left, 0, MPI_COMM_WORLD);
  }
  if (rank != size-1)
  { 
    if (rank == MASTER){
       
       MPI_Recv(&local_image[(local_ncols)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &status);
      //  printf(" ***The rank %d receive right successful ***** !\n",rank);
    }else{
      
       MPI_Recv(&local_image[(local_ncols+1)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &status);
      //  printf(" ***The rank %d receive right successful ***** !\n",rank); 
    }

  }
  // send right ,receive left
  if (rank != size-1)
  {
    if (rank == 0){
    MPI_Ssend(&local_image[(local_ncols-1)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD);
    // printf(" ***The rank %d send right successful ***** !\n",rank); 
    }else{
    MPI_Ssend(&local_image[(local_ncols)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD);
    // printf(" ***The rank %d send right successful ***** !\n",rank); 
    }
    
  }
  if (rank != 0)
  {
    MPI_Recv(&local_image[0], height, MPI_FLOAT, left, 0, MPI_COMM_WORLD, &status);
    // printf(" ***The rank %d receive left successful ***** !\n",rank); 
  }
  // printf(" The initial Halo successful \n!");


  /*********************************************************************
   * 
   * 
   */
  /////////////////////////Iterate begin////////////////////////////////////////////
  double tic = wtime();
  for (int t = 0; t < niters; ++t)
  { 
    // printf(" *** The rank %d 's %d interate begin ***** !\n",rank,t); 

    //////////////////////stencil first//////////////////
    stencil(local_ncole_cal, height, local_image, tmp_local_image);
    // printf(" *** The rank %d 's %d interate first calucate successful! ***** !\n",rank,t);
    /////////////////////exchange halos//////////////////////
    //send left and receive right
    if (rank != MASTER)
    {
      MPI_Send(&tmp_local_image[1*height], height, MPI_FLOAT, left, 0, MPI_COMM_WORLD);
    }
    if (rank != size-1)
    {
     if (rank == MASTER){
       MPI_Recv(&tmp_local_image[(local_ncols)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &status);
     }else{
       MPI_Recv(&tmp_local_image[(local_ncols+1)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &status);
     }      
       
     // MPI_Recv(&tmp_local_image[(local_ncols+1)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &status);
    }
    //send right and receive left
    if (rank != size-1)
    {
      if (rank == 0){
        MPI_Ssend(&tmp_local_image[(local_ncols-1)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD);
        }else{
        MPI_Ssend(&tmp_local_image[(local_ncols)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD);
      }
    //  MPI_Send(&tmp_local_image[(local_ncols)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD);
    }
    if (rank != 0)
    {
      MPI_Recv(&tmp_local_image[0], height, MPI_FLOAT, left, 0, MPI_COMM_WORLD, &status);
    }


    
    ///////////////////////stencil second/////////////////////////// 
      stencil(local_ncole_cal, height, tmp_local_image, local_image);
      // printf(" *** The rank %d 's %d interate second calucate successful! ***** !\n",rank,t);
      ////send left and receive right
      if (rank != MASTER)
      {
        // printf(" ***The rank %d send left successful ***** !\n",rank);
        MPI_Ssend(&local_image[1*height], height, MPI_FLOAT, left, 0, MPI_COMM_WORLD);
      }
      if (rank != size-1)
      { 
        if (rank == MASTER){
          
          MPI_Recv(&local_image[(local_ncols)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &status);
          // printf(" ***The rank %d receive right successful ***** !\n",rank);
        }else{
          
          MPI_Recv(&local_image[(local_ncols+1)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD, &status);
          // printf(" ***The rank %d receive right successful ***** !\n",rank); 
        }
      }
      //send right and receive left
      if (rank != size-1)
      {
        if (rank == MASTER){
        MPI_Ssend(&local_image[(local_ncols-1)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD);
        // printf(" ***The rank %d send right successful ***** !\n",rank); 
        }else{
        MPI_Ssend(&local_image[(local_ncols)*height], height, MPI_FLOAT, right, 0, MPI_COMM_WORLD);
        // printf(" ***The rank %d send right successful ***** !\n",rank); 
        }
      }
      if (rank != 0)
      {
        MPI_Recv(&local_image[0], height, MPI_FLOAT, left, 0, MPI_COMM_WORLD, &status);
        // printf(" ***The rank %d receive left successful ***** !\n",rank); 
      }
    // printf(" *** The rank %d 's %d interate successful ***** !\n",rank,t);

    
   
  //  if (rank == 0){
  //   for(int i = 0; i<130;i++)
  //     printf("The %d block has number: %f \n",i,tmp_local_image[135]);
  //   }
  }
  double toc = wtime();
  ///
  //Collect result
  //tmp_local_image[135] = 77;
  #ifdef GATHER
  int tmp_local_offset;
  tmp_local_offset = (local_ncols) * height;
  if (rank == MASTER)
      tmp_local_offset = (local_ncols - 1) * height;
  if (rank == size - 1)
      tmp_local_offset = (local_ncols - 1) * height;
  MPI_Gatherv(&local_image[height], tmp_local_offset, MPI_FLOAT,
                final_image, sendcounts, displs,
                MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  #else

  if (rank == MASTER) 
    { // need collect the left halo
    for(int j=0; j<local_ncols; ++j){
      for(int i=0;i<local_nrows; ++i) 
      {
        final_image[i+j*height] = local_image[i+j*height];
      }
    }

    for (int rank_number = 1; rank_number<size; ++rank_number)
    {
      int cols = calc_ncols_from_rank(rank_number, size, height);
      //if (rank_number != size - 1)
        for(int j=1; j<cols+1; ++j)
        {
          MPI_Recv(&final_image[((rank_number*(width/size)-1)+j)*height], height, MPI_FLOAT, rank_number, 0, MPI_COMM_WORLD, &status);
        }
    }
  }
  else
  { 
    for(int i=1; i<local_ncols+1; ++i)
      { //send back to the master processor

        MPI_Ssend(&local_image[i*height], height, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
      }
  }

  #endif

  if (rank == MASTER)
  {
    // Output
    printf("------------------------------------\n");
    printf(" runtime: %lf s\n", toc-tic);
    printf("------------------------------------\n");
    printf("With %d processes\n", size);

    output_image(OUTPUT_FILE, nx, ny, width, height,final_image);
    
  }

  /* don't forget to tidy up when we're done */
  MPI_Finalize();

  /* free up allocated memory */
  free(final_image);
  free(image);
  free(tmp_image);
  free(local_image);
  free(tmp_local_image);


  /* and exit the program */
  return EXIT_SUCCESS;
  // Call the stencil kernel
 // double tic = wtime();
 // for (int t = 0; t < niters; ++t) {
 //   stencil(nx, ny, width, height, image, tmp_image);
 //   stencil(nx, ny, width, height, tmp_image, image);
 // }
 // double toc = wtime();

  // Output
 // printf("------------------------------------\n");
 // printf(" runtime: %lf s\n", toc - tic);
 // printf("------------------------------------\n");

 // output_image(OUTPUT_FILE, nx, ny, width, height, image);
 // free(image);
 // free(tmp_image);
}

void stencil(const int nx, const int ny, float *  image, float *  tmp_image)
{ // local_ncols, height, local_image, tmp_local_image
  // first row
  //need condider the different ranks 0
  //rank0: 0 1 2 3 4
  //         1 2 3
  //rnak inter: 0 1 
  float tp_col, tp_row;
  for (int i = 1; i < nx-1 ; ++i) {
    for (int j = 1; j < ny-1 ; ++j) {
       tp_row = image[j     + i       * ny] * 3.0f/5.0f + image[j - 1 + i       * ny] * 0.5f/5.0f + image[j + 1 + i       * ny] * 0.5f/5.0f;
       tp_col = image[j     + (i - 1) * ny] * 0.5f/5.0f + image[j     + (i + 1) * ny] * 0.5f/5.0f;
      // tmp_image[j + i * ny] =  image[j     + i       * ny] * 3.0f / 5.0f;
      // tmp_image[j + i * ny] += image[j     + (i - 1) * ny] * 0.5f / 5.0f;
      // tmp_image[j + i * ny] += image[j     + (i + 1) * ny] * 0.5f / 5.0f;
      // tmp_image[j + i * ny] += image[j - 1 + i       * ny] * 0.5f / 5.0f;
      // tmp_image[j + i * ny] += image[j + 1 + i       * ny] * 0.5f / 5.0f;
     tmp_image[j + i * ny] = tp_row + tp_col;
    }
  }
}

/*
{ //2.4s
  // O3:1.19
  float tp_col, tp_row;
  for (int i = 1; i < ny + 1; ++i) {
    for (int j = 1; j < nx + 1; ++j) {
       tp_row = image[j     + i       * ny] * 3.0f / 5.0f + image[j - 1 + i       * ny] * 0.5f / 5.0f + image[j + 1 + i       * ny] * 0.5f / 5.0f;
       tp_col = image[j     + (i - 1) * ny] * 0.5f / 5.0f + image[j     + (i + 1) * ny] * 0.5f / 5.0f;
       tmp_image[j + i * ny] = tp_row + tp_col;
    }
  }
}
*/
/*             
{
  for (int j = 1; j < ny + 1; ++j) {
    for (int i = 1; i < nx + 1; ++i) {
      tmp_image[j + i * height] =  image[j     + i       * height] * 3.0 / 5.0;
      tmp_image[j + i * height] += image[j     + (i - 1) * height] * 0.5 / 5.0;
      tmp_image[j + i * height] += image[j     + (i + 1) * height] * 0.5 / 5.0;
      tmp_image[j + i * height] += image[j - 1 + i       * height] * 0.5 / 5.0;
      tmp_image[j + i * height] += image[j + 1 + i       * height] * 0.5 / 5.0;
    }
  }
}
*/
// Create the input image
void init_image(const int nx, const int ny, const int width, const int height,
                float* image, float* tmp_image)
{
  // Zero everything
  for (int j = 0; j < ny + 2; ++j) {
    for (int i = 0; i < nx + 2; ++i) {
      image[j + i * height] = 0.0f;
      tmp_image[j + i * height] = 0.0f;
    }
  }

  const int tile_size = 64;
  // checkerboard pattern
  for (int jb = 0; jb < ny; jb += tile_size) {
    for (int ib = 0; ib < nx; ib += tile_size) {
      if ((ib + jb) % (tile_size * 2)) {
        const int jlim = (jb + tile_size > ny) ? ny : jb + tile_size;
        const int ilim = (ib + tile_size > nx) ? nx : ib + tile_size;
        for (int j = jb + 1; j < jlim + 1; ++j) {
          for (int i = ib + 1; i < ilim + 1; ++i) {
            image[j + i * height] = 100.0f;
          }
        }
      }
    }
  }
}

// Routine to output the image in Netpbm grayscale binary image format
void output_image(const char* file_name, const int nx, const int ny,
                  const int width, const int height, float* image)
{
  // Open output file
  FILE* fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  double maximum = 0.0;
  for (int j = 1; j < ny + 1; ++j) {
    for (int i = 1; i < nx + 1; ++i) {
      if (image[j + i * height] > maximum) maximum = image[j + i * height];
    }
  }
   
  // Output image, converting to numbers 0-255
  for (int j = 1; j < ny + 1; ++j) {
    for (int i = 1; i < nx + 1; ++i) {
      fputc((char)(255.0f * image[j + i * height] / maximum), fp);
      // fputc((char)(image[j + i * height]), fp);
    }
  }

  // Close the file
  fclose(fp);
}

// Get the current time in seconds since the Epoch
double wtime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}
//MPI_Each——
int calc_ncols_from_rank(int rank, int size, int NCOLS)
{
  int ncols;

  ncols = NCOLS / size;       /* integer division */
  if ((NCOLS % size) != 0) {  /* if there is a remainder */
    if (rank == size - 1)
      ncols += NCOLS % size;  /* add remainder to last rank */
  }
  
  return ncols;
}