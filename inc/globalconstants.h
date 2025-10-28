#ifndef GLOBALCONSTANTS_H
#define GLOBALCONSTANTS_H

    // error cods //

#define ERR_OK 0
#define ERR_NULL -1
#define ERR_BADARGS -2
#define ERR_INNNER -3
    
    // constants // 

#define CIRCLE 360
#define DEFAULT_N_LINKS 6
#define DH_COLS 4
#define TRANSFORM_MATRIX_DIM 4
#define TRANSFORM_MATRIX_SIZE (TRANSFORM_MATRIX_DIM * TRANSFORM_MATRIX_DIM)
#define JACOBIAN_SIZE 36

    // matrix indices for 4x4 transformation matrix //
#define TX 3    // translation X index
#define TY 7    // translation Y index  
#define TZ 11   // translation Z index
#define XX 0    // X-axis X component
#define XY 4    // X-axis Y component
#define XZ 8    // X-axis Z component
#define YX 1    // Y-axis X component
#define YY 5    // Y-axis Y component
#define YZ 9    // Y-axis Z component
#define ZX 2    // Z-axis X component
#define ZY 6    // Z-axis Y component
#define ZZ 10   // Z-axis Z component
  
    //  default lenghts  //

#define DEFAULT_L1 1.2
#define DEFAULT_L2 0.6
#define DEFAULT_L3 0.55
#define DEFAULT_L4 0
#define DEFAULT_L5 0
#define DEFAULT_L6 0.2

#endif
