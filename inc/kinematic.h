#ifndef KINEMATIC_H
#define KINEMATIC_H

//---LIBRARIES---//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>

//---HEADERS---//

#include "globalconstants.h"
#include "mathutils.h"

// Error codes

#define ERR_NO_CONVERGE -4
#define ERR_SINGULAR -5

// FORW
int forward_kinematics(double *solution_matrix, const double *dh_array, size_t n_links, const double *joint_angles_deg);
int fill_transformation_matrix(double *trans_matrix, const size_t index, const double *dh_array, const double angle_coordinate, const size_t n_links);
int fill_dh_parameters(double *dh_array, size_t n_links, const double *links_lengths);

// INVR
int inverse_kinematics_dls(double *q_out, const double *q_init,
                          const double T_des[16], const double *dh_array,
                          size_t n_links, double tol_pos, double tol_ori,
                          int max_iters, double lambda0, double max_delta_q,
                          const double *joint_min, const double *joint_max,
                          int *iters_out);
int build_jacobian(double *J, const double *q_rad, const double *dh_array,
                  size_t n_links, const double *T_current);
int compute_orientation_error(const double *R_des, const double *R_cur, double *error);

// SERV
int mat4_mul(double *left_matrix, const double *right_matrix);
int mat6_mul(double *A, const double *B, double *result);
int mat6_transpose(const double *A, double *result);
int mat6_add(double *A, const double *B, double *result);
int mat6_scale(double *A, double scalar, double *result);
int mat6_identity(double *I);

int solve_dls_system(const double *JTJ_plus_lambda2I, const double *JT_error, 
                     double *delta_q, int n, double lambda);
int mat6_mul_vector(const double *A, const double *v, double *result, int rows, int cols);
int rotation_matrix_to_angle_axis(const double *R, double *angle, double *axis);
int mat3_transpose(const double *A, double *result);
int mat3_mul(const double *A, const double *B, double *result);

#endif
