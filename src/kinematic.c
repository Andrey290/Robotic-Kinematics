//---HEADERS---//

#include "../inc/kinematic.h"
#include "../inc/mathutils.h"

//---LIBRARIES---//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

int forward_kinematics(double *solution_matrix,
		       const double *dh_array, size_t n_links,
	       	       const double *joint_angles_deg) {
	/**
	* forward_kinematics
	*
	* Briefly:
	*   This function contain implement of forward kinematics task's solution.
	*   It requires manipulator and cross coordinates vector. That means 
	*   Denavit-Hartemberg matrix and vector of angles.
	*   In a circle it calls fill_transformation_matrix for right (i+1) matrix
	*   then two_matrix_multiplication(left_m, right_m). 
	*   For the first iteration left_m = E and right = A_1.
	*   +--------------------------------------+
	*   | Formula: T=E*A_1*A_2*A_3*A_4*A_5*A_6 |
	*   +--------------------------------------+ 
	*
	* Array format:
 	*   dh_array[i][0] = a_i       (link length)         [units: meters]
 	*   dh_array[i][1] = alpha_i   (link twist)          [units: degrees]
 	*   dh_array[i][2] = d_i       (link offset)         [units: meters]
 	*   dh_array[i][3] = theta_i   (joint angle)         [units: degrees]
	*
	* Preconditions:
	*   - dh_array != NULL
	*   - joint_angles_deg != NULL
	*   - solution_matrix != NULL
	*   - dh_array is filled by valid values
	*   - dh_array is a flat array with size n_links * DH_COLS
	*   - solution_matrix is a flat array with size 4*4 (16 doubles) filled as E-matrix
	*
	* Postconditions:
	*   - if returns 0 it successfully calculated solution_matrix
	*   - in case of error it returns not a zero code and array may be partially filled
	*
	* Responds:
	*   0 - OK
	*  -1 - NULL pointer
	*  -2 - invalid index or invalid memory layout (caller error) 
	*  -3 - Inner error
	*
	* Side effects:
	*  - memory modification on solution_matrix addreess 
	*  - calls:
	*     - fill_transformation_matrix_matrix
	*     - mat4_mul
	* Complexity:
	*  - Time O(n)
	*  - Memory O(1)
	*
	* Example:
	*  "int err = forward_kinematics(solution_matrix, (double*)dh, n_links, joint_angles_deg)"
	*
	*/
	
	if (!dh_array || !solution_matrix || !joint_angles_deg) return ERR_NULL;
	if (n_links < 1 || n_links > DEFAULT_N_LINKS) return ERR_BADARGS;

	double right_matrix[16];
	for (size_t i = 0; i < n_links; i++) {
		if (fill_transformation_matrix(right_matrix, i, dh_array, joint_angles_deg[i], n_links)) {
			return ERR_INNNER;
		}
		if (mat4_mul(solution_matrix, right_matrix)) {
			return ERR_INNNER;
		}
	}

	return ERR_OK;
}

int mat4_mul(double *left_matrix, const double *right_matrix) {
	/**
	* mat4_mul
	*
	* Briefly:
	*   This function do nothing unless calculation of "A_12=A_1*A_2".
	*   Where A_i is supposed to be 4x4 transformation matrix.
	*
	* Preconditions:
	*   - left_matrix != NULL
	*   - right_matrix != NULL
	*   - both left and right matrix is 4x4 flat array type double
	*
	* Postconditions:
	*   - if return 0 it means that now left matrix is the result of multiplication
	*   - in case of error it returns not a zero code and left matrix may be corrupted
	*
	* Responds:
	*   0 - OK
	*  -1 - NULL pointer
	*  -2 - invalid index or invalid memory layout (caller error) 
	*
	* Side effects:
	*   - creates result_matrix for techical uses
	*
	* Complexity:
	*   - Time O(1)
	*   - Memory O(1)
	*
	* Example:
	*  "int err = mat4_mul(A_1, A_2)"
	*/

	// VERIFICATION
	
	if (!left_matrix || !right_matrix) return ERR_NULL; 

	// MULTIPLICATE
	
	double result_matrix[16];
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			result_matrix[i * 4 + j] = 0;
			for (size_t k = 0; k < 4; k++) {
				result_matrix[i * 4 + j] += left_matrix[i * 4 + k] * right_matrix[k * 4 + j];
			}
		}
	}

	// COPY
	
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			left_matrix[i * 4 + j] = result_matrix[i * 4 + j];
		}
	}
	
	return ERR_OK;
}

int fill_transformation_matrix(double *trans_matrix, const size_t index,
		              const double *dh_array,
			      const double angle_coordinate,
			      const size_t n_links) {
	/**
	* straight_kinematics_task
	*
	* Briefly: 
	*   Fills tranformation matrix A_i that reduce coordinates of i-links end to i-1-links end.
	*   This matrix needs Denavit-Hartemberg parametrs of i-link and q_i angle coordinate.
	*
	* Preconditions:
	*   - dh_array != NULL
	*   - dh_array is a flat array n_links * DH_COLS (double)
	*   - index > -1
	*   - trans_matrix != NULL
	*   - trans_matrix is a flat array 4x4 (16 doubles)
	*   - angle_coordinate is a valid double number containing q_i angle coordinate in degrees
	*
	* Postconditions:
	*   - if returns 0 array is filled according to format
	*   - in case of error it returns not a zero code and array may be partially filled
	*
	* Responds:
	*   0 - OK
	*  -1 - NULL pointer
	*  -2 - invalid index or invalid memory layout (caller error) 
	*
	* Side effects:
	*  - memory modification on trans_matrix addreess
	*
	* Complexity:
	*  - Time: O(1)
	*  - Memory: O(1)
	*
	* Example:
	*  "int err = fill_transformation_matrix(trans_matrix, i, (double*)dh, q_i); "
	*
	*/
        
	if (!dh_array || !trans_matrix) return ERR_NULL;
	if (index >= n_links) return ERR_BADARGS;

	double q_i = dh_array[index * DH_COLS + 3] + DEG_TO_RAD(angle_coordinate);

        trans_matrix[0 * 4 + 0] = cos(q_i);
	trans_matrix[0 * 4 + 1] = -sin(q_i) * cos(dh_array[index * DH_COLS + 1]);
	trans_matrix[0 * 4 + 2] = sin(q_i) * sin(dh_array[index * DH_COLS + 1]);
	trans_matrix[0 * 4 + 3] = dh_array[index * DH_COLS + 0] * cos(q_i);

	trans_matrix[1 * 4 + 0] = sin(q_i);
	trans_matrix[1 * 4 + 1] = cos(q_i)  * cos(dh_array[index * DH_COLS + 1]);
	trans_matrix[1 * 4 + 2] = -cos(q_i) * sin(dh_array[index * DH_COLS + 1]);
	trans_matrix[1 * 4 + 3] = dh_array[index * DH_COLS + 0] * sin(q_i);

	trans_matrix[2 * 4 + 0] = 0;
	trans_matrix[2 * 4 + 1] = sin(dh_array[index * DH_COLS + 1]);
	trans_matrix[2 * 4 + 2] = cos(dh_array[index * DH_COLS + 1]);
	trans_matrix[2 * 4 + 3] = dh_array[index * DH_COLS + 2];

	trans_matrix[3 * 4 + 0] = 0;
	trans_matrix[3 * 4 + 1] = 0;
	trans_matrix[3 * 4 + 2] = 0;
	trans_matrix[3 * 4 + 3] = 1; 
	
	return ERR_OK;
}

int fill_dh_parameters(double *dh_array, size_t n_links,
		       const double *links_lengths) {
	/**
 	* fill_dh_parameters
 	*
 	* Briefly:
 	*   Fills table of Denavit-Hartenberg parameters for manipulator with N_LINKS (default 6).
 	*
 	* Array format:
 	   dh_array[i][0] = a_i       (link length)         [units: meters]
 	*   dh_array[i][1] = alpha_i   (link twist)          [units: radians]
 	*   dh_array[i][2] = d_i       (link offset)         [units: meters]
 	*   dh_array[i][3] = theta_i   (joint angle)         [units: radians]
	*
	* Preconditions:
	*   - dh_array != NULL
	*   - links_lengths != NULL
	*   - caller provide, taht array has at least (n_links x 4) elements
	*
	* Postconditions:
	*   - if returns 0 array is filled according to format
	*   - in case of error it returns not a zero code and array may be partially filled
	*
	* Responds:
	*   0 - OK
	*  -1 - NULL pointer
	*  -2 - invalid n_links or invalid memory layout (caller error)
	*
	* Side effects:
	*   - memory modification on address array
	*
	* Complexity:
	*   - Time: O(1)
	*   - Memory: O(1)
	*
	* Example:
 	*  "double dh[6][4];
 	*   int err = fill_dh_parameters((double*)dh, 6, links_lengths);"
 	*/

	if (!dh_array || !links_lengths) return ERR_NULL;
	if (n_links < 1) return ERR_BADARGS;
	if (n_links != DEFAULT_N_LINKS) return ERR_BADARGS;

	dh_array[0*DH_COLS + 0] = 0.0f;
	dh_array[1*DH_COLS + 0] = links_lengths[1];
	dh_array[2*DH_COLS + 0] = links_lengths[2];
	dh_array[3*DH_COLS + 0] = 0.0f;
	dh_array[4*DH_COLS + 0] = 0.0f;
	dh_array[5*DH_COLS + 0] = 0.0f;

	dh_array[0*DH_COLS + 1] = DEG_TO_RAD(90.0f);
	dh_array[1*DH_COLS + 1] = DEG_TO_RAD(0.0f);
	dh_array[2*DH_COLS + 1] = DEG_TO_RAD(0.0f);
	dh_array[3*DH_COLS + 1] = DEG_TO_RAD(-90.0f);
	dh_array[4*DH_COLS + 1] = DEG_TO_RAD(90.0f);
	dh_array[5*DH_COLS + 1] = DEG_TO_RAD(0.0f);

	dh_array[0*DH_COLS + 2] = links_lengths[0];
	dh_array[1*DH_COLS + 2] = 0.0f;
	dh_array[2*DH_COLS + 2] = 0.0f;
	dh_array[3*DH_COLS + 2] = 0.0f;
	dh_array[4*DH_COLS + 2] = 0.0f;
	dh_array[5*DH_COLS + 2] = links_lengths[5];

	dh_array[0*DH_COLS + 3] = DEG_TO_RAD(-90.0f);
	dh_array[1*DH_COLS + 3] = DEG_TO_RAD(0.0f);
	dh_array[2*DH_COLS + 3] = DEG_TO_RAD(0.0f);
	dh_array[3*DH_COLS + 3] = DEG_TO_RAD(0.0f);
	dh_array[4*DH_COLS + 3] = DEG_TO_RAD(90.0f);
	dh_array[5*DH_COLS + 3] = DEG_TO_RAD(0.0f);
	
	return ERR_OK;
}

int inverse_kinematics_dls(
	double *q_out,
	const double *q_init,
    	const double T_des[16],
    	const double *dh_array,
    	size_t n_links,
    	double tol_pos,
    	double tol_ori,
    	int max_iters,
    	double lambda0,
    	double max_delta_q,
    	const double *joint_min,
    	const double *joint_max,
    	int *iters_out /* optional: pass NULL if unused */) {

	/**
	* inverse_kinematics_dls
 	*
	* Решает обратную задачу кинематики методом Ньютона с демпфированным
	* псевдообращением якобиана (Damped Least Squares).
 	*
 	* Входы:
 	*   - q_init         : указатель на массив длины n_links — начальное приближение (degrees).
 	*   - q_out          : указатель на массив длины n_links — сюда запишется найденное решение (degrees).
 	*   - T_des          : целевая матрица преобразования 4x4 (row-major, double[16]).
 	*   - dh_array       : DH-параметры (flat, row per link): [a, alpha, d, theta_off] length n_links*4.
 	*   - n_links        : число звеньев.
 	*   - tol_pos        : позиционная точность (meters), e.g. 1e-4.
 	*   - tol_ori        : ориентировочная точность (degrees), e.g. 1e-3.
 	*   - max_iters      : максимальное число итераций, e.g. 100.
 	*   - lambda0        : начальное демпфирование (например 1e-2).
 	*   - max_delta_q    : максимальное изменение по норме за итерацию (degrees), e.g. 0.1.
 	*   - joint_min, joint_max : (опционально) массивы длины n_links с ограничениями суставов (degrees).
 	*
 	* Выходы:
 	*   - q_out содержит найденный вектор q
 	*   - *iters_out (опционально) — число итераций, которое потребовалось.
 	*
 	* Возвращаемые коды:
 	*   ERR_OK (0)          — успех (ошибка < tolerances).
 	*   ERR_NO_CONVERGE     — не сошлось.
 	*   ERR_BADARGS         — NULL/p	size etc.
 	*   ERR_SINGULAR        — сингулярность/нелинейная неустранимая ошибка.
 	*
 	* Замечания:
 	*   - Все углы во входе и выходе — в градусах и переводятся в радианы в коде.
 	*   - Внутри используется double и size_t.
 	*   - Функция сама вызывает forward_kinematics и build_jacobian на каждой итерации.
 	*/
        if (!q_out || !q_init || !T_des || !dh_array) return ERR_NULL;
    	if (n_links != DEFAULT_N_LINKS) return ERR_BADARGS;

    	double q_rad[DEFAULT_N_LINKS];
    	double T_cur[16];
    	double J[6 * DEFAULT_N_LINKS];
    	double error[6];
    	int iter;

    	// Convert initial guess to radians
    	for (size_t i = 0; i < n_links; i++) {
    	    q_rad[i] = DEG_TO_RAD(q_init[i]);
    	}
	
   	for (iter = 0; iter < max_iters; iter++) {
		// Reset T_cur to identity
		memset(T_cur, 0, 16 * sizeof(double));
		for (int i = 0; i < 4; i++) T_cur[i*4 + i] = 1.0;
	
 	        // Compute current transformation matrix
	        int err_fk = forward_kinematics(T_cur, dh_array, n_links, q_init);
 	        if (err_fk != ERR_OK) return err_fk;
	
	        // Compute position error
 	        error[0] = T_des[3] - T_cur[3];
	        error[1] = T_des[7] - T_cur[7];
	        error[2] = T_des[11] - T_cur[11];

	       	// Compute orientation error (simplified)
		error[3] = T_des[2] - T_cur[2];  // dx
        	error[4] = T_des[6] - T_cur[6];  // dy
      	        error[5] = T_des[10] - T_cur[10]; // dz

        	// Check convergence
        	double pos_error = sqrt(error[0]*error[0] + error[1]*error[1] + error[2]*error[2]);
        	double ori_error = sqrt(error[3]*error[3] + error[4]*error[4] + error[5]*error[5]);
	
        	if (pos_error < tol_pos && ori_error < tol_ori) {
        	    if (iters_out) *iters_out = iter;
        	    // Convert result to degrees
        	    for (size_t i = 0; i < n_links; i++) {
        	        q_out[i] = q_rad[i] * 180.0 / M_PI;
        	    }
        	    return ERR_OK;
        	}
	
	        // Build Jacobian
	        if (build_jacobian(J, q_rad, dh_array, n_links, T_cur) != ERR_OK) {
	            return ERR_INNNER;
	        }
	
	        // Compute JT
	        double JT[DEFAULT_N_LINKS * 6];
	        mat6_transpose(J, JT);
	
	        // Compute JTJ + lambda^2*I
	        double JTJ[36], LI[36], JTJ_plus_lambda2I[36];
	        mat6_mul(JT, J, JTJ);
	        mat6_identity(LI);
	        mat6_scale(LI, lambda0*lambda0, LI);
	        mat6_add(JTJ, LI, JTJ_plus_lambda2I);
	
	        // Solve DLS
	        double delta_q_rad[DEFAULT_N_LINKS];
	        if (solve_dls(JTJ_plus_lambda2I, JT, error, delta_q_rad, n_links) != ERR_OK) {
	            return ERR_SINGULAR;
	        }
	
	        // Update joint angles
	        for (size_t i = 0; i < n_links; i++) {
	            q_rad[i] += delta_q_rad[i];
	            // Apply joint limits
	            if (joint_min && joint_max) {
	                double deg_angle = q_rad[i] * 180.0 / M_PI;
	                if (deg_angle < joint_min[i]) q_rad[i] = DEG_TO_RAD(joint_min[i]);
	                if (deg_angle > joint_max[i]) q_rad[i] = DEG_TO_RAD(joint_max[i]);
		    }
		}
	}
	if (iters_out) *iters_out = max_iters;
	return ERR_NO_CONVERGE;
}


int build_jacobian(double *J, const double *q_rad, const double *dh_array,
                  size_t n_links, const double *T_current) {
    if (!J || !q_rad || !dh_array || !T_current) return ERR_NULL;

    // Simplified Jacobian calculation for 6DOF manipulator
    // This should be replaced with proper analytical Jacobian
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (i == j) {
                J[i*6 + j] = 1.0; // Identity approximation
            } else {
                J[i*6 + j] = 0.0;
            }
        }
    }
    return ERR_OK;
}

int mat6_mul(double *A, const double *B, double *result) {
    if (!A || !B || !result) return ERR_NULL;
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            result[i*6 + j] = 0;
            for (int k = 0; k < 6; k++) {
                result[i*6 + j] += A[i*6 + k] * B[k*6 + j];
            }
        }
    }
    return ERR_OK;
}

int mat6_transpose(const double *A, double *result) {
    if (!A || !result) return ERR_NULL;
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            result[j*6 + i] = A[i*6 + j];
        }
    }
    return ERR_OK;
}

int mat6_add(double *A, const double *B, double *result) {
    if (!A || !B || !result) return ERR_NULL;
    
    for (int i = 0; i < 36; i++) {
        result[i] = A[i] + B[i];
    }
    return ERR_OK;
}

int mat6_scale(double *A, double scalar, double *result) {
    if (!A || !result) return ERR_NULL;
    
    for (int i = 0; i < 36; i++) {
        result[i] = A[i] * scalar;
    }
    return ERR_OK;
}

int mat6_identity(double *I) {
    if (!I) return ERR_NULL;
    
    memset(I, 0, 36 * sizeof(double));
    for (int i = 0; i < 6; i++) {
        I[i*6 + i] = 1.0;
    }
    return ERR_OK;
}

int solve_dls(const double *JTJ_plus_lambda2I, const double *JT, 
              const double *error, double *delta_q, size_t n_links) {
    // Simplified solver - should be replaced with proper linear algebra
    for (size_t i = 0; i < n_links; i++) {
        delta_q[i] = error[i] * 0.01; // Simple scaling
    }
    return ERR_OK;
}
