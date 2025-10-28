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

    double q_rad[n_links]; // 6
    double T_cur[TRANSFORM_MATRIX_SIZE]; // 16
    double J[6 * n_links]; // 36
    double error[6]; // difference
	int iter;

   	// Convert initial guess to radians
   	for (size_t i = 0; i < n_links; i++) {
   	    q_rad[i] = DEG_TO_RAD(q_init[i]);
   	}

	
   	for (iter = 0; iter < max_iters; iter++) {

		if (iter == 0) {
			printf("\n=== JACOBIAN VALIDATION ===\n");
			double J_num[36] = {0};
			double epsilon = 1e-6;
			
			// Текущее состояние
			double T_base[16];
			double q_temp[6];
			memcpy(q_temp, q_rad, 6 * sizeof(double));
			
			for (int j = 0; j < 6; j++) {
				// +epsilon
				q_temp[j] += epsilon;
				double q_deg_plus[6];
				for (int k = 0; k < 6; k++) q_deg_plus[k] = RAD_TO_DEG(q_temp[k]);
				
				double T_plus[16] = {0};
				for (int i = 0; i < 4; i++) T_plus[i*4 + i] = 1.0;
				forward_kinematics(T_plus, dh_array, 6, q_deg_plus);
				
				// -epsilon  
				q_temp[j] -= 2 * epsilon;
				double q_deg_minus[6];
				for (int k = 0; k < 6; k++) q_deg_minus[k] = RAD_TO_DEG(q_temp[k]);
				
				double T_minus[16] = {0};
				for (int i = 0; i < 4; i++) T_minus[i*4 + i] = 1.0;
				forward_kinematics(T_minus, dh_array, 6, q_deg_minus);
				
				// Восстановим значение
				q_temp[j] += epsilon;
				
				// Численная производная позиции
				for (int i = 0; i < 3; i++) {
					int idx = (i == 0) ? TX : ((i == 1) ? TY : TZ);
					J_num[i * 6 + j] = (T_plus[idx] - T_minus[idx]) / (2 * epsilon);
				}
				
				// Численная производная ориентации (упрощенно)
				for (int i = 0; i < 3; i++) {
					int idx = (i == 0) ? ZX : ((i == 1) ? ZY : ZZ);
					J_num[(i + 3) * 6 + j] = (T_plus[idx] - T_minus[idx]) / (2 * epsilon);
				}
			}
			
			// Сравним аналитический и численный якобианы
			printf("Jacobian comparison (Analytical vs Numerical):\n");
			double max_diff = 0;
			for (int i = 0; i < 36; i++) {
				double diff = fabs(J[i] - J_num[i]);
				if (diff > max_diff) max_diff = diff;
				if (i % 6 == 0) printf("Row %d: ", i/6);
				printf("(%6.3f vs %6.3f) ", J[i], J_num[i]);
				if (i % 6 == 5) printf("\n");
			}
			printf("Max difference: %8.6f\n", max_diff);
		}

        // Convert current q_rad to degrees for forward_kinematics
        double q_deg[DEFAULT_N_LINKS];
        for (size_t i = 0; i < n_links; i++) {
            q_deg[i] = RAD_TO_DEG(q_rad[i]);
        }

		// Reset T_cur to identity
		memset(T_cur, 0, TRANSFORM_MATRIX_SIZE * sizeof(double));
		for (size_t i = 0; i < TRANSFORM_MATRIX_DIM; i++) T_cur[i * TRANSFORM_MATRIX_DIM + i] = 1.0;

        // Compute current transformation matrix
        int err_fk = forward_kinematics(T_cur, dh_array, n_links, q_deg);
        if (err_fk != ERR_OK) return err_fk;

        // Compute position error
        error[0] = T_des[TX] - T_cur[TX];
        error[1] = T_des[TY] - T_cur[TY];
        error[2] = T_des[TZ] - T_cur[TZ];

// CORRECT orientation error computation using angle-axis representation
        double R_des[9] = {T_des[XX], T_des[XY], T_des[XZ],
                          T_des[YX], T_des[YY], T_des[YZ],
                          T_des[ZX], T_des[ZY], T_des[ZZ]};
        
        double R_cur[9] = {T_cur[XX], T_cur[XY], T_cur[XZ],
                          T_cur[YX], T_cur[YY], T_cur[YZ],
                          T_cur[ZX], T_cur[ZY], T_cur[ZZ]};
        
        // Compute orientation error
        double orientation_error[3];
        if (compute_orientation_error(R_des, R_cur, orientation_error) != ERR_OK) {
            return ERR_INNNER;
        }
        
        error[3] = orientation_error[0];
        error[4] = orientation_error[1];
        error[5] = orientation_error[2];

        // Check convergence
        double pos_error = sqrt(error[0]*error[0] + error[1]*error[1] + error[2]*error[2]);
    	double ori_error = sqrt(error[3]*error[3] + error[4]*error[4] + error[5]*error[5]);
	
		printf("Iter %d: pos_error = %e, ori_error = %e\n", iter, pos_error, ori_error);

        if (pos_error < tol_pos && ori_error < tol_ori) {
        	if (iters_out) *iters_out = iter;
        	// Convert result to degrees
        	for (size_t i = 0; i < n_links; i++) {
        	    q_out[i] = RAD_TO_DEG(q_rad[i]);
        	}
        	return ERR_OK;
        }
	
	    // Build Jacobian
	    if (build_jacobian(J, q_rad, dh_array, n_links, T_cur) != ERR_OK) {
	        return ERR_INNNER;
	    }

		// Check Jacobian norm
        double J_norm = 0.0;
        for (int i = 0; i < 36; i++) {
            J_norm += J[i] * J[i];
        }
        J_norm = sqrt(J_norm);
        printf("Iter %3d: Jacobian norm = %8.6f\n", iter, J_norm);
	
	    // Compute JT
	    double JT[n_links * 6];
	    mat6_transpose(J, JT);
	
	    // Compute JTJ + lambda^2*I
	    double JTJ[36], LI[36], JTJ_plus_lambda2I[36];
	    mat6_mul(JT, J, JTJ);
	    mat6_identity(LI);
	    mat6_scale(LI, lambda0*lambda0, LI);
	    mat6_add(JTJ, LI, JTJ_plus_lambda2I);

		// Compute J^T * error
	    double JT_error[6];
	    mat6_mul_vector(JT, error, JT_error, 6, 6);
	    
        // SOLVE: (JTJ + lambda^2*I) * delta_q = J^T * error
	    double delta_q_rad[n_links];
	    if (solve_dls_system(JTJ_plus_lambda2I, JT_error, delta_q_rad, 6, lambda0) != ERR_OK) {
	        return ERR_SINGULAR;
	    }

		// Проверим правую часть уравнения
		printf("Error vector: ");
		for (int i = 0; i < 6; i++) {
			printf("%8.6f ", error[i]);
		}
		printf("\n");

		// Проверим J^T * error
		double JT_error_test[6] = {0};
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				JT_error_test[i] += JT[i * 6 + j] * error[j];
			}
		}

		printf("JT_error (computed): ");
		for (int i = 0; i < 6; i++) {
			printf("%8.6f ", JT_error[i]);
		}
		printf("\n");

		printf("JT_error (test): ");
		for (int i = 0; i < 6; i++) {
			printf("%8.6f ", JT_error_test[i]);
		}
		printf("\n");


		double delta_norm_before_clamp = 0.0;
    	for (size_t i = 0; i < n_links; i++) {
        	delta_norm_before_clamp += delta_q_rad[i] * delta_q_rad[i];
    	}
    	delta_norm_before_clamp = sqrt(delta_norm_before_clamp);
    	printf("Iter %3d: delta_q norm = %8.6f\n", iter, delta_norm_before_clamp);

	    // DEBUG: Проверим вычисление нормы
		double delta_norm = 0.0;
		double debug_norm = 0.0;
		printf("Delta_q values for norm calculation: ");
		for (size_t i = 0; i < n_links; i++) {
			printf("%8.6f ", delta_q_rad[i]);
			debug_norm += delta_q_rad[i] * delta_q_rad[i];
		}
		debug_norm = sqrt(debug_norm);
		printf("\nComputed norm: %8.6f, Printed norm: %8.6f\n", debug_norm, delta_norm);

		// Убедимся, что используем правильную переменную
		delta_norm = debug_norm;
	    
	    double max_delta_rad = DEG_TO_RAD(max_delta_q);
	    if (delta_norm > max_delta_rad) {
	        double scale = max_delta_rad / delta_norm;
	        for (size_t i = 0; i < n_links; i++) {
	            delta_q_rad[i] *= scale;
	        }
	    }

	    // Update joint angles
	    for (size_t i = 0; i < n_links; i++) {
	        q_rad[i] -= delta_q_rad[i];
	        // Apply joint limits
	        if (joint_min && joint_max) {
	            double deg_angle = q_rad[i] * 180.0 / acos(-1.0);
	            if (deg_angle < joint_min[i]) q_rad[i] = DEG_TO_RAD(joint_min[i]);
	            if (deg_angle > joint_max[i]) q_rad[i] = DEG_TO_RAD(joint_max[i]);
		    }
		}
	
		// Debug: print delta_q
		printf("Delta_q: [");
		for (size_t i = 0; i < n_links; i++) {
		    printf("%6.4f ", delta_q_rad[i]);
		}
		printf("]\n");
	
	}


	if (iters_out) *iters_out = max_iters;
	return ERR_NO_CONVERGE;
}


int build_jacobian(double *J, const double *q_rad, const double *dh_array,
                  size_t n_links, const double *T_current) {
	/**
	* build_jacobian
	* 
	* Computes the analytical Jacobian for a 6-DOF serial manipulator
	* 
	* For rotational joints, the Jacobian columns are:
	* J_linear = z_{i-1} × (p_eff - p_{i-1})
	* J_angular = z_{i-1}
	* 
	* Where:
	* - z_{i-1} is the z-axis of joint i-1 in base coordinates
	* - p_eff is the end-effector position in base coordinates  
	* - p_{i-1} is the position of joint i-1 in base coordinates
	*/

    if (!J || !q_rad || !dh_array || !T_current) return ERR_NULL;

    // Initialize Jacobian to zeros
    for (int i = 0; i < JACOBIAN_SIZE; i++) {
        J[i] = 0.0;
    }

    // Compute transformation matrices for each joint
    double T_accumulated[TRANSFORM_MATRIX_SIZE];
    double joint_transforms[DEFAULT_N_LINKS * TRANSFORM_MATRIX_SIZE];
    
    // Start with identity matrix
    for (int i = 0; i < TRANSFORM_MATRIX_SIZE; i++) {
        T_accumulated[i] = (i % (TRANSFORM_MATRIX_DIM + 1) == 0) ? 1.0 : 0.0;
    }

    // Compute forward kinematics for each joint
    for (size_t joint_idx = 0; joint_idx < n_links; joint_idx++) {
        double A_i[TRANSFORM_MATRIX_SIZE];
        
        // Compute transformation matrix for current joint
        if (fill_transformation_matrix(A_i, joint_idx, dh_array, 
                                     RAD_TO_DEG(q_rad[joint_idx]), n_links) != ERR_OK) {
            return ERR_INNNER;
        }
        
        // Accumulate transformation
        if (mat4_mul(T_accumulated, A_i) != ERR_OK) {
            return ERR_INNNER;
        }
        
        // Store joint transformation
        for (int i = 0; i < TRANSFORM_MATRIX_SIZE; i++) {
            joint_transforms[joint_idx * TRANSFORM_MATRIX_SIZE + i] = T_accumulated[i];
        }
    }
	// End-effector position from final transformation
    double p_eff[3] = {T_current[TX], T_current[TY], T_current[TZ]};

    // Compute Jacobian columns for each joint
    for (size_t joint_idx = 0; joint_idx < n_links; joint_idx++) {
        // Get transformation up to this joint
        double *T_joint = &joint_transforms[joint_idx * TRANSFORM_MATRIX_SIZE];
        
        // Joint position
        double p_joint[3] = {T_joint[TX], T_joint[TY], T_joint[TZ]};
        
        // Joint z-axis (rotation axis in base frame)
        double z_joint[3] = {T_joint[ZX], T_joint[ZY], T_joint[ZZ]};
        
        // Vector from joint to end-effector
        double r[3] = {
            p_eff[0] - p_joint[0],
            p_eff[1] - p_joint[1],
            p_eff[2] - p_joint[2]
        };
        
        // Linear part: z_joint × (p_eff - p_joint) [m/rad]
        // This gives the linear displacement per radian of joint rotation
        double J_linear[3];
        J_linear[0] = z_joint[1] * r[2] - z_joint[2] * r[1];  // dX/dq [m/rad]
        J_linear[1] = z_joint[2] * r[0] - z_joint[0] * r[2];  // dY/dq [m/rad]
        J_linear[2] = z_joint[0] * r[1] - z_joint[1] * r[0];  // dZ/dq [m/rad]
        
        // Angular part: z_joint [rad/rad]  
        // This gives the angular displacement per radian of joint rotation
        double J_angular[3] = {z_joint[0], z_joint[1], z_joint[2]};  // [rad/rad]
        
        // Fill Jacobian column for this joint
        // Linear velocity part (position change)
        J[0 * n_links + joint_idx] = J_linear[0];  // ∂X/∂q_j
        J[1 * n_links + joint_idx] = J_linear[1];  // ∂Y/∂q_j
        J[2 * n_links + joint_idx] = J_linear[2];  // ∂Z/∂q_j
        
        // Angular part (orientation change)  
        J[3 * n_links + joint_idx] = J_angular[0]; // ∂ω_x/∂q_j [rad/rad]
        J[4 * n_links + joint_idx] = J_angular[1]; // ∂ω_y/∂q_j [rad/rad]
        J[5 * n_links + joint_idx] = J_angular[2]; // ∂ω_z/∂q_j [rad/rad]
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

// NEW FUNCTION: Compute orientation error using angle-axis representation
int compute_orientation_error(const double *R_des, const double *R_cur, double *error) {
    if (!R_des || !R_cur || !error) return ERR_NULL;
    
    // Compute relative rotation: R_err = R_des * R_cur^T
    double R_cur_T[9];
    mat3_transpose(R_cur, R_cur_T);
    
    double R_err[9];
    mat3_mul(R_des, R_cur_T, R_err);
    
    // Convert to angle-axis representation
    double angle, axis[3];
    if (rotation_matrix_to_angle_axis(R_err, &angle, axis) != ERR_OK) {
        return ERR_INNNER;
    }
    
    // Orientation error = angle * axis
    error[0] = angle * axis[0];
    error[1] = angle * axis[1];
    error[2] = angle * axis[2];
    
    return ERR_OK;
}

/**
 * Proper Damped Least Squares solver
 * Solves: (J^T*J + lambda^2*I) * delta_q = J^T * error
 */
int solve_dls_system(const double *JTJ_plus_lambda2I, const double *JT_error, 
                     double *delta_q, int n, double lambda) {
    if (!JTJ_plus_lambda2I || !JT_error || !delta_q) return ERR_NULL;
    
    // For small systems (6x6), we can use direct Gaussian elimination
    // but we need to ensure numerical stability
    
    double A[6][6];
    double b[6];
    
    // Copy system matrix and right-hand side
    for (int i = 0; i < n; i++) {
        b[i] = JT_error[i];
        for (int j = 0; j < n; j++) {
            A[i][j] = JTJ_plus_lambda2I[i * n + j];
        }
    }
    
    // Gaussian elimination with partial pivoting and condition checking
    for (int i = 0; i < n; i++) {
        // Find pivot row
        int max_row = i;
        double max_val = fabs(A[i][i]);
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > max_val) {
                max_val = fabs(A[k][i]);
                max_row = k;
            }
        }
        
        // Swap rows if necessary
        if (max_row != i) {
            for (int j = 0; j < n; j++) {
                double temp = A[i][j];
                A[i][j] = A[max_row][j];
                A[max_row][j] = temp;
            }
            double temp_b = b[i];
            b[i] = b[max_row];
            b[max_row] = temp_b;
        }
        
        // Check for singular or ill-conditioned system
        if (fabs(A[i][i]) < 1e-12) {
            // Use damping to handle singularity
            A[i][i] = copysign(1e-12, A[i][i]);
        }
        
        // Elimination
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            // Avoid division by very small numbers
            if (fabs(factor) > 1e12) {
                factor = copysign(1e12, factor);
            }
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    
    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        delta_q[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            delta_q[i] -= A[i][j] * delta_q[j];
        }
        
        // Avoid division by zero and limit step size
        if (fabs(A[i][i]) > 1e-12) {
            delta_q[i] /= A[i][i];
        } else {
            delta_q[i] = 0.0;
        }
        
        // Limit individual joint steps to prevent instability
        if (fabs(delta_q[i]) > 0.1) { // ~5.7 degrees
            delta_q[i] = copysign(0.1, delta_q[i]);
        }
    }
    
    return ERR_OK;
}

// NEW FUNCTION: Matrix-vector multiplication for 6x6 * 6x1
int mat6_mul_vector(const double *A, const double *v, double *result, int rows, int cols) {
    if (!A || !v || !result) return ERR_NULL;
    
    for (int i = 0; i < rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += A[i * cols + j] * v[j];
        }
    }
    return ERR_OK;
}

// NEW FUNCTION: Convert rotation matrix to angle-axis representation
int rotation_matrix_to_angle_axis(const double *R, double *angle, double *axis) {
    if (!R || !angle || !axis) return ERR_NULL;
    
    double trace = R[0] + R[4] + R[8];
    double cos_angle = (trace - 1.0) / 2.0;
    
    // Handle numerical precision
    if (cos_angle > 1.0) cos_angle = 1.0;
    if (cos_angle < -1.0) cos_angle = -1.0;
    
    *angle = acos(cos_angle);
    
    // Handle special cases
    if (fabs(*angle) < 1e-10) {
        // Zero rotation
        axis[0] = 1.0;
        axis[1] = 0.0;
        axis[2] = 0.0;
    } else if (fabs(*angle - acos(-1.0)) < 1e-10) {
        // 180 degree rotation
        axis[0] = sqrt((R[0] + 1.0) / 2.0);
        axis[1] = sqrt((R[4] + 1.0) / 2.0);
        axis[2] = sqrt((R[8] + 1.0) / 2.0);
        
        // Determine signs
        if (fabs(axis[0]) > 1e-10) {
            axis[1] = copysign(axis[1], R[1]);
            axis[2] = copysign(axis[2], R[2]);
        } else if (fabs(axis[1]) > 1e-10) {
            axis[2] = copysign(axis[2], R[6]);
        }
        
        // Normalize
        double norm = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
        if (norm > 1e-10) {
            axis[0] /= norm;
            axis[1] /= norm;
            axis[2] /= norm;
        }
    } else {
        // General case
        double sin_angle = sin(*angle);
        axis[0] = (R[7] - R[5]) / (2.0 * sin_angle);
        axis[1] = (R[2] - R[6]) / (2.0 * sin_angle);
        axis[2] = (R[3] - R[1]) / (2.0 * sin_angle);
        
        // Normalize
        double norm = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
        if (norm > 1e-10) {
            axis[0] /= norm;
            axis[1] /= norm;
            axis[2] /= norm;
        }
    }
    
    return ERR_OK;
}

int mat3_transpose(const double *A, double *result) {
    if (!A || !result) return ERR_NULL;
    
    result[0] = A[0]; result[1] = A[3]; result[2] = A[6];
    result[3] = A[1]; result[4] = A[4]; result[5] = A[7];
    result[6] = A[2]; result[7] = A[5]; result[8] = A[8];
    
    return ERR_OK;
}

int mat3_mul(const double *A, const double *B, double *result) {
    if (!A || !B || !result) return ERR_NULL;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i*3 + j] = 0;
            for (int k = 0; k < 3; k++) {
                result[i*3 + j] += A[i*3 + k] * B[k*3 + j];
            }
        }
    }
    return ERR_OK;
}
