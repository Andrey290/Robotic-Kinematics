//---LIBRARIES---//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//---HEADEDRS---//

#include "../inc/globalconstants.h"
#include "../inc/kinematic.h"

int main(void) {

	// forward kinematics task
	
	double links_lengths[DEFAULT_N_LINKS] = {
		DEFAULT_L1,
		DEFAULT_L2,
		DEFAULT_L3,
		DEFAULT_L4,
		DEFAULT_L5,
		DEFAULT_L6};
	double joint_angles_deg[DEFAULT_N_LINKS] = {
		30.0f, 30.0f, -30.0f, 10.0f, 10.0f, 10.0f
	};

	double dh[DEFAULT_N_LINKS][DH_COLS];
 	int err_dh = fill_dh_parameters((double*)dh, DEFAULT_N_LINKS, links_lengths);
	
	for (size_t i = 0; i < DEFAULT_N_LINKS; ++i) {
		printf("Link %zu: a=%.3f, alpha=%.3f, d=%.3f, theta=%.3f\n",
				i,
				dh[i][0], dh[i][1], dh[i][2], dh[i][3]);
	}
        
	double solution_matrix[16] = {0};
	for (int i = 0; i < 4; i++) solution_matrix[i*4 + i] = 1.0;

	int err_fk = forward_kinematics(solution_matrix, (double*)dh, (size_t)DEFAULT_N_LINKS, joint_angles_deg);
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			printf(" %.3f |", solution_matrix[i * 4 + j]);
		}
		printf("\n");
	}

	// Inverse kinematics task

	double q_init[DEFAULT_N_LINKS] = {0};
    	double q_result[DEFAULT_N_LINKS] = {0};
    	double joint_min[DEFAULT_N_LINKS] = {-180, -180, -180, -180, -180, -180};
    	double joint_max[DEFAULT_N_LINKS] = {180, 180, 180, 180, 180, 180};
    
    	int iters;
    	int err_ik = inverse_kinematics_dls(
        	q_result, q_init, solution_matrix,
        	(double*)dh, DEFAULT_N_LINKS,
        	1e-4, 1e-3, 100, 1e-2, 0.1,
        	joint_min, joint_max, &iters
    	);

	// Вывод матрицы
	printf("For inverse kinematics task:");
	printf("q = {");
	for (size_t i = 0; i < DEFAULT_N_LINKS; ++i) {
		printf(" %.3f ", q_result[i]);
	}
	printf("}\n");

	return EXIT_SUCCESS;
}


