//---LIBRARIES---//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//---HEADEDRS---//

#include "../inc/globalconstants.h"
#include "../inc/kinematic.h"

int main(void) {
    // Forward kinematics task
    double links_lengths[DEFAULT_N_LINKS] = {
        DEFAULT_L1, DEFAULT_L2, DEFAULT_L3,
        DEFAULT_L4, DEFAULT_L5, DEFAULT_L6
    };
    
    double joint_angles_deg[DEFAULT_N_LINKS] = {
        30.0, 30.0, -30.0, 10.0, 10.0, 10.0
    };

    double dh[DEFAULT_N_LINKS][DH_COLS];
    int err_dh = fill_dh_parameters((double*)dh, DEFAULT_N_LINKS, links_lengths);
    
    if (err_dh != ERR_OK) {
        printf("Error filling DH parameters: %d\n", err_dh);
        return EXIT_FAILURE;
    }
    
    printf("DH Parameters:\n");
    for (size_t i = 0; i < DEFAULT_N_LINKS; ++i) {
        printf("Link %zu: a=%.3f, alpha=%.3f, d=%.3f, theta=%.3f\n",
                i, dh[i][0], dh[i][1], dh[i][2], dh[i][3]);
    }
        
    double solution_matrix[TRANSFORM_MATRIX_SIZE] = {0};
    for (int i = 0; i < 4; i++) solution_matrix[i*4 + i] = 1.0;

    int err_fk = forward_kinematics(solution_matrix, (double*)dh, 
                                   DEFAULT_N_LINKS, joint_angles_deg);
    
    if (err_fk != ERR_OK) {
        printf("Error in forward kinematics: %d\n", err_fk);
        return EXIT_FAILURE;
    }
    
    printf("\nForward Kinematics Result:\n");
    for (size_t i = 0; i < 4; i++) {
        for (size_t j = 0; j < 4; j++) {
            printf(" %.3f |", solution_matrix[i * 4 + j]);
        }
        printf("\n");
    }

    // Inverse kinematics task
    double q_init[DEFAULT_N_LINKS] = {20.0, 20.0, -20.0, 5.0, 5.0, 5.0}; // Close initial guess
    double q_result[DEFAULT_N_LINKS] = {0};
    double joint_min[DEFAULT_N_LINKS] = {-180, -180, -180, -180, -180, -180};
    double joint_max[DEFAULT_N_LINKS] = {180, 180, 180, 180, 180, 180};
    
    int iters;
    printf("\nStarting Inverse Kinematics...\n");
    int err_ik = inverse_kinematics_dls(
        q_result, q_init, solution_matrix,
        (double*)dh, DEFAULT_N_LINKS,
        1e-3, 2.0, 100, 0.5, 2.0,  // More relaxed tolerances for testing
        joint_min, joint_max, &iters
    );

    if (err_ik == ERR_OK) {
        printf("Inverse Kinematics SUCCESS after %d iterations\n", iters);
        printf("Result angles: {");
        for (size_t i = 0; i < DEFAULT_N_LINKS; ++i) {
            printf(" %.3f", q_result[i]);
        }
        printf(" }\n");
        
        // Verify by doing forward kinematics with result
        double verification_matrix[TRANSFORM_MATRIX_SIZE] = {0};
        for (int i = 0; i < 4; i++) verification_matrix[i*4 + i] = 1.0;
        
        forward_kinematics(verification_matrix, (double*)dh, DEFAULT_N_LINKS, q_result);
        
        printf("\nVerification - FK with IK result:\n");
        for (size_t i = 0; i < 4; i++) {
            for (size_t j = 0; j < 4; j++) {
                printf(" %.3f |", verification_matrix[i * 4 + j]);
            }
            printf("\n");
        }
    } else {
        printf("Inverse Kinematics FAILED with error: %d\n", err_ik);
    }

    return EXIT_SUCCESS;
}

