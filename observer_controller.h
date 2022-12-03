#ifndef OBSERVER_CONTROLLER_H_
#define OBSERVER_CONTROLLER_H_

#include <stdbool.h>
#include "project_specific.h"

typedef struct
{
	float A[N_STATES][N_STATES];
    float B[N_STATES][N_STATES];
    float C[N_STATES][N_STATES];
    float K[N_STATES][N_STATES];

    float Q;
    float R;

    float timestep;
} kf_input_S;

typedef struct
{
    float L[N_STATES][N_STATES];
    float F[N_STATES][N_STATES];
    float F_transpose[N_STATES][N_STATES];
    float Q_matrix[N_STATES][N_STATES];
    float R_matrix[N_STATES][N_STATES];
    float B_transpose[N_STATES][N_STATES];
    float C_transpose[N_STATES][N_STATES];
    float P[N_STATES][N_STATES];
    float A_minus_BK[N_STATES][N_STATES];
    float I[N_STATES][N_STATES];

    float x_hat[N_STATES];
} kf_states_S;

extern float control_output_process(const float computed_output, const float x_hat[N_STATES], const float timestep);

void observer_init(kf_input_S* kf_input, kf_states_S* kf_states);

bool covariance_matrix_step(kf_input_S* kf_input, kf_states_S* kf_states);

void kf_a_priori_state_estimate(kf_input_S* kf_input, const bool enable, kf_states_S* kf_states);

void kf_a_posteriori_state_estimate(const float measurement[N_STATES], kf_input_S* kf_input, kf_states_S* kf_states);

void observer_step(const float measurement[N_STATES], const bool enable, kf_input_S* kf_input, kf_states_S* kf_states);

float control_output(const float x_hat[N_STATES], const float timestep, kf_input_S* kf_input);

#endif // OBSERVER_CONTROLLER_H_
