#ifndef OBSERVER_CONTROLLER_H_
#define OBSERVER_CONTROLLER_H_

#include <stdbool.h>
#include "project_specific.h"

extern const float A[N_STATES][N_STATES];
extern const float B[N_STATES][N_STATES];
extern const float C[N_STATES][N_STATES];

extern const float K[N_STATES][N_STATES];

extern const float Q;
extern const float R;

extern float control_output_process(const float computed_output, const float x_hat[N_STATES], const float timestep);

void observer_init(float x_hat[N_STATES], const float timestep);

bool covariance_matrix_step(void);

void kf_a_priori_state_estimate(const float A_minus_BK[N_STATES][N_STATES], const float timestep, const bool enable, float x_hat[N_STATES]);

void kf_a_posteriori_state_estimate(const float measurement[N_STATES], const float L[N_STATES][N_STATES], const float C[N_STATES][N_STATES], float x_hat[N_STATES]);

void observer_step(const float measurement[N_STATES], const bool enable, float x_hat[N_STATES]);

float control_output(const float x_hat[N_STATES], const float timestep);

#endif // OBSERVER_CONTROLLER_H_
