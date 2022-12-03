#include <math.h>
#include <stdio.h>
#include "observer_controller.h"
#include "matrix_operations.h"

#ifndef UNIT_TEST
#define STATIC
#else
#define STATIC
#endif

void observer_init(kf_input_S* kf_input, kf_states_S* kf_states)
{
	vector_initialize(kf_states->x_hat, 0.0f);
	
	matrix_initialize(kf_states->L, 0.0f);
	matrix_initialize(kf_states->F, 0.0f);
	matrix_initialize(kf_states->F_transpose, 0.0f);
	matrix_initialize(kf_states->Q_matrix, 0.0f);
	matrix_initialize(kf_states->R_matrix, 0.0f);
	matrix_initialize(kf_states->C_transpose, 0.0f);
	matrix_initialize(kf_states->P, 0.0f);
	matrix_initialize(kf_states->A_minus_BK, 0.0f);
	
	identity(kf_states->I);
	
	// Discretize state transition matrix, ie., F = I + dt*A
	float A_discrete[N_STATES][N_STATES] = { { 0.0f } };
	matrix_scale((const float (*)[N_STATES])kf_input->A, kf_input->timestep, A_discrete);
	
	matrix_sum((const float (*)[N_STATES])A_discrete, 
			   (const float (*)[N_STATES])kf_states->I, 
			   kf_states->F);
	
	matrix_transpose((const float (*)[N_STATES])kf_states->F, kf_states->F_transpose);
	matrix_transpose((const float (*)[N_STATES])kf_input->C, kf_states->C_transpose);
	
	float B_transpose[N_STATES][N_STATES] = { { 0.0f } };
	matrix_transpose((const float (*)[N_STATES])kf_input->B, B_transpose);
	
	float B_B_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])kf_input->B, 
						   (const float (*)[N_STATES])B_transpose, 
						   B_B_transpose);
	
	matrix_scale((const float (*)[N_STATES])B_B_transpose, kf_input->Q, kf_states->Q_matrix);
	
	matrix_scale((const float (*)[N_STATES])kf_states->Q_matrix, (kf_input->timestep * kf_input->timestep), kf_states->Q_matrix);
	matrix_scale((const float (*)[N_STATES])kf_states->I, kf_input->R, kf_states->R_matrix);
	
	float B_K[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])kf_input->B, 
						   (const float (*)[N_STATES])kf_input->K, 
						   B_K);
	
	matrix_diff((const float (*)[N_STATES])kf_input->A, 
				(const float (*)[N_STATES])B_K, 
				kf_states->A_minus_BK);
				
	matrix_scale((const float (*)[N_STATES])kf_states->A_minus_BK, kf_input->timestep, kf_states->A_minus_BK);
}

bool covariance_matrix_step(kf_input_S* kf_input, kf_states_S* kf_states)
{
	bool inverse_valid = true;

	// P = (F * P * F') + Q
	float F_P[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])kf_states->F, 
						   (const float (*)[N_STATES])kf_states->P, 
						   F_P);
	
	float F_P_F_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])F_P, 
						   (const float (*)[N_STATES])kf_states->F_transpose, 
						   F_P_F_transpose);
	
	matrix_sum((const float (*)[N_STATES])F_P_F_transpose, 
			   (const float (*)[N_STATES])kf_states->Q_matrix, 
			   kf_states->P);
	
	// S = (C * P * C') + R
	float C_P[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])kf_input->C, 
						   (const float (*)[N_STATES])kf_states->P, 
						   C_P);
	
	float C_P_C_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])C_P, 
						   (const float (*)[N_STATES])kf_states->C_transpose, 
						   C_P_C_transpose);
	
	float S[N_STATES][N_STATES] = {{0.0f}};
	matrix_sum((const float (*)[N_STATES])C_P_C_transpose, 
			   (const float (*)[N_STATES])kf_states->R_matrix, 
			   S);
	
	// L = P * C' * S_inverse
	float S_inverse[N_STATES][N_STATES] = {{0.0f}};
	
	inverse_valid = matrix_inverse_cholesky((const float (*)[N_STATES])S, S_inverse);
	
	float P_C_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])kf_states->P, 
						   (const float (*)[N_STATES])kf_states->C_transpose, 
						   P_C_transpose);
	
	matrix_matrix_multiply((const float (*)[N_STATES])P_C_transpose, 
						   (const float (*)[N_STATES])S_inverse, 
						   kf_states->L);
	
	// P_next = (I - (L * C)) * P;
	float L_C[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])kf_states->L, 
						   (const float (*)[N_STATES])kf_input->C, 
						   L_C);
	
	float I_minus_L_C[N_STATES][N_STATES] = {{0.0f}};
	matrix_diff((const float (*)[N_STATES])kf_states->I, 
				(const float (*)[N_STATES])L_C, 
				I_minus_L_C);
	
	float P_next[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply((const float (*)[N_STATES])I_minus_L_C, 
						   (const float (*)[N_STATES])kf_states->P, 
						   P_next);
	
	matrix_assign((const float (*)[N_STATES])P_next, kf_states->P);

	return inverse_valid;
}

void kf_a_priori_state_estimate(kf_input_S* kf_input, const bool enable, kf_states_S* kf_states)
{
	float prediction[N_STATES] = {0.0f};
	matrix_vector_multiply((const float (*)[N_STATES])kf_states->A_minus_BK, 
						   (const float*)kf_states->x_hat, 
						   prediction);
   
	vector_sum((const float*)kf_states->x_hat, (const float*)prediction, kf_states->x_hat);
	
	vector_scale((const float*)kf_states->x_hat, (enable ? 1.0f : 0.0f), kf_states->x_hat);
}

void kf_a_posteriori_state_estimate(const float measurement[N_STATES], kf_input_S* kf_input, kf_states_S* kf_states)
{	
	// error = y - C * x_hat_prev
	float C_times_x_hat[N_STATES] = {0.0f};
	matrix_vector_multiply((const float (*)[N_STATES])kf_input->C, 
						   (const float*)kf_states->x_hat, 
						   C_times_x_hat);

	float error[N_STATES] = {0.0f};
	vector_diff((const float*)measurement, 
				(const float*)C_times_x_hat, 
				error);
    
	// correction = L * error
	float correction[N_STATES] = {0.0f};
	matrix_vector_multiply((const float (*)[N_STATES])kf_states->L, 
						   (const float*)error, 
						   correction);
						   
	vector_sum((const float*)kf_states->x_hat, correction, kf_states->x_hat);
}

void observer_step(const float measurement[N_STATES], const bool enable, kf_input_S* kf_input, kf_states_S* kf_states)
{	
	kf_a_priori_state_estimate(kf_input, enable, kf_states);
	
	kf_a_posteriori_state_estimate(measurement, kf_input, kf_states);
}

float control_output(const float x_hat[N_STATES], const float timestep, kf_input_S* kf_input)
{	
	const float control_output = -1.0f * dot_product(kf_input->K[0], x_hat);
	
	const float control_output_final = control_output_process(control_output, x_hat, timestep);
	
	return control_output_final;
}
