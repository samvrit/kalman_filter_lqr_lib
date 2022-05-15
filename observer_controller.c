#include <math.h>
#include <stdio.h>
#include "observer_controller.h"
#include "matrix_operations.h"

#ifndef UNIT_TEST
#define STATIC static
#else
#define STATIC
#endif

STATIC float L[N_STATES][N_STATES] = {{0.0f}};
STATIC float F[N_STATES][N_STATES] = {{0.0f}};
STATIC float F_transpose[N_STATES][N_STATES] = {{0.0f}};
STATIC float Q_matrix[N_STATES][N_STATES] = {{0.0f}};
STATIC float R_matrix[N_STATES][N_STATES] = {{0.0f}};
STATIC float B_transpose[N_STATES][N_STATES] = {{0.0f}};
STATIC float C_transpose[N_STATES][N_STATES] = {{0.0f}};
STATIC float P[N_STATES][N_STATES] = {{0.0f}};
STATIC float A_minus_BK[N_STATES][N_STATES] = {{0.0f}};
STATIC float I[N_STATES][N_STATES] = {{0.0f}};

void observer_init(float x_hat[N_STATES], const float timestep)
{
	vector_initialize(x_hat, 0.0f);
	
	matrix_initialize(L, 0.0f);
	matrix_initialize(F, 0.0f);
	matrix_initialize(F_transpose, 0.0f);
	matrix_initialize(Q_matrix, 0.0f);
	matrix_initialize(R_matrix, 0.0f);
	matrix_initialize(B_transpose, 0.0f);
	matrix_initialize(C_transpose, 0.0f);
	matrix_initialize(P, 0.0f);
	matrix_initialize(A_minus_BK, 0.0f);
	
	identity(I);
	
	// Discretize state transition matrix, ie., F = I + dt*A
	matrix_scale(A, timestep, F);
	matrix_sum(F, I, F);
	
	matrix_transpose(F, F_transpose);
	matrix_transpose(B, B_transpose);
	matrix_transpose(C, C_transpose);
	
	float B_B_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(B, B_transpose, B_B_transpose);
	
	matrix_scale(B_B_transpose, Q, Q_matrix);
	matrix_scale(Q_matrix, timestep * timestep, Q_matrix);
	matrix_scale(I, R, R_matrix);
	
	float B_K[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(B, K, B_K);
	
	matrix_diff(A, B_K, A_minus_BK);
	matrix_scale(A_minus_BK, timestep, A_minus_BK);
}

void covariance_matrix_step(void)
{
	// P = (F * P * F') + Q
	float F_P[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(F, P, F_P);
	
	float F_P_F_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(F_P, F_transpose, F_P_F_transpose);
	
	matrix_sum(F_P_F_transpose, Q_matrix, P);
	
	// S = (C * P * C') + R
	float C_P[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(C, P, C_P);
	
	float C_P_C_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(C_P, C_transpose, C_P_C_transpose);
	
	float S[N_STATES][N_STATES] = {{0.0f}};
	matrix_sum(C_P_C_transpose, R_matrix, S);
	
	// L = P * C' * S_inverse
	float S_inverse[N_STATES][N_STATES] = {{0.0f}};
	const bool cholesky_inverse = matrix_inverse_cholesky(S, S_inverse);
	
	// if cholesky decomposition fails because of numerical error, use the general inversion technique
	if (!cholesky_inverse)
	{
		matrix_inverse(S, S_inverse);
	}
	
	float P_C_transpose[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(P, C_transpose, P_C_transpose);
	
	matrix_matrix_multiply(P_C_transpose, S_inverse, L);
	
	// P_next = (I - (L * C)) * P;
	float L_C[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(L, C, L_C);
	
	float I_minus_L_C[N_STATES][N_STATES] = {{0.0f}};
	matrix_diff(I, L_C, I_minus_L_C);
	
	float P_next[N_STATES][N_STATES] = {{0.0f}};
	matrix_matrix_multiply(I_minus_L_C, P, P_next);
	
	matrix_assign(P_next, P);
}	

void observer_step(const float measurement[N_STATES], const bool enable, float x_hat[N_STATES])
{	
	// x_hat[k] = x_hat[k-1] + Ts * ((A-B*K)*x_hat[k-1] + L*(y - C*x_hat[k-1]))
	
	// error = y - C * x_hat_prev
	float C_times_x_hat[N_STATES] = {0.0f};
	matrix_vector_multiply(C, x_hat, C_times_x_hat);

	float error[N_STATES] = {0.0f};
	vector_diff(measurement, C_times_x_hat, error);
	
	// correction = L * error
	float correction[N_STATES] = {0.0f};
	matrix_vector_multiply(L, error, correction);
	
	// prediction = (A-B*K) * x_hat
	float prediction[N_STATES] = {0.0f};
	matrix_vector_multiply(A_minus_BK, x_hat, prediction);

	float prediction_plus_correction[N_STATES] = {0.0f};
	vector_sum(prediction, correction, prediction_plus_correction);
	
	for (int i = 0; i < N_STATES; i++)
	{
		x_hat[i] += prediction_plus_correction[i];
		x_hat[i] = enable ? x_hat[i] : 0.0f;
	}
}

float control_output(const float x_hat[N_STATES], const float timestep)
{	
	const float control_output = -1.0f * dot_product(K[0], x_hat);
	
	const float control_output_final = control_output_process(control_output, x_hat, timestep);
	
	return control_output_final;
}
