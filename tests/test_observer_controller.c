#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "test_observer_controller.h"
#include "observer_controller.h"
#include "matrix_operations.h"

const float A[4][4] = {	{0,    1.0000,         0,         0},
						{0,         0,         0,         0},
						{0,         0,         0,    1.0000},
						{0,         0,    1.9620,         0}};

const float B[4][4] = {	{0.0f, 0.0f, 0.0f, 0.0f}, 
						{1.0f, 0.0f, 0.0f, 0.0f},
						{0.0f, 0.0f, 0.0f, 0.0f}, 
						{0.2f, 0.0f, 0.0f, 0.0f}};

const float C[4][4] = {	{1.0f, 0.0f, 0.0f, 0.0f},
						{0.0f, 0.0f, 0.0f, 0.0f},
						{0.0f, 0.0f, 1.0f, 0.0f},
						{0.0f, 0.0f, 0.0f, 0.0f}};

const float K[4][4] = {	{-1.0000, -5.2796, 97.4833, 73.0332},
						{ 0.0000,  0.0000,  0.0000,  0.0000},
						{ 0.0000,  0.0000,  0.0000,  0.0000},
						{ 0.0000,  0.0000,  0.0000,  0.0000}};

const float Q = 1000.0f;
const float R = 1.0f;

float L_expected[N_STATES][N_STATES] = {	{0.000775f, 0.0f, 0.000143f, 0.0f},
											{0.003109f, 0.0f, 0.000574f, 0.0f},
											{0.000143f, 0.0f, 0.000092f, 0.0f},
											{0.000640f, 0.0f, 0.000216f, 0.0f}};

void display_matrix(const float matrix[N_STATES][N_STATES], const char * matrix_name)
{
	printf("\n%s:\n", matrix_name);
	for(int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			printf("%f\t", matrix[i][j]);
		}
		printf("\n");
	}
}

int main(void)
{
	const float Ts = 1e-4f;
	const int N = 20000;
	const float tolerance = 0.00001f;
	
	kf_input_S kf_input;
	kf_states_S kf_states;
	
	memcpy(kf_input.A, A, N_STATES*N_STATES*sizeof(float));
	memcpy(kf_input.B, B, N_STATES*N_STATES*sizeof(float));
	memcpy(kf_input.C, C, N_STATES*N_STATES*sizeof(float));
	memcpy(kf_input.K, K, N_STATES*N_STATES*sizeof(float));

	kf_input.Q = Q;
	kf_input.R = R;

	kf_input.timestep = Ts;
	
	kf_observer_init(&kf_input, &kf_states);
	
	for(int sim_step = 0; sim_step < N; sim_step++)
	{
		kf_covariance_matrix_step(&kf_input, &kf_states);
	}
	
	const bool equal = matrix_equal_check(kf_states.L, L_expected, tolerance);
	
	printf("\n\nMatrices equal: %d\n\n", equal);
	
	display_matrix(kf_states.L, "Kalman gain");
	
	return 0;
}
