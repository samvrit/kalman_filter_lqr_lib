#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include "test_observer_controller.h"
#include "observer_controller.h"
#include "matrix_operations.h"

extern float L[N_STATES][N_STATES];
extern float F[N_STATES][N_STATES];
extern float F_transpose[N_STATES][N_STATES];
extern float Q_matrix[N_STATES][N_STATES];
extern float R_matrix[N_STATES][N_STATES];
extern float B_transpose[N_STATES][N_STATES];
extern float C_transpose[N_STATES][N_STATES];
extern float P[N_STATES][N_STATES];
extern float A_minus_BK[N_STATES][N_STATES];
extern float I[N_STATES][N_STATES];

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
	float x_hat[N_STATES] = {0.0f};
	const float Ts = 1e-4f;
	const int N = 20000;
	const float tolerance = 0.00001f;
	
	observer_init(x_hat, Ts);
	
	for(int sim_step = 0; sim_step < N; sim_step++)
	{
		covariance_matrix_step();
	}
	
	const bool equal = matrix_equal_check(L, L_expected, tolerance);
	
	printf("\n\nMatrices equal: %d\n\n", equal);
	
	display_matrix(L, "Kalman gain");
	
	return 0;
}
