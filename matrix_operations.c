#include <math.h>
#include "matrix_operations.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))

void matrix_vector_multiply(const float A[N_STATES][N_STATES], const float x[N_STATES], float output[N_STATES])
{
    for(int i = 0; i < N_STATES; i++)
	{
	    float sum = 0.0f;
	    for(int j = 0; j < N_STATES; j++)
	    {
	        sum += A[i][j] * x[j];
	    }
	    output[i] = sum;
	}
}

void matrix_matrix_multiply(const float A[N_STATES][N_STATES], const float B[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = 0.0f;
			for (int k = 0; k < N_STATES; k++)
			{
				output[i][j] += (A[i][k] * B[k][j]);
			}
		}
	}
}

float dot_product(const float K[N_STATES], const float x[N_STATES])
{
    float sum = 0.0f;
    for(int i = 0; i < N_STATES; i++)
    {
        sum += K[i] * x[i];
    }
    return sum;
}

void vector_sum(const float x1[N_STATES], const float x2[N_STATES], float output[N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		output[i] = x1[i] + x2[i];
	}
}

void matrix_sum(const float A[N_STATES][N_STATES], const float B[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = A[i][j] + B[i][j];
		}
	}
}

void matrix_diff(const float A[N_STATES][N_STATES], const float B[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = A[i][j] - B[i][j];
		}
	}
}

void vector_diff(const float x1[N_STATES], const float x2[N_STATES], float output[N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		output[i] = x1[i] - x2[i];
	}
}

void vector_scale(const float x[N_STATES], const float a, float output[N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		output[i] = x[i] * a;
	}
}

void matrix_scale(const float A[N_STATES][N_STATES], const float a, float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = A[i][j] * a;
		}
	}
}

void matrix_transpose(const float A[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = A[j][i];
		}
	}
}

void matrix_assign(const float input[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = input[i][j];
		}
	}
}

void vector_initialize(float vector[N_STATES], const float init_val)
{
	for (int i = 0; i < N_STATES; i++)
	{
		vector[i] = init_val;
	}
}

void vector_assign(const float input[N_STATES], float output[N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		output[i] = input[i];
	}
}

void matrix_initialize(float matrix[N_STATES][N_STATES], const float init_val)
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			matrix[i][j] = init_val;
		}
	}
}

void identity(float matrix[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			matrix[i][j] = (i == j) ? 1.0f : 0.0f;
		}
	}
}

bool matrix_equal_check(float matrix1[N_STATES][N_STATES], float matrix2[N_STATES][N_STATES], const float tolerance)
{
	bool equal = true;
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			if ( fabsf(matrix1[i][j] - matrix2[i][j]) > tolerance )
			{
				equal = false;
			}
		}
	}
	return equal;
}

/*-------------------------- Inverse functions -----------------------------*/

static bool cholesky_decomposition(const float A[N_STATES][N_STATES], float R[N_STATES][N_STATES])
{
	bool result = true;
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = i; j < N_STATES; j++)
		{
			R[i][j] = 0.0f;
			float sum = 0.0f;
			for (int k = 0; k < i; k++)
			{
				sum += R[k][i] * R[k][j];
			}
			if (i == j)
			{
				const float diff = A[i][i] - sum;
				if (diff < 0.0f)
				{
					result = false;
					break;
				}
				else
				{
					R[i][i] = sqrtf(diff);
				}
				
			}
			else
			{
				if (R[i][i] == 0.0f)
				{
					result = false;
					break;
				}
				else
				{
					R[i][j] = (1.0f / R[i][i]) * (A[i][j] - sum);
				}
			}
		}
	}
	return result;
}

static bool upper_triangular_inverse(const float R[N_STATES][N_STATES], const float S[N_STATES][N_STATES], float A_inv[N_STATES][N_STATES])
{
	bool result = true;
	for (int i = N_STATES-1; i >= 0; i--)
	{
		for (int j = N_STATES-1; j >= i; j--)
		{
			float sum = 0.0f;
			for(int k = N_STATES-1; k > i; k--)
			{
				sum += R[i][k] * A_inv[k][j];
			}
			if (R[i][i] == 0.0f)
			{
				result = false;
				break;
			}
			else
			{
				A_inv[i][j] = (1.0f / R[i][i]) * (S[i][j] - sum);
			}			
		}
	}

	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < i; j++)
		{
			A_inv[i][j] = A_inv[j][i];
		}
	}
	
	return result;
}

// Based on the paper "Matrix Inversion Using Cholesky Decomposition" by Aravindh Krishnamoorthy and Deepak Menon
bool matrix_inverse_cholesky(const float A[N_STATES][N_STATES], float A_inv[N_STATES][N_STATES])
{
	float R[N_STATES][N_STATES] = {{0.0f}};
    const bool cholesky = cholesky_decomposition(A, R);
	
	if (!cholesky)
		return false;
	
	float S[N_STATES][N_STATES] = {{0.0f}};
	for (int i = 0; i < N_STATES; i++)
	{
		S[i][i] = (1.0f / R[i][i]);
	}
	
	const bool upper_inverse = upper_triangular_inverse((const float (*)[N_STATES])R,
	                                                    (const float (*)[N_STATES])S,
	                                                    A_inv);
    
	return (cholesky && upper_inverse);
}
