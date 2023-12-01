#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "utility.h"
#include "my_fir.h"
#include "log.h"
#include "moving_average.h"
#include "kalmanfilter.h"
#define EXTREMA_TEST 1
#define FIR_FILTER_ENABLE 1
extern kalman_filter_data_s kalman_data;
void my_find_extrema(int16_t const *restrict x, size_t N, int32_t *restrict maxx, int16_t *restrict maxy, int32_t *nmax)
{
	// Set the number of extrema and zero crossings to zero initially
	*nmax = 0;
	// Handle empty array as a special case
	if (N == 0)
	{
		return;
	}
	// Add the ends of the data as both local minima and maxima. These
	// might be changed later by linear extrapolation.
	maxx[0] = 0;
	maxy[0] = x[0];
	(*nmax)++;
	// If we had only one data point this is it
	if (N == 1)
	{
		return;
	}
	// printf("data n is %d", N);
	//  Now starts the main extrema-finding loop. The loop detects points where
	//  the slope of the data changes sign. In the case of flat regions at the
	//  extrema, the center point of the flat region will be considered the
	//  extremal point. While detecting extrema,
	enum slope
	{
		UP,
		DOWN,
		NONE
	};
	enum slope previous_slope = NONE;
	int flat_counter = 0;

	for (size_t i = 0; i < N - 1; i++)
	{
		if (x[i + 1] > x[i])
		{ // Going up
			// if (previous_slope == DOWN)
			// {
			// 	// Was going down before -> local minimum found
			// }
			previous_slope = UP;
			flat_counter = 0;
		}
		else if (x[i + 1] < x[i])
		{ // Going down

			if (previous_slope == UP)
			{
				// Was going up before -> local maximum found
				// if ((last_index - i) > 2)
				{
					// if (((x[i] - last_y) > 5) || ((last_y - x[i]) > 5))
					if (x[i] > 0)
					{
						/* code */
						maxx[*nmax] = (int32_t)(i) - (int32_t)(flat_counter >> 1);
						// maxx[*nmax] = (int32_t)(i);
						maxy[*nmax] = x[i];
						// printf("%d  %d \r\n", maxx[*nmax], x[i]);
						(*nmax)++;
					}
				}
			}

			previous_slope = DOWN;
			flat_counter = 0;
		}
		else
		{ // Staying flat
			flat_counter++;
		}
	}
	// // Add the other end of the data as extrema as well.
	maxx[*nmax] = N - 1;
	maxy[*nmax] = x[N - 1];
	(*nmax)++;

	return;
}
// 寻找数组最小值的下标
int argmin(int *index, int index_len)
{
	int min_index = 0;
	int min = index[0];
	for (int i = 1; i < index_len; i++)
	{
		if (index[i] < min)
		{
			min = index[i];
			min_index = i;
		}
	}
	return min_index;
}
// 寻找极值点函数
//  data是存放数据的数组
// index是存放峰值点下标的数组
// len_index是峰值个数，即index数组长度
void my_AMPD(int16_t const *restrict data, int sizeN, int *index, int *len_index)
{
	int16_t *p_data = (int16_t *)malloc(sizeof(int16_t) * sizeN); // size可以最大为数组长度
	int *arr_rowsum = (int *)malloc(sizeof(int) * sizeN);
	int min_index, max_window_length;
	for (int i = 0; i < sizeN; i++)
	{
		p_data[i] = 0;
	}
	for (int k = 0; k < (int)(sizeN / 2) + 1; k++)
	{
		int row_sum = 0;
		for (int i = k; i <= sizeN - k; i++)
		{
			if ((data[i] > data[i - k]) && (data[i] > data[i + k]))
				row_sum -= 1;
		}
		// *(arr_rowsum + k - 1) = row_sum; //
		arr_rowsum[k - 1] = row_sum;
	}
	// for (int i = 0; i < sizeN/2; i++)
	// {
	// 	printf("%d ", arr_rowsum[i]);
	// }
	// printf("\r\n");
	min_index = argmin(arr_rowsum, sizeN >> 1); // 此处为最大的窗口
	printf("min_index:%d\n", min_index);
	max_window_length = (min_index);
	for (int k = 0; k < max_window_length + 1; k++)
	{
		for (int i = 0; i < sizeN - k; i++)
		{
			if ((data[i] > data[i - k]) && (data[i] > data[i + k]))
				p_data[i] += 1;
		}
	}
	for (int i_find = 0; i_find < sizeN; i_find++)
	{
		if (p_data[i_find] == max_window_length)
		{
			index[*len_index] = i_find;
			(*len_index) += 1;
			// printf("%d ", i_find);
		}
	}
	printf("\r\nindex len %d \r\n", *len_index);

	free(p_data);
	free(arr_rowsum);
}
// int16_t my_abs(int16_t num)
// {
// 	if (num < 0)
// 	{
// 		/* 负数补码先减一, 再取反, 得到原码, 即正数 */
// 		return ~(--num);
// 	}
// 	return num;
// }
int16_t my_abs(int16_t num)
{
	return num < 0 ? -num : num;
}
void mypd(int16_t const *restrict xdata, int data_size_N, int32_t *restrict maxx, int32_t *nmax)
{
	int16_t *temp1 = malloc(sizeof(int16_t) * data_size_N);
	int16_t *output = malloc(sizeof(int16_t) * data_size_N);
	memset(output, 0, data_size_N);
	memset(temp1, 0, data_size_N);

	int max_len = 0;
	for (int i = 0; i < data_size_N - 1; i++)
	{
		temp1[i] = xdata[i + 1] - xdata[i]; /*first difference*/
	}
	printf("%d  %d \r\n", xdata[2], xdata[22]);
	for (int i = 0; i < data_size_N - 2; i++)
	{
		output[i] = (int16_t)((temp1[i + 1] - temp1[i]) * (temp1[i + 1] - temp1[i])); /*second difference and squaring*/
	}
	*nmax = 0;
	printf("%d  %d \r\n", xdata[2], xdata[22]);
	for (int i = 2; i < data_size_N - 2; i++)
	{
		if ((temp1[i] == 0) && (temp1[i - 1] == 0))
		{
			if ((temp1[i - 2] > 0) && (temp1[i + 1] < 0))
			{
				// printf("%d  %d \r\n", i, xdata[i]);
				maxx[*nmax] = i;
				(*nmax)++;
			}
		}
	}
	// for (int i = 0; i < data_size_N - 2; i++)
	// { /*determine peak value in each window*/
	// 	if (output[i] > output[max])
	// 	{
	// 		max = i;
	// 		printf("%d  %d \r\n", i, xdata[i]);
	// 	}
	// }
	free(temp1);
	free(output);
	return;
}

/********************************************
 *  Fuction : FindPV
 *  Note    :
 *******************************************/
void FindPV(int16_t *Sample, int samp_len, int32_t *Pos_Peak, int32_t *pcnt)
{
	int i = 0;
	int Vcnt = 0;
	int16_t *SampleDiff = malloc(samp_len * sizeof(int16_t));
	memset((uint8_t *)SampleDiff, 0, samp_len * sizeof(int16_t));
	*pcnt = 0;
	// step 1 :
	for (i = 0; i < samp_len - 1; i++)
	{
		if (Sample[i + 1] - Sample[i] > 0)
		{
			SampleDiff[i] = 1;
		}
		else if (Sample[i + 1] - Sample[i] < 0)
		{
			SampleDiff[i] = -1;
		}
		else
		{
			SampleDiff[i] = 0;
		}
	}

	// step 2 :
	for (i = 0; i < samp_len - 1; i++)
	{
		if (SampleDiff[i] == 0)
		{
			if (i == (samp_len - 2))
			{
				if (SampleDiff[i - 1] >= 0)
				{
					SampleDiff[i] = 1;
				}
				else
				{
					SampleDiff[i] = -1;
				}
			}
			else
			{
				if (SampleDiff[i + 1] >= 0)
				{
					SampleDiff[i] = 1;
				}
				else
				{
					SampleDiff[i] = -1;
				}
			}
		}
	}
	printf("\n");
	// for (int i = 0; i < SAMPLE_MAX; i++)
	// {
	//     printf("diff2[%d] = %d \t", i, SampleDiff[i]);
	//     if ( ((i + 1) & 0x03) == 0)
	//     {
	//         printf("\n");
	//     }
	// }
	// step 3 :
	int last_index = 0;
	for (i = 0; i < samp_len - 1; i++)
	{
		if (SampleDiff[i + 1] - SampleDiff[i] == -2) //
		{
			// if (Sample[i + 1] - Sample[i] > 0)
			{
				Pos_Peak[*pcnt] = i + 1;
				(*pcnt)++;
			}
		}
		else if (SampleDiff[i + 1] - SampleDiff[i] == 2) //
		{
			Vcnt++;
		}
	}
}

#if 0
//模糊熵
float FuzzyEntropy( int16_t* data_test,int16_t N,int16_t m,int16_t dim)
{
	float s = 0;
	for (int i = 0; i < N; i++) {
		s += data_test[i];
	}
	float data_mean = s / N;
	float var_data = 0;
	for (int i = 0; i < N; i++) {
		var_data += pow((data_test[i] - data_mean), 2);
	}
	float data_std = sqrt(var_data / N);
	float r = 0.15 * data_std;
	float R[2];
	for (int v = m; v < m + 2; v++) {
		float** a;
		float** b;
		int i, j;
		a = (float**)fh_malloc(sizeof(float*) * N);//为二维数组分配3行 
		b = (float**)fh_malloc(sizeof(float*) * N);
		for (i = 0; i < N; i++) {//为每列分配4个大小空间 
			a[i] = (float*)fh_malloc(sizeof(float)*4);
			b[i] = (float*)fh_malloc(sizeof(float)*4);
		}
		//初始化 
		for (i = 0; i < N - v + 1; i++) {
			for (j = 0; j < v; j++) {
				a[i][j] = data_test[i + j];
				b[i][j] = 0;
			}
		}
		for (int i = 0; i < N - v + 1; i++) {
			float sum = 0;
			for (int j = 0; j < v; j++) {
				sum += a[i][j];
			}
			for (int z = 0; z < v; z++) {
				float s_1 = sum / v;
				float s = a[i][z];
				b[i][z] = s - s_1;
			}
		}
		float* K;              //数据长度
		K = (float*)fh_malloc(sizeof(float*) * (N - v + 1));
		for (int i = 0; i < N - v + 1; i++) {
			float* cha;
			cha = (float*)fh_malloc(sizeof(float*) * (N - v + 1));
			for (int j = 0; j < N - v + 1; j++) {
				float M_1 = 0, M = 0;
				for (int z = 0; z < v; z++) {
					M_1 = fabsf(b[i][z] - b[j][z]);
					M = M > M_1 ? M : M_1;
					cha[j] = fh_exp(-fh_pow(M, dim) / r);
				}
			}
			float result_cha = 0;
			for (int i = 0; i < N - v + 1; i++) {
				result_cha += cha[i];
			}
			K[i] = (result_cha - 1) / (N - v);
			free(cha);
		}
		for (i = 0; i < N; i++) {//为每列分配4个大小空间 
			free(a[i]);
			free(b[i]);
		}
		free(b);
		free(a);
		float result_1 = 0;
		for (int i = 0; i < N - v + 1; i++) {
			result_1 += K[i];
		}
		R[v - m] = result_1 / (N - v + 1);
		free(K);
	}
	float FE_test = log(R[0]) - log(R[1]);
	return FE_test;
}
#endif

#if 1

/*
 * 函数:  findPeaks
 * 参数:  *src        源数据数组
 *          src_lenth   源数据数组长度
 *          distance    峰与峰,谷与谷的搜索间距
 *          *indMax     找到的峰的index数组
 *          *indMax_len 数组长度
 *          *indMin     找到的谷的index数组
 *          *indMin_len 数组长度
 */

void findPeaks(int16_t *src, int32_t src_lenth, int16_t distance, int16_t *indMax, int32_t *indMax_len, int16_t *indMin, int32_t *indMin_len)
{
	int16_t *sign = (int16_t *)malloc(src_lenth * sizeof(int16_t));
	int32_t max_index = 0,
			min_index = 0;
	*indMax_len = 0;
	*indMin_len = 0;
	int32_t i, j, k;

	log_info("goto  %d", __LINE__);

	for (i = 1; i < src_lenth; i++)
	{
		int16_t diff = src[i] - src[i - 1];
		if (diff > 0)
			sign[i - 1] = 1;
		else if (diff < 0)
			sign[i - 1] = -1;
		else
			sign[i - 1] = 0;
	}
	for (j = 1; j < src_lenth - 1; j++)
	{
		int16_t diff = sign[j] - sign[j - 1];
		if (diff < 0)
			indMax[max_index++] = j;
		else if (diff > 0)
			indMin[min_index++] = j;
	}
	log_info("goto  %d", __LINE__);
	log_info("max index: %d min index: %d", max_index, min_index);
	int16_t *flag_max_index = (int16_t *)malloc(sizeof(int16_t) * (max_index > min_index ? max_index : min_index));
	int16_t *idelete = (int16_t *)malloc(sizeof(int16_t) * (max_index > min_index ? max_index : min_index));
	int16_t *temp_max_index = (int16_t *)malloc(sizeof(int16_t) * (max_index > min_index ? max_index : min_index));
	int32_t bigger = 0;
	int16_t tempvalue = 0;
	// 波峰
	for (int16_t i = 0; i < max_index; i++)
	{
		flag_max_index[i] = 0;
		idelete[i] = 0;
	}
	log_info("goto  %d", __LINE__);

	// for (i = 0; i < max_index; i++)
	{
		tempvalue = -1;
		for (j = 0; j < max_index; j++)
		{
			if (!flag_max_index[j])
			{
				if (src[indMax[j]] > tempvalue)
				{
					bigger = j;
					tempvalue = src[indMax[j]];
				}
			}
		}
		flag_max_index[bigger] = 1;
		if (!idelete[bigger])
		{
			for (k = 0; k < max_index; k++)
			{
				idelete[k] |= (indMax[k] - distance <= indMax[bigger] & indMax[bigger] <= indMax[k] + distance);
			}
			idelete[bigger] = 0;
		}
	}
	log_info("goto  %d", __LINE__);

	for (i = 0, j = 0; i < max_index; i++)
	{
		if (!idelete[i])
			temp_max_index[j++] = indMax[i];
	}
	for (i = 0; i < max_index; i++)
	{
		if (i < j)
			indMax[i] = temp_max_index[i];
		else
			indMax[i] = 0;
	}
	max_index = j;
	log_info("goto  %d", __LINE__);

	// 波谷
	for (int16_t i = 0; i < min_index; i++)
	{
		flag_max_index[i] = 0;
		idelete[i] = 0;
	}
	for (i = 0; i < min_index; i++)
	{
		tempvalue = 1;
		for (j = 0; j < min_index; j++)
		{
			if (!flag_max_index[j])
			{
				if (src[indMin[j]] < tempvalue)
				{
					bigger = j;
					tempvalue = src[indMin[j]];
				}
			}
		}
		flag_max_index[bigger] = 1;
		if (!idelete[bigger])
		{
			for (k = 0; k < min_index; k++)
			{
				idelete[k] |= (indMin[k] - distance <= indMin[bigger] & indMin[bigger] <= indMin[k] + distance);
			}
			idelete[bigger] = 0;
		}
	}
	log_info("goto  %d", __LINE__);

	for (i = 0, j = 0; i < min_index; i++)
	{
		if (!idelete[i])
			temp_max_index[j++] = indMin[i];
	}
	for (i = 0; i < min_index; i++)
	{
		if (i < j)
			indMin[i] = temp_max_index[i];
		else
			indMin[i] = 0;
	}
	min_index = j;
	log_info("goto  %d", __LINE__);

	*indMax_len = max_index;
	*indMin_len = min_index;

	free(sign);
	free(flag_max_index);
	free(temp_max_index);
	free(idelete);
}

#endif
#if (USE_IIR_FILTER_TYPE == 0)
#if 1

/* Kalman Structure Initialization */
// kalman_filter_data_s kalman_data =
// {
// 	/* Transition matrix: 2x2 */
// 	/* float Phi_matrix[4]; */
// 	{1.0, 0.25e-3, 0.0, 1.0},
// 	/* Q covariance plant noise matrix: 2x2 */
// 	/* float Q_matrix[4]; */
// 	{0.0, 0.0, 0.0, 1082.323},
// 	/* Sensitivity matrix: 1X2 */
// 	/* float H_matrix[2]; */
// 	{1.0, 0.0},
// 	/* Observation noise: R covariance matrix 1x1 */
// 	/* float R_matrix; */
// 	0.04,
// 	/* P plus current covariance matrix 2x2: estimate error */
// 	/* float P_plus[4]; */
// 	{0.04, 160.0, 160.0, 641082.323},
// 	/* x plus current state vector 2x1: value, speed */
// 	/* float x_plus[2]; */
// 	{0.0, 0.0},
// };

// 	/*
//  Scalar Kalman Filter
//  Input
//  float input - input measured signal
//  kalman_filter_data_s* kalman_data - Kalman filter data
//  Return
//  float x[0] estimate value
// */
// float KalmanFilterExample(float input, kalman_filter_data_s* kalman_data)
// {
// 	float P_minus[4]; /* matrix 2x2 */
// 	float x_minus[2]; /* vector 2x1 */
// 	float K_gain[2];  /* matrix 2x1 */
// 	float temp_help;

// 	/* Prediction Step */
// 	x_minus[0] = kalman_data->Phi_matrix[0]*kalman_data->x_plus[0] + kalman_data->Phi_matrix[1]*kalman_data->x_plus[1];
// 	x_minus[1] = kalman_data->Phi_matrix[2]*kalman_data->x_plus[0] + kalman_data->Phi_matrix[3]*kalman_data->x_plus[1];
// 	P_minus[0] = (kalman_data->Phi_matrix[0]*kalman_data->P_plus[0] + kalman_data->Phi_matrix[1]*kalman_data->P_plus[2])*kalman_data->Phi_matrix[0];
// 	P_minus[0] += (kalman_data->Phi_matrix[0]*kalman_data->P_plus[1] + kalman_data->Phi_matrix[1]*kalman_data->P_plus[3])*kalman_data->Phi_matrix[1];
// 	P_minus[0] += kalman_data->Q_matrix[0];
// 	P_minus[1] = (kalman_data->Phi_matrix[0]*kalman_data->P_plus[0] + kalman_data->Phi_matrix[1]*kalman_data->P_plus[2])*kalman_data->Phi_matrix[2];
// 	P_minus[1] += (kalman_data->Phi_matrix[0]*kalman_data->P_plus[1] + kalman_data->Phi_matrix[1]*kalman_data->P_plus[3])*kalman_data->Phi_matrix[3];
// 	P_minus[1] += kalman_data->Q_matrix[1];
// 	P_minus[2] = (kalman_data->Phi_matrix[2]*kalman_data->P_plus[0] + kalman_data->Phi_matrix[3]*kalman_data->P_plus[2])*kalman_data->Phi_matrix[0];
// 	P_minus[2] += (kalman_data->Phi_matrix[2]*kalman_data->P_plus[1] + kalman_data->Phi_matrix[3]*kalman_data->P_plus[3])*kalman_data->Phi_matrix[1];
// 	P_minus[2] += kalman_data->Q_matrix[2];
// 	P_minus[3] = (kalman_data->Phi_matrix[2]*kalman_data->P_plus[0] + kalman_data->Phi_matrix[3]*kalman_data->P_plus[2])*kalman_data->Phi_matrix[2];
// 	P_minus[3] += (kalman_data->Phi_matrix[2]*kalman_data->P_plus[1] + kalman_data->Phi_matrix[3]*kalman_data->P_plus[3])*kalman_data->Phi_matrix[3];
// 	P_minus[3] += kalman_data->Q_matrix[3];
// 	/* Kalman Gain */
// 	temp_help = (kalman_data->H_matrix[0]*P_minus[0] + kalman_data->H_matrix[1]*P_minus[2])*kalman_data->H_matrix[0];
// 	temp_help += (kalman_data->H_matrix[0]*P_minus[1] + kalman_data->H_matrix[1]*P_minus[3])*kalman_data->H_matrix[1];
// 	temp_help += kalman_data->R_matrix;
// 	K_gain[0] = (kalman_data->H_matrix[0]*P_minus[0] + kalman_data->H_matrix[1]*P_minus[1])/temp_help; /* temp_help shall be !=0 */
// 	K_gain[1] = (kalman_data->H_matrix[0]*P_minus[2] + kalman_data->H_matrix[1]*P_minus[3])/temp_help;
// 	/* Correction Step */
// 	kalman_data->P_plus[0] = (1.0 - K_gain[0]*kalman_data->H_matrix[0])*P_minus[0] - K_gain[0]*kalman_data->H_matrix[1]*P_minus[2];
// 	kalman_data->P_plus[1] = (1.0 - K_gain[0]*kalman_data->H_matrix[0])*P_minus[1] - K_gain[0]*kalman_data->H_matrix[1]*P_minus[3];
// 	kalman_data->P_plus[2] = -K_gain[1]*kalman_data->H_matrix[0]*P_minus[0] + (1.0 - K_gain[1]*kalman_data->H_matrix[1])*P_minus[2];
// 	kalman_data->P_plus[3] = -K_gain[1]*kalman_data->H_matrix[0]*P_minus[1] + (1.0 - K_gain[1]*kalman_data->H_matrix[1])*P_minus[3];
// 	kalman_data->x_plus[0] = x_minus[0] + K_gain[0]*(input - x_minus[0]);
// 	kalman_data->x_plus[1] = x_minus[1] + K_gain[1]*(input - x_minus[0]);

// 	return kalman_data->x_plus[0];
// }

#endif
#if 0
typedef  struct{
	double filterValue;  //k-1时刻的滤波值，即是k-1时刻的值
	double kalmanGain;   //   Kalamn增益
	double A;   // x(n)=A*x(n-1)+u(n),u(n)~N(0,Q)
	double H;   // z(n)=H*x(n)+w(n),w(n)~N(0,R)
	double Q;   //预测过程噪声偏差的方差
	double R;   //测量噪声偏差，(系统搭建好以后，通过测量统计实验获得)
	double P;   //估计误差协方差
}  KalmanInfo;
KalmanInfo k_info;

/**
* @brief Init_KalmanInfo   初始化滤波器的初始值
* @param info  滤波器指针
* @param Q 预测噪声方差 由系统外部测定给定
* @param R 测量噪声方差 由系统外部测定给定
*/
void Init_KalmanInfo(KalmanInfo* info, double Q, double R)
{
	info->A = 1;  //标量卡尔曼
	info->H = 1;  //
	info->P = 0;  //后验状态估计值误差的方差的初始值（不要为0问题不大）
	info->Q = Q;    //预测（过程）噪声方差 影响收敛速率，可以根据实际需求给出
	info->R = R;    //测量（观测）噪声方差 可以通过实验手段获得
	info->filterValue = 0;// 测量的初始值
}
double KalmanFilter(KalmanInfo* kalmanInfo, double lastMeasurement)
{
	//预测下一时刻的值
	double predictValue = kalmanInfo->A* kalmanInfo->filterValue;   //x的先验估计由上一个时间点的后验估计值和输入信息给出，此处需要根据基站高度做一个修改
	
	//求协方差
	kalmanInfo->P = kalmanInfo->A*kalmanInfo->A*kalmanInfo->P + kalmanInfo->Q;  //计算先验均方差 p(n|n-1)=A^2*p(n-1|n-1)+q
	double preValue = kalmanInfo->filterValue;  //记录上次实际坐标的值
 
	//计算kalman增益
	kalmanInfo->kalmanGain = kalmanInfo->P*kalmanInfo->H / (kalmanInfo->P*kalmanInfo->H*kalmanInfo->H + kalmanInfo->R);  //Kg(k)= P(k|k-1) H’ / (H P(k|k-1) H’ + R)
	//修正结果，即计算滤波值
	kalmanInfo->filterValue = predictValue + (lastMeasurement - predictValue)*kalmanInfo->kalmanGain;  //利用残余的信息改善对x(t)的估计，给出后验估计，这个值也就是输出  X(k|k)= X(k|k-1)+Kg(k) (Z(k)-H X(k|k-1))
	//更新后验估计
	kalmanInfo->P = (1 - kalmanInfo->kalmanGain*kalmanInfo->H)*kalmanInfo->P;//计算后验均方差  P[n|n]=(1-K[n]*H)*P[n|n-1]
 
	return  kalmanInfo->filterValue;
}
#endif
#if 0
typedef struct
{
	float Last_P; // 上次估算协方差 不可以为0 ! ! ! ! !
	float Now_P;  // 当前估算协方差
	float out;	  // 卡尔曼滤波器输出
	float Kg;	  // 卡尔曼增益
	float Q;	  // 过程噪声协方差
	float R;	  // 观测噪声协方差
} KalmanInfo;
KalmanInfo k_info;

void Kalman_Init()
{
	k_info.Last_P = 100;
	k_info.Now_P = 0;
	k_info.out = 0;
	k_info.Kg = 1;
	k_info.Q = 0.1;
	k_info.R = 200;
}

/**
 *卡尔曼滤波器
 *@param 	Kalman *kfp 卡尔曼结构体参数
 *   			float input 需要滤波的参数的测量值（即传感器的采集值）
 *@return 滤波后的参数（最优值）
 */
float KalmanFilter(KalmanInfo *kfp, float input)
{
	// 预测协方差方程：k时刻系统估算协方差 = k-1时刻的系统协方差 + 过程噪声协方差
	kfp->Now_P = kfp->Last_P + kfp->Q;
	// 卡尔曼增益方程：卡尔曼增益 = k时刻系统估算协方差 / （k时刻系统估算协方差 + 观测噪声协方差）
	kfp->Kg = kfp->Now_P / (kfp->Now_P + kfp->R);
	// 更新最优值方程：k时刻状态变量的最优值 = 状态变量的预测值 + 卡尔曼增益 * （测量值 - 状态变量的预测值）
	kfp->out = kfp->out + kfp->Kg * (input - kfp->out); // 因为这一次的预测值就是上一次的输出值
	// 更新协方差方程: 本次的系统协方差付给 kfp->LastP 威下一次运算准备。
	kfp->Last_P = (1 - kfp->Kg) * kfp->Now_P;
	return kfp->out;
}
#endif
// // IIR low pass filter to tracking DC level
// // y(n) = (x(n) - y(n-1))/(2^K) + y(n-1)
// // To ensure there is sufficient accuracy,
// // the input is left shift by K
#define DC_K 13
int32_t dc_estimator(int32_t *p, uint16_t x)
{
	*p += ((((int32_t)x << DC_K) - *p) >> DC_K);
	return (*p >> DC_K);
}
//  y += b * (x - y)
// int16_t dc_estimator(register int32_t *p, register int16_t x)
// {
// 	/* Noise shaped DC estimator. */
// 	*p += ((((int32_t)x << 16) - *p) >> 11);
// 	//*p += ((((int32_t)x << 16) - *p) >> 11);

// 	return (*p >> 16);
// }

#endif
int main(void)
{
	FILE *pFile;
	char fw_buf[16];
	struct timeval begin, end;
	long microseconds = 0;
	srand((unsigned)time(NULL)); // 初始化随机数
	// Start measuring time

	gettimeofday(&begin, 0);
	gettimeofday(&end, 0);
	microseconds = end.tv_usec - begin.tv_usec;
	printf("elapsed %ld us\r\n", microseconds);
	// 获取源数据
	pFile = fopen("raw_data.txt", "r"); // test_con  raw_data
	if (pFile == NULL)
	{
		printf(" ERROR\r\n");
	}
	fseek(pFile, SEEK_SET, SEEK_END); // 定位到文件尾，偏移量为0

	int file_len = ftell(pFile); // 返回当前定位的文件位置（也就是文件总长度）
	char *data_str = malloc(file_len + 1);
	fseek(pFile, SEEK_SET, SEEK_SET); // 定位到文件首，偏移量为0
	int data_len = (int)(file_len / 2) + 1;
	int *data_buf = malloc(data_len * sizeof(int));
	fread(data_str, file_len + 1, 1, pFile);
	fclose(pFile);
	int len;
	int *raw_data = data_buf;

	istrToarray(data_str, data_buf, &len);
	int *filter_data = malloc(len * sizeof(int));
	int *fp_data = filter_data;
	FILTER sig_filter;
	sig_filter.offset = 0;
	int fir_data = 0;
#if (USE_IIR_FILTER_TYPE == 0)
	int lp_data = 0;

	int32_t ir_2nd_dc_register = 0;
	// Kalman_Init();
	//        Init_KalmanInfo(&k_info, 100, 200);
#endif
	FilterTypeDef avg_filter;
	Moving_Average_Init(&avg_filter);
	gettimeofday(&begin, 0);
	smoth(raw_data, len, fp_data);

	// 	for (int i = 0; i < len; i++)
	// 	{
	// 		//	fir_data = Moving_Average_Compute(raw_data[i],&avg_filter);
	// 		//fir_data = fir_filter(raw_data[i], &sig_filter);		fir_data = (int)KalmanFilterExample(fir_data, &kalman_data);
	// 		// fir_data = KalmanFilter(&k_info, raw_data[i]);

	// 		 fir_data =	(int)KalmanFilterExample(raw_data[i], &kalman_data);  fir_data = fir_filter(fir_data, &sig_filter);

	// #if (USE_IIR_FILTER_TYPE == 0)

	// 		lp_data = 0;
	// // lp_data = dc_estimator(&ir_2nd_dc_register, fir_data);
	// #else
	// 		lp_data = 0;
	// #endif
	// 		if ((fir_data - lp_data) < 0)
	// 		{
	// 			*fp_data++ = 0;
	// 		}
	// 		else
	// 		{
	// 			*fp_data++ = fir_data - lp_data;
	// 		}
	// 	}
	gettimeofday(&end, 0);
	microseconds = end.tv_usec - begin.tv_usec;
	printf("filter %d-elapsed %ld us\r\n", __LINE__, microseconds);

	log_info("write filtered data to file");
	pFile = fopen("filter_data.txt", "wb");
	fp_data = filter_data;
	for (size_t i = 0; i < len; i++)
	{
		/* code */
		sprintf(fw_buf, "%d \0", *fp_data++);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);

#if (EXTREMA_TEST)
	int16_t *p_xdata = malloc(len * sizeof(int16_t));

	for (int i = 0; i < len; i++)
	{
		p_xdata[i] = (int16_t)(raw_data[i]);
	}
	int32_t *pmaxx_data = malloc(len * sizeof(int32_t));
	int16_t *pmaxy_data = malloc(len * sizeof(int16_t));
	int32_t nmax;
	gettimeofday(&begin, 0);

	my_find_extrema(p_xdata, len, pmaxx_data, pmaxy_data, &nmax);
	gettimeofday(&end, 0);
	microseconds = end.tv_usec - begin.tv_usec;
	printf("%d-elapsed %ld us\r\n", __LINE__, microseconds);

	log_debug("my_find_extrema n max is %d ", nmax);
	log_info("write extrema data");
	pFile = fopen("m_xdata.txt", "wb");
	for (int i = 0; i < nmax; i++)
	{
		sprintf(fw_buf, "%d \0", pmaxx_data[i]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);
	pFile = fopen("m_ydata.txt", "wb");
	for (int i = 0; i < nmax; i++)
	{
		sprintf(fw_buf, "%d \0", pmaxy_data[i]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);
	free(pmaxx_data);
	free(pmaxy_data);

#endif

#if (1)
	// 获取源数据

	int16_t *p_xdata_f = malloc(len * sizeof(int16_t));
	for (int i = 0; i < len; i++)
	{
		p_xdata_f[i] = (int16_t)(filter_data[i]);
	}

	int32_t *pmaxx_data_f = malloc(len * sizeof(int32_t));
	int16_t *pmaxy_data_f = malloc(len * sizeof(int16_t));
	nmax = 0;
	gettimeofday(&begin, 0);

	my_find_extrema(p_xdata_f, len, pmaxx_data_f, pmaxy_data_f, &nmax);
	gettimeofday(&end, 0);
	microseconds = end.tv_usec - begin.tv_usec;
	printf("%d-elapsed %ld us\r\n", __LINE__, microseconds);

	log_debug("my_find_extrema filtered n max is %d ", nmax);
	pFile = fopen("m_xdata_f.txt", "wb");
	for (int i = 0; i < nmax; i++)
	{
		sprintf(fw_buf, "%d \0", pmaxx_data_f[i]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);
	pFile = fopen("m_ydata_f.txt", "wb");
	for (int i = 0; i < nmax; i++)
	{
		sprintf(fw_buf, "%d \0", pmaxy_data_f[i]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);
	free(pmaxx_data_f);
	free(pmaxy_data_f);
#endif
#if 0
	int16_t *peakFs = malloc(len * sizeof(int16_t));
	int16_t *peakFs2 = malloc(len * sizeof(int16_t));

	int32_t peakFs_len = 0;
	int32_t peakFs2_len = 0;
	uint16_t peak_count = 0;
	uint16_t peak_len = 0;
	log_info("goto findPeaks");
	findPeaks(p_xdata, len, 3, peakFs, &peakFs_len, peakFs2, &peakFs2_len);
	log_debug("findPeaks n max is %d ", peakFs_len);

	pFile = fopen("m2_xdata.txt", "wb");
	for (int i = 0; i < peakFs_len; i++)
	{
		sprintf(fw_buf, "%d \0", peakFs[i]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);
	pFile = fopen("m2_ydata.txt", "wb");
	for (int i = 0; i < peakFs_len; i++)
	{
		sprintf(fw_buf, "%d \0", p_xdata[peakFs[i]]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);
	free(peakFs);
	free(peakFs2);
#endif
#if 0
	int32_t *px_index_data = malloc(len * sizeof(int32_t));
	int32_t len_max = 0;
	// mypd(p_xdata, len,px_index_data, &len_max);
	FindPV(p_xdata, len, px_index_data, &len_max);

	log_debug("mypd n max is %d ", len_max);
	log_info("write mypd data");
	pFile = fopen("pd_xdata.txt", "wb");
	for (int i = 0; i < len_max; i++)
	{
		sprintf(fw_buf, "%d \0", px_index_data[i]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);

	pFile = fopen("pd_ydata.txt", "wb");
	for (int i = 0; i < len_max; i++)
	{
		sprintf(fw_buf, "%d \0", p_xdata[px_index_data[i]]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);

#endif
#if 0
	int *p_index = malloc(len * sizeof(int));
	memset(p_index, 0, len * sizeof(int));
	int len_index = 0;
	log_info("befor ampd");
	my_AMPD(p_xdata, len, p_index, &len_index);
	log_info("my_AMPD len index is %d", len_index);
	log_info("write extrema data");
	pFile = fopen("x_ampd.txt", "wb");
	for (size_t i = 0; i < len_index; i++)
	{
		sprintf(fw_buf, "%d \0", p_index[i]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);
	pFile = fopen("y_ampd.txt", "wb");

	for (size_t i = 0; i < len_index; i++)
	{
		/* code */

		sprintf(fw_buf, "%d \0", p_xdata[p_index[i]]);
		fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	}
	fclose(pFile);
	free(p_index);
#endif
	//   pFile = fopen_s("x_data.txt", "wb");
	//    for (size_t i = 0; i < nmax; i++)
	//   {
	//       /* code */
	//	 sprintf_s(fw_buf, "%d \0", pmaxx_data[i]);
	//       fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	//   }
	//   fclose(pFile);
	// pFile = fopen_s("y_data.txt", "wb");
	//	for (size_t i = 0; i < nmax; i++)
	//{
	//	/* code */
	//	sprintf_s(fw_buf, "%d \0", pmaxy_data[i]);
	//	fwrite(fw_buf, strlen(fw_buf), 1, pFile);
	//}
	// fclose(pFile);

	// THE FOLLOWING CODE WILL PLOT THE TIME DOMAIN SIGNAL AFTER BASELINE INTERFERENCE REMOVAL
	// FILE *gnuplot = popen("gnuplot -persistent","w");
	// fprintf(gnuplot,"plot'-'\n");

	// for (int j=0;j<=len;j++)
	// {
	// 	fprintf(gnuplot,"%g \n", p_xdata[j]);
	// }

	// fprintf(gnuplot,"e\n");
	// fflush(gnuplot);

	free(filter_data);
	free(data_buf);
	free(data_str);

	return 0;
}