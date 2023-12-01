
/* Includes ------------------------------------------------------------------*/
#include "moving_average.h"
#include <stdio.h>
/**
  * @brief  This function initializes filter's data structure.
	* @param  filter_struct : Data structure
  * @retval None.
  */
void Moving_Average_Init(FilterTypeDef* filter_struct)
{
	filter_struct->Sum = 0;
	filter_struct->WindowPointer = 0;
	
	for(uint32_t i=0; i<WindowLength; i++)
	{
		filter_struct->History[i] = 0;
	}
}

/**
  * @brief  This function filters data with moving average filter.
	* @param  raw_data : input raw sensor data.
	* @param  filter_struct : Data structure
  * @retval Filtered value.
  */
uint16_t Moving_Average_Compute(uint16_t raw_data, FilterTypeDef* filter_struct)
{
	filter_struct->Sum += raw_data;
	filter_struct->Sum -= filter_struct->History[filter_struct->WindowPointer];
	filter_struct->History[filter_struct->WindowPointer] = raw_data;
	if(filter_struct->WindowPointer < WindowLength - 1)
	{
		filter_struct->WindowPointer += 1;
	}
	else
	{
		filter_struct->WindowPointer = 0;
	}
	return filter_struct->Sum/WindowLength;
}

/**
 * 函 数 名：smoth(x,n,y)
 * 功能描述：五点三次平滑
 * 输入参数：x（输入数据），n（系数个数）、y（平滑后序列）
 * 返 回 值：整型数字。计算成功则返回1，否则返回0
 */
int smoth(int *x,int n,int *y)
{
    int i;

    if(n < 5)
    {
        printf("n should be at least 5.\n");
        return 0;
    }

    y[0] = (69*x[0]+4*(x[1]+x[3])-6*x[2]-x[4])/70;        /* 前两点*/
    y[1] = (2*x[0]+27*x[1]+12*x[2]-8*x[3]+2*x[4])/35;

    for(i=2; i<n-2; i++)                                          /* 中间点*/
    {
        y[i] = (-3*x[i-2]+12*x[i-1]+17*x[i]+12*x[i+1]-3*x[i+2])/35;
    }

    y[n-2] = (2*x[n-5]-8*x[n-4]+12*x[n-3]+27*x[n-2]+2*x[n-1])/35;
    y[n-1] = (x[n-5]+4*x[n-4]-6*x[n-3]+4*x[n-2]+69*x[n-1])/70;
    return 1;
}

// /**
//  * 函 数 名：smoth(x,n,y)
//  * 功能描述：五点三次平滑
//  * 输入参数：x（输入数据），n（系数个数）、y（平滑后序列）
//  * 返 回 值：整型数字。计算成功则返回1，否则返回0
//  */
// int smoth(double *x,int n,double *y)
// {
//     int i;

//     if(n < 5)
//     {
//         printf("n should be at least 5.\n");
//         return 0;
//     }

//     y[0] = (69.0*x[0]+4.0*(x[1]+x[3])-6.0*x[2]-x[4])/70.0;        /* 前两点*/
//     y[1] = (2.0*x[0]+27.0*x[1]+12.0*x[2]-8.0*x[3]+2.0*x[4])/35.0;

//     for(i=2; i<n-2; i++)                                          /* 中间点*/
//     {
//         y[i] = (-3.0*x[i-2]+12.0*x[i-1]+17.0*x[i]+12.0*x[i+1]-3.0*x[i+2])/35.0;
//     }

//     y[n-2] = (2.0*x[n-5]-8.0*x[n-4]+12.0*x[n-3]+27.0*x[n-2]+2.0*x[n-1])/35.0;
//     y[n-1] = (x[n-5]+4.0*x[n-4]-6.0*x[n-3]+4.0*x[n-2]+69.0*x[n-1])/70.0;
//     return 1;
// }
