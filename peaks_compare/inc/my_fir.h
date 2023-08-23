#ifndef _MY_FIR_H
#define _MY_FIR_H

#include <stdio.h>
#include <stdint.h>

typedef struct
{
    int16_t offset;
    int16_t buf[32];
} FILTER;
int16_t fir_filter(int16_t input, FILTER *param);

#endif