
#include "my_fir.h"
static const int16_t coeffs[12] = {
    210, 254, 385, 592,
    861, 1170, 1498, 1819,
    2111, 2351, 2522, 2611};

int32_t mul16(register int16_t x, register int16_t y)
{

    return x * y;
}
int64_t mul32(register int32_t x, register int32_t y)
{

    return x * y;
}

int16_t fir_filter(int16_t input, FILTER *param)
{
    int i;
    int32_t z;
    int16_t *buf = param->buf;

    buf[param->offset] = input;
    z = mul16(coeffs[11], buf[(param->offset - 11) & 0x1F]);
    for (i = 0; i < 11; i++)
        z += mul16(coeffs[i], buf[(param->offset - i) & 0x1F] + buf[(param->offset - 22 + i) & 0x1F]);
    param->offset = (param->offset + 1) & 0x1F;

    return z >> 15;
}
 