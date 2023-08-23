#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

  int istrToarray(char *str_buf, int *data, int *dlen)
{
    char *token;
    char tmpbuffer[6];
    int num = 0;
    token = strtok(str_buf, " ");
    while (token != NULL)
    {
        *data++ = atoi(token);
        num++;
        token = strtok(NULL, " ");
    }
    printf("num is %d \r\n", num);
    *dlen = num;
    return 0;
}