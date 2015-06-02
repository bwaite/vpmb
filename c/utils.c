#include "utils.h"

double min(double i, double j)
{
    if(i < j) {
        return i;
    }
    return j;
}

double max(double i, double j)
{
    if(i > j) {
        return i;
    }
    return j;
}

void lowercase_string(char *str)
{
    size_t i;
    for(i = 0; i < strlen(str); i++) {
        str[i] = tolower(str[i]);
    }
}
