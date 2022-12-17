#include "stdlib.h"
#include "stdio.h"
#include "time.h"

void log_warn(const char* logging) {
    time_t now = time(NULL);
    struct tm *tm_date = localtime(&now);
    int hour = tm_date->tm_hour;
    printf("\033[0;33m");
    printf("%s", logging);
}