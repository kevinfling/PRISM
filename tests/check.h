#ifndef PRISM_TEST_CHECK_H
#define PRISM_TEST_CHECK_H
#include <stdio.h>
#include <stdlib.h>
#define CHECK(c)                                                                                   \
    do {                                                                                           \
        if (!(c)) {                                                                                \
            fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #c);                                \
            exit(1);                                                                               \
        }                                                                                          \
    } while (0)
#endif
