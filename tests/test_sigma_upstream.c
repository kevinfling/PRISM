/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling
 * Reproduce pinned upstream sequence policies independently of PRISM's adapter.
 * These checks explain why the adapter deliberately uses scalar cached draws. */
#include "check.h"
#include <sigma/rng.h>
#include <string.h>
int main(void) {
    sigma_state_t scalar, buffered;
    sigma_rng_init(&scalar, 42);
    scalar.vtable = &sigma_vtable_xoshiro256ss;
    buffered = scalar;
    double expected[3], actual[3];
    for (size_t i = 0; i < 3; ++i)
        expected[i] = sigma_rng_next_normal(&scalar);
    actual[0] = sigma_rng_next_normal(&buffered);
    sigma_rng_fill_normal(&buffered, actual + 1, 2);
    CHECK(actual[0] == expected[0]);
    CHECK(actual[1] != expected[1] || actual[2] != expected[2]);
    printf("pinned upstream cached normal scalar=[%.17g,%.17g,%.17g] mixed=[%.17g,%.17g,%.17g]\n",
           expected[0], expected[1], expected[2], actual[0], actual[1], actual[2]);
    sigma_rng_init(&scalar, 42);
    sigma_rng_init(&buffered, 42);
    scalar.vtable = &sigma_vtable_xoshiro256ss;
    double a[8], b[8];
    sigma_rng_fill_d(&scalar, a, 8);
    sigma_rng_fill_d(&buffered, b, 8);
    int simd_selected = buffered.vtable != &sigma_vtable_xoshiro256ss;
    int equal = memcmp(a, b, sizeof(a)) == 0;
    CHECK(simd_selected ? !equal : equal);
    printf("pinned upstream selected_uniform_simd=%d scalar_buffer_equal=%d\n", simd_selected,
           equal);
    if (simd_selected) {
        sigma_state_t all, chunks;
        sigma_rng_init(&all, 42);
        chunks = all;
        sigma_rng_fill_d(&all, a, 8);
        sigma_rng_fill_d(&chunks, b, 4);
        sigma_rng_fill_d(&chunks, b + 4, 4);
        CHECK(memcmp(a, b, sizeof(a)) != 0);
        puts("pinned upstream selected SIMD uniform chunks produce a different sequence");
    }
    return 0;
}
