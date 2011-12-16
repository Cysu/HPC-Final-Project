#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<lin_kernighan.h>

void lin_kernighan(int size, int* x_of_cities, int* y_of_cities, int* opt_tour, double* opt_tour_len) {
	gen_rand_tour(size, opt_tour);
	*opt_tour_len = tour_len(size, opt_tour, x_of_cities, y_of_cities);
}
