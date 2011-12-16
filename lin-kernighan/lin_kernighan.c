#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>
#include "lin_kernighan.h"
#include "my_conf.h"

void lin_kernighan(int size, double* x_of_cities, double* y_of_cities, int* opt_tour, double* opt_tour_len) {
	/*
	 * declare and init variable
	 */
	int* x1;
	int* x2;
	int* y1;
	int* y2;
	int* T;
	int* y12;
	int* y22;
	int **near;
	int xs, y11, y21;
	int G_opt, G;
	int i, j, k;
	int btflag = 1;

	T = (int*)malloc(size * sizeof(int));
	/*
	 * 1 generate a random tour
	 */
	gen_rand_tour(size, opt_tour);

	/* 2-6
	 * x1 = (t1, t2)
	 * t1 = T[size - 1]
	 * t2 = T[0]
	 */
	memcpy(T, opt_tour, size * sizeof(int));
	G_opt = 0;
	G = 0;
	for (i = 0; i < size; i++) {
		for (xs = 0; xs < 2; xs++) {
			y11 = T[0];
			y12 = gen_inc_array(y11, x_of_cities, y_of_cities);
			for (j = 1; j < size; j++) {
				/*
				 * y1 = (y11, y12[j])
				 */
				if (y12[j] != T[1] && y12[j] != T[size - 1] &&
					dist(y11, y12[j], x_of_cities, y_of_cities) <
					dist(T[0], T[size - 1], x_of_cities, y_of_cities)) {
					G +=  dist(T[0], T[size - 1], x_of_cities, y_of_cities) -
							dist(y11, y12[j], x_of_cities, y_of_cities);
					for (k = 0; k < size && T[k] != y12[j]; k++);
					reverse(T, T + k);
				} else
					continue;

				//6 backtrack (c) alternate y1 in an inc order of length
				if (!btflag)
					break;
				G = 0;
				reverse(T, T + k);
			}
			//6 backtrack (d) alternate x1
			if (!btflag)
				break;
			reverse(T, T + size - 1);
		}
		//6 backtrack (e) choose a different t1
		if (!btflag)
			break;
		shift(size, 1, T);
	}

	/*
	 * return result
	 */
	*opt_tour_len = tour_len(size, opt_tour, x_of_cities, y_of_cities);

	free(x1);
	free(x2);
	free(y1);
	free(y2);
	free(T);
}

void swap(int* a, int* b) {
	int temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

void gen_rand_tour(int size, int* tour) {
	int i, pos;
	for (i = 0; i < size; i++) {
		tour[i] = i;
	}
	srand(time(NULL));
	for (i = 0; i < size - 1; i++) {
		pos = i + rand() % (size - i);
		swap(tour + i, tour + pos);
	}
}

double dist(int x, int y, double* x_of_cities, double* y_of_cities) {
	return sqrt((x_of_cities[x] - x_of_cities[y]) * (x_of_cities[x] - x_of_cities[y])
			+ (y_of_cities[x] - y_of_cities[y]) * (y_of_cities[x] - y_of_cities[y]));
}

double tour_len(int size, int* tour, double* x_of_cities, double* y_of_cities) {
	int i;
	double len = 0.0;
	for(i = 0; i < size - 1; i++) {
		len += dist(tour[i], tour[i + 1], x_of_cities, y_of_cities);
	}
	len += dist(tour[size - 1], tour[0], x_of_cities, y_of_cities);
	return len;
}

void shift(int size, int x, int** tour) {
	int* new_tour;
	int i;
	new_tour = (int*)malloc(size * sizeof(int));
	for(i = 0; i < size; i++) {
		new_tour[i] = (*tour)[(i + x) % size];
	}
	free((*tour));
	*tour = new_tour;
}

void reverse(int* a, int* b) {
	int* temp;
	int* p, q;

	temp = (int*)malloc((b - a) * sizeof(int));
	memcpy(temp, a, (b - a) * sizeof(int));
	for (p = a, q = temp + b - a - 1; p < b; p++, q--) {
		*p = *q;
	}

	free(temp);
}
