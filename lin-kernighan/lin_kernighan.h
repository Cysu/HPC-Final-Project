#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>

#include<my_conf.h>
/*
 *
 */
void lin_kernighan(int size, int* x_of_cities, int* y_of_cities, int* opt_tour, double* opt_tour_len);

/*
 * swap two int
 */
void swap(int* a, int* b) {
	int temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

/*
 * generate a random array of 0..size-1
 */
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

/*
 * compute the dist of two cities
 */
inline double dist(int x, int y, int* x_of_cities, int* y_of_cities) {
	return sqrt((x_of_cities[x] - x_of_cities[y]) * (x_of_cities[x] - x_of_cities[y])
			+ (y_of_cities[x] - y_of_cities[y]) * (y_of_cities[x] - y_of_cities[y]));
}

/*
 * compute the length of a tour
 */
double tour_len(int size, int* tour, int* x_of_cities, int* y_of_cities) {
	int i;
	double len = 0.0;
	for(i = 0; i < size - 1; i++) {
		len += dist(tour[i], tour[i + 1], x_of_cities, y_of_cities);
	}
	len += dis(tour[size - 1], tour[0], x_of_cities, y_of_cities);
	return len;
}
