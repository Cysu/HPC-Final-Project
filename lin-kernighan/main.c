#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include "my_conf.h"
#include "lin_kernighan.h"

int main(int argc, char** argv) {
	/*
	 * Generate or read data set
	 */
	double x_of_cities[MAX_SIZE];
	double y_of_cities[MAX_SIZE];
	int size;
	int i;
	int temp;
	int opt_tour[MAX_SIZE];
	double opt_tour_len;
	FILE* fin;
	if (argc < 2 || (strcmp(argv[1], "-r") != 0 && strcmp(argv[1], "-f") != 0)) {
		printf("-r size # randomly generate test data of the given size\n");
		printf("-f filename # read test data from file\n");
	} else if (strcmp(argv[1], "-r") == 0) {
		size = atoi(argv[2]);
		srand(time(NULL));
		for (i = 0; i < size; i++) {
			x_of_cities[i] = (double)(rand() % MAX_X);
			y_of_cities[i] = (double)(rand() % MAX_Y);
		}
	} else {
		scanf("%d", &size);
		for(i = 0; i < 16; i++) {
			scanf("%d %lf %lf", &temp, x_of_cities + i, y_of_cities + i);
		}
	}
	/*
	 * Find the optimized tour
	 */
	lin_kernighan(size, x_of_cities, y_of_cities, opt_tour, &opt_tour_len);

	/*
	 * output
	 */
	printf("Optimized tour:");
	for(i = 0; i < size; i++) {
		printf(" %d", opt_tour[i]);
	}
	printf("\nOptimized tour length:%lf", opt_tour_len);

	getchar();
	return 0;
}
