#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAXN 1000

int n;
double x[MAXN], y[MAXN];
double dist[MAXN][MAXN];
int route[MAXN];

void readln(FILE* file) {
	char c;
	while (1) {
		fscanf(file, "%c", &c);
		if (c == '\n') break;
	}
}

int main(int argc, char* argv[]) {
	FILE* dataFile = fopen(argv[1], "r");
	FILE* ansFile = fopen(argv[2], "r");

	int i, j;
	fscanf(dataFile, "%d", &n);
	for (i = 0; i < n; i ++) {
		int id;
		double a, b;
		fscanf(dataFile, "%d %lf %lf", &id, &a, &b);
		x[id] = a;
		y[id] = b;
	}
	for (i = 1; i <= n; i ++)
		for (j = 1; j <= n; j ++) 
			dist[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
	fclose(dataFile);

	for (i = 0; i < 5; i ++) readln(ansFile);
	for (i = 1; i <= n; i ++)
		fscanf(ansFile, "%d", &route[i]);
	fclose(ansFile);

	double len = 0;
	route[0] = route[n];
	for (i = 1; i <= n; i ++) {
		printf("%d ", route[i]);
		len += dist[route[i]][route[i-1]];
	}
	printf("\nlen = %lf\n", len);
	return 0;
}
