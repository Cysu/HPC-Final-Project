#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAXN 100

int n;
double x[MAXN], y[MAXN];
double dist[MAXN][MAXN];

int i, j, k;

typedef struct Tour Tour;
struct Tour {
	int p[MAXN];
	double length;
};

Tour bestTour;

int nearest[MAXN][MAXN];
int t[MAXN];
int inX[MAXN][MAXN];
int inY[MAXN][MAXN];

void shift(int* p, int j) {
	int tmp[MAXN];
	memcpy(tmp, p, sizeof(tmp));
	int i;
	for (i = 1; i <= n; i ++) {
		int t = i + j - 1;
		if (t > n) t -= n;
		p[i] = tmp[t];
	}
}

void reorg(int* p, int j) {
	int tmp[MAXN];
	memcpy(tmp, p, sizeof(tmp));
	int i;
	for (i = 2; i <= j - 1; i ++) {
		p[i] = tmp[j - 1 - (i - 2)];
	}
}


void sort(int u, int l, int r) {
	int i = l, j = r;
	double x = dist[u][nearest[u][(l + r) / 2]];
	while (i <= j) {
		while (dist[u][nearest[u][i]] < x) i ++;
		while (dist[u][nearest[u][j]] > x) j --;
		if (i <= j) {
			int t = nearest[u][i];
			nearest[u][i] = nearest[u][j];
			nearest[u][j] = t;
			i ++;
			j --;
		}
	}
	if (i < r) sort(u, i, r);
	if (l < j) sort(u, l, j);
}

int findIndex(int u) {
	int i;
	for (i = 1; i <= n; i ++)
		if (bestTour.p[i] == u) return i;
	return 0;
}

double calcPlan(int* p) {
	double sum = 0;
	p[0] = p[n];
	for (i = 1; i <= n; i ++) {
		sum += dist[p[i]][p[i-1]];
	}
	return sum;
}

int main() {
	// Get input.
	scanf("%d", &n);
	for (i = 0; i < n; i ++) {
		int id;
		double a, b;
		scanf("%d %lf %lf", &id, &a, &b);
		x[id] = a;
		y[id] = b;
	}
	for (i = 1; i <= n; i ++)
		for (j = 1; j <= n; j ++) {
			dist[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
			nearest[i][j] = j;
		}

	int u;
	for (u = 1; u <= n; u ++) sort(u, 1, n);

	// Initialize
	bestTour.length = 0;
	bestTour.p[0] = n;
	for (i = 1; i <= n; i ++) {
		bestTour.p[i] = i;
		bestTour.length += dist[bestTour.p[i]][bestTour.p[i-1]];
	}
	double cntLength = bestTour.length;

	while (1) {
		// Choose t[1].
		int tourUpdate = 0;
		printf("---------------------Update tour---------------------\n");
		for (t[1] = 1; t[1] <= n; t[1] ++) {
			printf("---------------------select another t1---------------------\n");
			printf("t[1] = %d\n", t[1]);

			// Initialize {X}, {Y}
			memset(inX, 0, sizeof(inX));
			memset(inY, 0, sizeof(inY));

			// Find and shift
			j = findIndex(t[1]);
			shift(bestTour.p, j);

			// TODO: Currently make t[2] = bestTour.p[2].
			// t[2] can equal to bestTour.p[n], too.
			t[2] = bestTour.p[2];
			inX[t[1]][t[2]] = inX[t[2]][t[1]] = 1;

			printf("t[2] = %d\n", t[2]);
			printf("%0.5lf\n", bestTour.length);
			for (i = 1; i <= n; i ++) printf("%d ", bestTour.p[i]);
			printf("\nrealLen = %0.5lf\n", calcPlan(bestTour.p));

			int inXbuf[MAXN][MAXN], inYbuf[MAXN][MAXN];
			Tour bestBuf = bestTour;
			memcpy(inXbuf, inX, sizeof(inX));
			memcpy(inYbuf, inY, sizeof(inY));

			// Choose y1(t[3]).
			for (j = 2; j <= n; j ++) {
				memcpy(inX, inXbuf, sizeof(inX));
				memcpy(inY, inYbuf, sizeof(inY));
				bestTour = bestBuf;
				cntLength = bestTour.length;

				t[3] = nearest[t[2]][j];
				if (t[3] == t[2] || t[3] == t[1] || t[3] == bestTour.p[3]) continue;

				printf("t[3] = %d\n", t[3]);

				if (dist[t[1]][t[2]] - dist[t[2]][t[3]] <= 0) continue;
				inY[t[2]][t[3]] = inY[t[3]][t[2]] = 1;
				t[4] = bestTour.p[findIndex(t[3]) - 1];
				inX[t[3]][t[4]] = inX[t[4]][t[3]] = 1;
				cntLength -= (dist[t[1]][t[2]] - dist[t[2]][t[3]]);

				printf("t[4] = %d\n", t[4]);

				// Reorganize tour.
				reorg(bestTour.p, findIndex(t[3]));
				if (cntLength - dist[t[3]][t[4]] + dist[t[4]][t[1]] < bestTour.length) {
					bestTour.length = cntLength - dist[t[3]][t[4]] + dist[t[4]][t[1]];
					cntLength = bestTour.length;
					tourUpdate = 1;
					printf("%0.5lf\n", bestTour.length);
					for (i = 1; i <= n; i ++) printf("%d ", bestTour.p[i]);
					printf("\nrealLen = %0.5lf\n", calcPlan(bestTour.p));
					break;
				}

				int inXbuf[MAXN][MAXN], inYbuf[MAXN][MAXN];
				Tour bestBuf = bestTour;
				memcpy(inXbuf, inX, sizeof(inX));
				memcpy(inYbuf, inY, sizeof(inY));

				// Choose y2(t[5]).
				for (k = 2; k <= n; k ++) {
					memcpy(inX, inXbuf, sizeof(inX));
					memcpy(inY, inYbuf, sizeof(inY));
					bestTour = bestBuf;
					cntLength = bestTour.length;

					t[5] = nearest[t[4]][k];
					if (t[5] == t[4] || t[5] == t[3] || t[5] == t[2] || t[5] == t[1] || t[5] == bestTour.p[3]) continue;

					printf("t[5] = %d\n", t[5]);

					if (dist[t[3]][t[4]] - dist[t[4]][t[5]] <= 0 || inX[t[4]][t[5]]) continue;
					inY[t[4]][t[5]] = inY[t[5]][t[4]] = 1;
					int lt = findIndex(t[5]);
					cntLength -= (dist[t[3]][t[4]] - dist[t[4]][t[5]]);

					// Unique choice.
					for (i = 3; i <= n; i ++) {
						// Choose xi(t[2i]).
						t[2*i] = bestTour.p[lt - 1];
						if (inY[t[2*i-1]][t[2*i]]) break;
						printf("t[%d] = %d\n", 2*i, t[2*i]);
						inX[t[2*i-1]][t[2*i]] = 1;
						inX[t[2*i]][t[2*i-1]] = 1;
						reorg(bestTour.p, lt);
						if (cntLength - dist[t[2*i-1]][t[2*i]] + dist[t[2*i]][t[1]] < bestTour.length) {
							bestTour.length = cntLength - dist[t[2*i-1]][t[2*i]] + dist[t[2*i]][t[1]]; 
							cntLength = bestTour.length;
							tourUpdate = 1;
							printf("%0.5lf\n", bestTour.length);
							for (i = 1; i <= n; i ++) printf("%d ", bestTour.p[i]);
							printf("\nrealLen = %0.5lf\n", calcPlan(bestTour.p));
							break;
						}
						// Choose yi(t[2i+1]).
						int r, find = 0;
						for (r = 2; r <= n; r ++) {
							t[2*i+1] = nearest[t[2*i]][r];
							if (dist[t[2*i-1]][t[2*i]] - dist[t[2*i]][t[2*i+1]] > 0 && inX[t[2*i]][t[2*i+1]]) {
								find = 1;
								break;
							}
						}
						if (!find) break;
						printf("t[%d] = %d\n", 2*i+1, t[2*i+1]);
						inY[t[2*i]][t[2*i+1]] = 1;
						inY[t[2*i+1]][t[2*i]] = 1;
						lt = findIndex(t[2*i+1]);
						cntLength -= (dist[t[2*i-1]][t[2*i]] - dist[t[2*i]][t[2*i+1]]);
					}
					if (tourUpdate) break;
				}
				if (tourUpdate) break;
			}
			if (tourUpdate) break;
		}
		if (!tourUpdate) break;
	}

	printf("%0.5lf\n", bestTour.length);
	for (i = 1; i <= n; i ++)
		printf("%d ", bestTour.p[i]);
	printf("\n");
	return 0;
}

