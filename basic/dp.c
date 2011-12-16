#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int n;
double x[25], y[25];
double dist[25][25];

int i, j, k;
int mark[1048577][25];
int prev[1048577][25];
double f[1048577][25];
int head, tail;
int tag[1048577*25][2];

int main() {
	scanf("%d", &n);
	for (i = 1; i <= n; i ++) {
		int id;
		double a, b;
		scanf("%d %lf %lf", &id, &a, &b);
		x[id] = a;
		y[id] = b;
	}
	for (i = 1; i <= n; i ++)
		for (j = 1; j <= n; j ++) {
			dist[i][j] = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
			dist[i][j] = (int)(dist[i][j] * 100);
		}

	head = 0;
	tail = 0;
	tag[0][0] = 1;
	tag[0][1] = 1;
	f[1][1] = 0;
	mark[1][1] = 1;
	prev[1][1] = -1;

	while (head <= tail) {
		int S = tag[head][0], u = tag[head][1];
		head++;
		int v;
		for (v = 2; v <= n; v ++)
			if ((S & (1 << (v-1))) == 0) {
				int R = (S | (1 << (v-1)));
				if (f[R][v] == 0 || f[R][v] > f[S][u] + dist[u][v]) {
					f[R][v] = f[S][u] + dist[u][v];
					prev[R][v] = head - 1;
				}
				if (!mark[R][v]) {
					mark[R][v] = 1;
					tail++;
					tag[tail][0] = R;
					tag[tail][1] = v;
				}
			}
	}

	double bestLen = 1e20;
	int V = (1 << n) - 1, v;
	for (i = 2; i <= n; i ++)
		if (f[V][i] + dist[1][i] < bestLen) {
			bestLen = f[V][i] + dist[1][i];
			v = i;
		}
	printf("%0.6lf\n", bestLen);
	int R = V;
	while (prev[R][v] != -1) {
		printf("%d ", v);
		int t = prev[R][v];
		R = tag[t][0];
		v = tag[t][1];
	}
	printf("%d\n", 1);
#ifdef DEBUG
	for (i = 1; i <= n; i ++) {
		for (j = 1; j <= n; j ++) printf("%0.3lf ", dist[i][j]);
		printf("\n");
	}


	int ans[] = {1, 14, 13, 12, 7, 6, 15, 5, 11, 9, 10, 16, 3, 2, 4, 8};
	double sum = 0;
	for (i = 1; i < 16; i ++) {
		sum += dist[ans[i-1]][ans[i]];
		printf("%lf\n", sum);
	}
	printf("%lf\n", sum + dist[ans[0]][ans[15]]);
#endif
		

	return 0;
}
