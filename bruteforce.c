#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAXN 100

int n;
double x[MAXN], y[MAXN];
double e[MAXN * MAXN];

int i, j, k;
int mark[MAXN * MAXN];
int plan[MAXN];
int vis[MAXN];

double bestLen;
int bestPlan[MAXN];

int totCut;

int compare(const void* a, const void* b) {
	return (*(double*)a - *(double*)b);
}

int getEdgeId(int i, int j) {
	if (j == i) return -1;
	if (j < i) {
		int t = i;
		i = j;
		j = t;
	}
	return (2*n-i)*(i-1)/2+j-i-1;
}

void search(int i, double len) {
	if (i > n) {
		len = len + e[getEdgeId(plan[1], plan[n])];
		if (len < bestLen) {
			bestLen = len;
			memcpy(bestPlan, plan, sizeof(plan));
		}
		return;
	}
	int v;
	for (v = 1; v <= n; v ++) {
		if (!vis[v]) {
			// heuristic
			double h = 0;
			int count = 0;
			int eId = getEdgeId(plan[i-1], v);
			for (j = 0; j < n*(n-1)/2 && count < n-i+1; j ++)
				if (!mark[j] && j != eId) {
					count ++;
					h += e[j];
				}
			if (len + e[eId] + h >= bestLen) {
				totCut ++;
				continue;
			}
			// search
			vis[v] = 1;
			mark[eId] = 1;
			plan[i] = v;
			search(i + 1, len + e[eId]);
			mark[eId] = 0;
			vis[v] = 0;
		}
	}
}


int main() {
	scanf("%d", &n);
	for (i = 1; i <= n; i ++) {
		int id;
		double a, b;
		scanf("%d %lf %lf", &id, &a, &b);
		x[id] = a;
		y[id] = b;
	}
	int id = 0;
	for (i = 1; i <= n; i ++)
		for (j = i + 1; j <= n; j ++)
			e[id++] = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
	qsort(e, id, sizeof(double), compare);
#ifdef DEBUG
	for (i = 1; i <= n; i ++) {
		for (j = 1; j <= n; j ++) {
			int eId = getEdgeId(i, j);
			if (eId == -1) printf("- ");
			else printf("%0.3lf ", e[eId]);
		}
		printf("\n");
	}
#endif
	memset(mark, 0, sizeof(mark));
	plan[1] = 1;
	vis[1] = 1;
	bestLen = 1e20;
	totCut = 0;
	search(2, 0);

	printf("totCut = %d\n", totCut);
	printf("bestLen = %0.6lf\n", bestLen);
	for (i = 1; i <= n; i ++) printf("%d ", bestPlan[i]);
	printf("\n");
	return 0;
}




