#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAXN 1000

int n;
double x[MAXN], y[MAXN];
double dist[MAXN][MAXN];

int i, j, k;

typedef struct Tour Tour;
struct Tour {
	int p[MAXN];
	double length;
};

Tour bestTour, cntTour;

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

int findIndex(int* p, int u) {
	int i;
	for (i = 1; i <= n; i ++)
		if (p[i] == u) return i;
	return 0;
}

void getInput() {
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
	for (i = 1; i <= n; i ++) sort(i, 1, n);
}

void init() {
	bestTour.length = 0;
	bestTour.p[0] = n;
	for (i = 1; i <= n; i ++) {
		bestTour.p[i] = i;
		bestTour.length += dist[bestTour.p[i]][bestTour.p[i-1]];
	}
}

int choose(int i, int lt) {
	if (i > n) return 0;
	// Choose xi(t[2i]).
	t[2*i] = cntTour.p[lt - 1];
	if (inY[t[2*i-1]][t[2*i]]) return 0;
	inX[t[2*i-1]][t[2*i]] = 1;
	inX[t[2*i]][t[2*i-1]] = 1;
	reorg(cntTour.p, lt);
	if (cntTour.length - dist[t[2*i-1]][t[2*i]] + dist[t[2*i]][t[1]] < bestTour.length) {
		cntTour.length -= (dist[t[2*i-1]][t[2*i]] - dist[t[2*i]][t[1]]); 
		bestTour = cntTour;
		return 1;
	}
	// Choose yi(t[2i+1]).
	int r, find = 0;
	for (r = 2; r <= n; r ++) {
		t[2*i+1] = nearest[t[2*i]][r];
		if (dist[t[2*i-1]][t[2*i]] - dist[t[2*i]][t[2*i+1]] > 0 && !inX[t[2*i]][t[2*i+1]]) {
			find = 1;
			break;
		}
	}
	if (!find) {
		inX[t[2*i-1]][t[2*i]] = 0;
		inX[t[2*i]][t[2*i-1]] = 0;
		return 0;
	}
	inY[t[2*i]][t[2*i+1]] = 1;
	inY[t[2*i+1]][t[2*i]] = 1;
	lt = findIndex(cntTour.p, t[2*i+1]);
	cntTour.length -= (dist[t[2*i-1]][t[2*i]] - dist[t[2*i]][t[2*i+1]]);
	if (choose(i + 1, lt)) return 1;
	else {
		inY[t[2*i]][t[2*i+1]] = 0;
		inY[t[2*i+1]][t[2*i]] = 0;
		return 0;
	}
}

int search(int i) {
	if (i == 1) {
		for (t[1] = 1; t[1] <= n; t[1] ++) {
			memset(inX, 0, sizeof(inX));
			memset(inY, 0, sizeof(inY));
			j = findIndex(bestTour.p, t[1]);
			shift(bestTour.p, j);
			t[2] = bestTour.p[2];
			inX[t[1]][t[2]] = inX[t[2]][t[1]] = 1;
			if (search(i + 1)) return 1;
		}
	} else if (i == 2) {
		for (j = 2; j <= n; j ++) {
			cntTour = bestTour;
			t[3] = nearest[t[2]][j];
			if (t[3] == t[2] || t[3] == t[1] || t[3] == cntTour.p[3]) continue;
			if (dist[t[1]][t[2]] - dist[t[2]][t[3]] <= 0) continue;
			t[4] = cntTour.p[findIndex(cntTour.p, t[3]) - 1];
			inY[t[2]][t[3]] = inY[t[3]][t[2]] = 1;
			inX[t[3]][t[4]] = inX[t[4]][t[3]] = 1;
			cntTour.length -= (dist[t[1]][t[2]] - dist[t[2]][t[3]]);
			reorg(cntTour.p, findIndex(cntTour.p, t[3]));
			if (cntTour.length - dist[t[3]][t[4]] + dist[t[4]][t[1]] < bestTour.length) {
				cntTour.length -= (dist[t[3]][t[4]] - dist[t[4]][t[1]]);
				bestTour = cntTour;
				return 1;
			}
			if (search(i + 1)) return 1;
			inY[t[2]][t[3]] = inY[t[3]][t[2]] = 0;
			inX[t[3]][t[4]] = inX[t[4]][t[3]] = 0;
		}
	} else if (i == 3) {
		Tour tourBuf = cntTour;
		for (k = 2; k <= n; k ++) {
			cntTour = tourBuf;
			t[5] = nearest[t[4]][k];
			if (t[5] == t[4] || t[5] == t[3] || t[5] == t[2] || t[5] == t[1] || t[5] == cntTour.p[3]) continue;
			if (dist[t[3]][t[4]] - dist[t[4]][t[5]] <= 0 || inX[t[4]][t[5]]) continue;
			inY[t[4]][t[5]] = inY[t[5]][t[4]] = 1;
			int lt = findIndex(cntTour.p, t[5]);
			cntTour.length -= (dist[t[3]][t[4]] - dist[t[4]][t[5]]);
			if (choose(3, lt)) return 1;
			inY[t[4]][t[5]] = inY[t[5]][t[4]] = 0;
		}
	}
	return 0;
}

int main() {
	// Get input.
	getInput();

	// Initialize
	init();

	while (1) {
		int update = search(1);
		if (!update) break;
	}

	printf("%0.5lf\n", bestTour.length);
	bestTour.p[0] = bestTour.p[n];
	double sum = 0;
	for (i = 1; i <= n; i ++) {
		printf("%d ", bestTour.p[i]);
		sum += dist[bestTour.p[i]][bestTour.p[i-1]];
	}
	printf("\n%0.5lf\n", sum);
	return 0;
}

