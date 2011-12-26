#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define MAXN 226
#define LOOPTIME 100

int n;
double x[MAXN], y[MAXN];
double dist[MAXN][MAXN];

typedef struct Tour Tour;
struct Tour {
	int p[MAXN];
	double length;
};

Tour bestTour, cntTour, recvTour;
Tour bestTours[LOOPTIME];

int nearest[MAXN][MAXN];
int t[MAXN];
int inX[MAXN][MAXN];
int inY[MAXN][MAXN];
int initY[MAXN][MAXN];
int fixed[MAXN][MAXN];
int fixedEdges[MAXN][2];

int loopRound;
int myid, numprocs;
MPI_Status status;

void swap(int* a, int* b) {
	int temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

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
	int i;
	memcpy(tmp, p, sizeof(tmp));
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
	int i, j;
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
	int i, j;
	int pos;
	bestTour.length = 0;
	for (i = 1; i <= n; i ++) {
		bestTour.p[i] = i;
	}
	srand(time(NULL)+myid+loopRound*16);
	for (i = 1; i < n; i++) {
		pos = i + rand() % (n - i + 1);
		swap(bestTour.p + i, bestTour.p + pos);
	}
	bestTour.p[0] = bestTour.p[n];
	for (i = 1; i <= n; i++) {
		bestTour.length += dist[bestTour.p[i]][bestTour.p[i-1]];
	}
	memset(initY, 0, sizeof(initY));
	for (i = 0; i < MAXN; i ++) {
		if (fixedEdges[i][0] == 0 && fixedEdges[i][1] == 0) break;
		int u = fixedEdges[i][0], v= fixedEdges[i][1];
		initY[u][v] = initY[v][u] = 1;
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
	int j, k;
	if (i == 1) {
		for (t[1] = 1; t[1] <= n; t[1] ++) {
			memset(inX, 0, sizeof(inX));
			memcpy(inY, initY, sizeof(initY));
			j = findIndex(bestTour.p, t[1]);
			shift(bestTour.p, j);
			t[2] = bestTour.p[2];
			if (inY[t[1]][t[2]]) continue;
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
			if (inY[t[3]][t[4]]) continue;
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

int inTour(int* p, int u, int v) {
	int i;
	p[0] = p[n];
	for (i = 1; i <= n; i ++)
		if ((p[i] == u && p[i-1] == v) || (p[i] == v && p[i-1] == u))
			return 1;
	return 0;
}

int main(int argc, char** argv) {
	double begin, end;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	//init
	if (myid == 0) {
		// Get input.
		getInput();
		begin = MPI_Wtime();
	}

	//brodcast n
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//brodcast dist
	MPI_Bcast(dist, MAXN * MAXN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//brodcast nearest
	MPI_Bcast(nearest, MAXN * MAXN, MPI_INT, 0, MPI_COMM_WORLD);

	Tour gBestTour;

	int i, j, k;
	for (loopRound = 0; loopRound < LOOPTIME; loopRound ++) {
		Tour localBestTour;
		for (i = 1; i <= 1; i ++) {
			init();
			if (myid == 0 && loopRound == 0 && i == 1) gBestTour = bestTour;
			if (i == 1) localBestTour = bestTour;
			while (1) {
				int update = search(1);
				if (!update) break;
			}
			if (bestTour.length < localBestTour.length) {
				localBestTour = bestTour;
			}
		}
		// Reduce
		int changed = 0;
		if (myid != 0) {
			MPI_Send(localBestTour.p, MAXN, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&localBestTour.length, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			MPI_Recv(&changed, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
		} else {
			if (localBestTour.length < gBestTour.length)
				gBestTour = localBestTour;
			for (i = 1; i <= n; i ++)
				for (j = i + 1; j <= n; j ++)
					if (inTour(localBestTour.p, i, j)) {
						fixed[i][j] ++;
						fixed[j][i] ++;
					}
			for (k = 0; k < numprocs - 1; k++) {
				MPI_Recv(recvTour.p, MAXN, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&recvTour.length, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
				if (recvTour.length < gBestTour.length)	gBestTour = recvTour;
				if (recvTour.length < localBestTour.length) localBestTour = recvTour;
				for (i = 1; i <= n; i ++)
					for (j = i + 1; j <= n; j ++)
						if (inTour(recvTour.p, i, j)) {
							fixed[i][j] ++;
							fixed[j][i] ++;
						}
			}
			printf("loopRound %d, gBestTour.length = %lf\n", loopRound, gBestTour.length);
			bestTours[loopRound] = gBestTour;
			if (loopRound > 0 && fabs(bestTours[loopRound].length - bestTours[loopRound-1].length) > 1e-3) changed = 1;
			for (j = 1; j < numprocs; j++)
				MPI_Send(&changed, 1, MPI_INT, j, 2, MPI_COMM_WORLD);
			k = 0;
			memset(fixedEdges, 0, sizeof(fixedEdges));
			for (i = 1; i <= n; i ++)
				for (j = i + 1; j <= n; j ++) {
					if (fixed[i][j] == (loopRound + 1) * numprocs) {
						fixedEdges[k][0] = i;
						fixedEdges[k][1] = j;
						k ++;
					}
				}
		}
		MPI_Bcast(fixedEdges, MAXN * 2, MPI_INT, 0, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (myid == 0) {
		end = MPI_Wtime();
		printf("\nresult:\n");
		printf("time = %lf\n", end - begin);
		printf("len = %lf\n", gBestTour.length);
		for (i = 1; i <= n; i ++) {
			printf("%d ", gBestTour.p[i]);
		}
		printf("\n");
	}

	MPI_Finalize();
	return 0;
}

