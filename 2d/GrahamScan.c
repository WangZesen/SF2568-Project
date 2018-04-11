#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)
#define RAND_MAX 2000000000
#define PI 3.14159265359
#define UPPER_MODE 1
#define LOWER_MODE 0
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct Index {
	int index;
	float angle;
};

struct Points {
	float x, y;
};

// For sorting positions
int compare_point(const void *a, const void *b) {
	if ((*(struct Points *)a).x < (*(struct Points *)b).x) {
		return 0;
	} else if ((*(struct Points *)a).x > (*(struct Points *)b).x) {
		return 1;
	} else if ((*(struct Points *)a).y < (*(struct Points *)b).y) {
		return 0;
	} else {
		return 1;
	}
}

// For sorting positions clockwisely
int compare_index(const void *a, const void *b) {
	return (*(struct Index *)a).angle > (*(struct Index *)b).angle; 
}

// Return the distance between two points
float dist(const struct Points *x, const struct Points *y) {
	return sqrt(((*x).x - (*y).x) * ((*x).x - (*y).x) + ((*x).y - (*y).y) * ((*x).y - (*y).y));
}

// Different version of cross product
float cross_product(const struct Points *point, const struct Index *localPoints, int x, int y, int z) {
	float x0 = point[localPoints[y].index].x - point[localPoints[x].index].x;
	float y0 = point[localPoints[y].index].y - point[localPoints[x].index].y;
	float x1 = point[localPoints[z].index].x - point[localPoints[y].index].x;
	float y1 = point[localPoints[z].index].y - point[localPoints[y].index].y;
	return x0 * y1 - y0 * x1;
}

float cross_product_index(const struct Points *point, int x, const int *storage, int y, int z, int len) {
	if ((y >= len) || (y < 0) || (z >= len) || (z < 0)) {
		return 1;
	} 
	float x0 = point[storage[y]].x - point[x].x;
	float y0 = point[storage[y]].y - point[x].y;
	float x1 = point[storage[z]].x - point[storage[y]].x;
	float y1 = point[storage[z]].y - point[storage[y]].y;
	return x0 * y1 - y0 * x1;
}

float cross_product_tri_index(const struct Points *point, int x, int y, int z) {
	float x0 = point[y].x - point[x].x;
	float y0 = point[y].y - point[x].y;
	float x1 = point[z].x - point[y].x;
	float y1 = point[z].y - point[y].y;
	return x0 * y1 - y0 * x1;
}

int left_evaluate(const struct Points *point, int index, int preIndex, const int *storage, int len, int mode) {
	int left, right, mid, ans;
	left = 0;
	right = len - 1;
	if (mode == UPPER_MODE) {
		while (left + 1 < right) {
			mid = (left + right) / 2;
			if (cross_product_index(point, index, storage, mid + 1, mid, len) <= 0) {
				right = mid;
			} else {
				left = mid;
			}
		}
		ans = len - 1;
		if (cross_product_index(point, index, storage, left + 1, left, len) <= 0) {
			ans = MIN(ans, left);
		}
		if (cross_product_index(point, index, storage, left + 2, left + 1, len) <= 0) {
			ans = MIN(ans, left + 1);
		}
		if (cross_product_index(point, index, storage, right + 1, right, len) <= 0) {
			ans = MIN(ans, right);
		}
		if (cross_product_tri_index(point, index, preIndex, storage[ans]) <= 0) {
			return 1;
		} else {
			return 0;
		}
	} else {
		while (left + 1 < right) {
			mid = (left + right) / 2;
			if (cross_product_index(point, index, storage, mid + 1, mid, len) >= 0) {
				right = mid;
			} else {
				left = mid;
			}
		}
		ans = len - 1;
		if (cross_product_index(point, index, storage, left + 1, left, len) >= 0) {
			ans = MIN(ans, left);
		}
		if (cross_product_index(point, index, storage, left + 2, left + 1, len) >= 0) {
			ans = MIN(ans, left + 1);
		}
		if (cross_product_index(point, index, storage, right + 1, right, len) >= 0) {
			ans = MIN(ans, right);
		}
		if (cross_product_tri_index(point, index, preIndex, storage[ans]) >= 0) {
			return 1;
		} else {
			return 0;
		}
	}
}

int right_evaluate(const struct Points *point, int index, int nxtIndex, const int *storage, int len, int mode) {
	int left, right, mid, ans;
	left = 0;
	right = len - 1;
	if (mode == UPPER_MODE) {
		while (left + 1 < right) {
			mid = (left + right) / 2;
			if (cross_product_index(point, index, storage, mid - 1, mid, len) <= 0) {
				right = mid;
			} else {
				left = mid;
			}
		}
		ans = 0;
		if (cross_product_index(point, index, storage, left - 1, left, len) >= 0) {
			ans = MAX(ans, left);
		}
		if (cross_product_index(point, index, storage, left, left + 1, len) >= 0) {
			ans = MAX(ans, left + 1);
		}
		if (cross_product_index(point, index, storage, right - 1, right, len) >= 0) {
			ans = MAX(ans, right);
		}
		if (cross_product_tri_index(point, index, nxtIndex, storage[ans]) >= 0) {
			return 1;
		} else {
			return 0;
		}
	} else {
		while (left + 1 < right) {
			mid = (left + right) / 2;
			if (cross_product_index(point, index, storage, mid - 1, mid, len) >= 0) {
				right = mid;
			} else {
				left = mid;
			}
		}
		ans = 0;
		if (cross_product_index(point, index, storage, left - 1, left, len) <= 0) {
			ans = MAX(ans, left);
		}
		if (cross_product_index(point, index, storage, left, left + 1, len) <= 0) {
			ans = MAX(ans, left + 1);
		}
		if (cross_product_index(point, index, storage, right - 1, right, len) <= 0) {
			ans = MAX(ans, right);
		}
		if (cross_product_tri_index(point, index, nxtIndex, storage[ans]) <= 0) {
			return 1;
		} else {
			return 0;
		}
	}
}

int main(int argc, char *argv[]) {

	MPI_Status status;
	int i, j, I, N, n, p, P, R, L, rc, shape;

	// Set Random Seed
	srand(1000);
	// Get Parameter N and Shape
	N = atoi(argv[1]);
	shape = atoi(argv[2]);
	// Initialization
	rc = MPI_Init(&argc, &argv);
	rc = MPI_Comm_size(MPI_COMM_WORLD, &P);
	rc = MPI_Comm_rank(MPI_COMM_WORLD, &p);
	
	if (N < P) {
		fprintf(stdout, "Too few discretization points...\n");
		exit(1);
	}
	// Create Pseudo-Random Positions
	struct Points *point = (struct Points *) malloc(N * sizeof(struct Points));
	switch (shape) {
		case 0: { // Square Shape
			
			for (i = 0; i < N; i++) {
				point[i].x = (rand() % RAND_MAX) / 100.0;
				point[i].y = (rand() % RAND_MAX) / 100.0;
			}
			break;
		}
		case 1: { // Disk Shape
			float r = RAND_MAX / 100.0 / 2.0;
			for (i = 0; i < N; i++) {
				float x = (rand() % RAND_MAX) / 100.0;
				float y = (rand() % RAND_MAX) / 100.0;
				while ((x - r) * (x - r) + (y - r) * (y - r) > r * r) {
					x = (rand() % RAND_MAX) / 100.0;
					y = (rand() % RAND_MAX) / 100.0;
				}
				point[i].x = x;
				point[i].y = y;
			}
			break;
		}
		case 2: { // Circle Shape
			float r = RAND_MAX / 100.0 / 2.0;
			for (i = 0; i < N; i++) {
				float angle = (rand() % 36000000) / 18000000.0 * PI;
				point[i].x = r + cos(angle) * r;
				point[i].y = r + sin(angle) * r;
			}
			break;
		}
	}
	qsort(point, N, sizeof(struct Points), compare_point);
	// Write Positions to File (For Small Cases)
	if ((p == 0) && (N < 3000)) {
		FILE *f;
		f = fopen("positions.txt", "w");
		for (i = 0; i < N; i++) {
			fprintf(f, "%.2f, %.2f\n", point[i].x, point[i].y);
		}
		fclose(f);
	}

	// Linear Load-balance Data Distribution
	I = (N + P - p - 1) / P;
	R = N % P;
	L = N / P;
	n = p * L + MIN(p, R);

	// Algorithm Starts Here
	MPI_Barrier(MPI_COMM_WORLD);
	double startTime = MPI_Wtime();

	// Find the Lowest Point
	float minY = RAND_MAX;
	int minNode = -1;
	for (i = n; i < I + n; i++) {
		if (point[i].y < minY) {
			minY = point[i].y;
			minNode = i;
		}
	}
	// Calculate Angle and Sort
	struct Index *localPoints = (struct Index *) malloc((I + 1) * sizeof(struct Index));
	for (i = n; i < I + n; i++) {
		localPoints[i - n].index = i;
		if (i == minNode) {
			localPoints[i - n].angle = -1;
		} else {
			if (point[i].x > point[minNode].x) {
				localPoints[i - n].angle = (point[i].y - minY) / dist(point + i, point + minNode);
			} else {
				localPoints[i - n].angle = 2 - (point[i].y - minY) / dist(point + i, point + minNode);
			}
		}
	}
	localPoints[I].index = minNode;
	qsort(localPoints, I, sizeof(struct Index), compare_index);
	// Graham Scan
	int *q = (int *) malloc(I * sizeof(int));
	q[0] = 0;
	q[1] = 1;
	int top = 1;
	for (i = 2; i < I + 1; i++) {
		while ((top > 0) && (cross_product(point, localPoints, q[top - 1], q[top], i) <= 0)) {
			top--;
		}
		q[++top] = i;
	}

	// Split Basic Convex Hull into Upper Bound and Lower Bound
	int *lower = (int *) malloc(top * sizeof(int));
	int *upper = (int *) malloc(top * sizeof(int));
	int lowerCount = 0;
	int upperCount = 0;
	int leftNode = 0;
	int rightNode = 0;
	for (i = 1; i < top; i++) {
		if (compare_point(localPoints[q[leftNode]].index + point, localPoints[q[i]].index + point)) {
			leftNode =  i;
		}
		if (!compare_point(localPoints[q[rightNode]].index + point, localPoints[q[i]].index + point)) {
			rightNode = i;
		}
	}
	i = leftNode;
	while (i != rightNode) {
		lower[lowerCount++] = localPoints[q[i]].index;
		i = (i + 1) % top;
	}
	lower[lowerCount++] = localPoints[q[i]].index;
	i = leftNode;
	while (i != rightNode) {
		upper[upperCount++] = localPoints[q[i]].index;
		i = (i - 1 + top) % top;
	}
	upper[upperCount++] = localPoints[q[i]].index;
	
	int maxLowerCount, maxUpperCount;
	MPI_Reduce(&lowerCount, &maxLowerCount, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&upperCount, &maxUpperCount, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(&maxLowerCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&maxUpperCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	lower = realloc(lower, maxLowerCount * sizeof(int));
	upper = realloc(upper, maxUpperCount * sizeof(int));
	
	int *lowerStorage[P];
	int *upperStorage[P];
	int *lowerSizes = (int *) malloc(P * sizeof(int));
	int *upperSizes = (int *) malloc(P * sizeof(int));
	
	for (i = 0; i < P; i++) {
		if (i != p) {
			lowerStorage[i] = (int *) malloc(maxLowerCount * sizeof(int));
			upperStorage[i] = (int *) malloc(maxUpperCount * sizeof(int));
		} else {
			lowerStorage[i] = lower;
			upperStorage[i] = upper;
			lowerSizes[i] = lowerCount;
			upperSizes[i] = upperCount;
		}
	}
	
	MPI_Request *requests = (MPI_Request *) malloc((P - 1) * 8 * sizeof(MPI_Request));
	int requestCount = 0;
	for (i = 0; i < P; i++) {
		if (i != p) {
			MPI_Isend(&lowerCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD, requests + (requestCount++));
			MPI_Irecv(lowerSizes + i, 1, MPI_INT, i, 0, MPI_COMM_WORLD, requests + (requestCount++));
			MPI_Isend(&upperCount, 1, MPI_INT, i, 1, MPI_COMM_WORLD, requests + (requestCount++));
			MPI_Irecv(upperSizes + i, 1, MPI_INT, i, 1, MPI_COMM_WORLD, requests + (requestCount++));

			MPI_Isend(lowerStorage[p], maxLowerCount, MPI_INT, i, 0, MPI_COMM_WORLD, requests + (requestCount++));
			MPI_Irecv(lowerStorage[i], maxLowerCount, MPI_INT, i, 0, MPI_COMM_WORLD, requests + (requestCount++));
			MPI_Isend(upperStorage[p], maxUpperCount, MPI_INT, i, 1, MPI_COMM_WORLD, requests + (requestCount++));
			MPI_Irecv(upperStorage[i], maxUpperCount, MPI_INT, i, 1, MPI_COMM_WORLD, requests + (requestCount++));
		}
	}
	// Wait Until All Requests Finish
	for (i = 0; i < requestCount; i++) {
		MPI_Wait(requests + i, &status);
	}

	// Construct Upper Bound on Complete Data
	int lp = 0;
	int rp = upperCount - 1;
	
	for (i = 0; i < P; i++) {
		if (i < p) {
			int left, right, mid, ans;
			left = 0; 
			right = upperCount - 1;
			while (left + 1 < right) {
				mid = (left + right) / 2;
				if ((mid == 0) || (left_evaluate(point, upperStorage[p][mid], upperStorage[p][mid - 1], upperStorage[i], upperSizes[i], UPPER_MODE))) {
					left = mid;
				} else {
					right = mid;
				}
			}
			ans = 0;
			if ((left == 0) || (left_evaluate(point, upperStorage[p][left], upperStorage[p][left - 1], upperStorage[i], upperSizes[i], UPPER_MODE))) {
				ans = MAX(ans, left);
			}
			if (left_evaluate(point, upperStorage[p][left + 1], upperStorage[p][left], upperStorage[i], upperSizes[i], UPPER_MODE)) {
				ans = MAX(ans, left + 1);
			}
			if ((right == 0) || (left_evaluate(point, upperStorage[p][right], upperStorage[p][right - 1], upperStorage[i], upperSizes[i], UPPER_MODE))) {
				ans = MAX(ans, right);
			}
			lp = MAX(lp, ans);
		} else if (i > p) {
			int left, right, mid, ans;
			left = 0;
			right = upperCount - 1;
			while (left + 1 < right) {
				mid = (left + right) / 2;
				if ((mid > upperCount - 2) || (right_evaluate(point, upperStorage[p][mid], upperStorage[p][mid + 1], upperStorage[i], upperSizes[i], UPPER_MODE))) {
					right = mid;
				} else {
					left = mid;
				}
			}
			ans = upperSizes[i] - 1;
			if ((left > upperCount - 2) || (right_evaluate(point, upperStorage[p][left], upperStorage[p][left + 1], upperStorage[i], upperSizes[i], UPPER_MODE))) {
				ans = MIN(ans, left);
			}
			if ((left > upperCount - 3) || (right_evaluate(point, upperStorage[p][left + 1], upperStorage[p][left + 2], upperStorage[i], upperSizes[i], UPPER_MODE))) {
				ans = MIN(ans, left + 1);
			}
			if ((right > upperCount - 2) || (right_evaluate(point, upperStorage[p][right], upperStorage[p][right + 1], upperStorage[i], upperSizes[i], UPPER_MODE))) {
				ans = MIN(ans, right);
			}
			rp = MIN(rp, ans);
		}
	}

	int *upperLP;
	int *upperRP;
	if (p == 0) {
		upperLP = (int *) malloc(P * sizeof(int));
		upperRP = (int *) malloc(P * sizeof(int));
	}
	MPI_Gather(&lp, 1, MPI_INT, upperLP, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&rp, 1, MPI_INT, upperRP, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int finalUpperCount = 0;
	int *finalUpperHull;

	if (p == 0) {
		
		int upperNodeCount = 1;
		finalUpperHull = (int *) malloc(sizeof(int));
	
		for (i = 0; i < P; i++) {
			upperNodeCount = finalUpperCount + upperRP[i] - upperLP[i] + 1;
			finalUpperHull = (int *) realloc(finalUpperHull, upperNodeCount * sizeof(int));
			for (j = upperLP[i]; j <= upperRP[i]; j++) {
				while ((finalUpperCount > 1) && (cross_product_tri_index(point, upperStorage[i][j], finalUpperHull[finalUpperCount - 1], finalUpperHull[finalUpperCount - 2]) <= 0)) {
					finalUpperCount--;
				}
				finalUpperHull[finalUpperCount++] = upperStorage[i][j];
			}
		}
	}

	// Construct Lower Bound on Complete Data
	lp = 0;
	rp = lowerCount - 1;
	
	for (i = 0; i < P; i++) {
		if (i < p) {
			int left, right, mid, ans;
			left = 0; 
			right = lowerCount - 1;
			while (left + 1 < right) {
				mid = (left + right) / 2;
				if ((mid == 0) || (left_evaluate(point, lowerStorage[p][mid], lowerStorage[p][mid - 1], lowerStorage[i], lowerSizes[i], LOWER_MODE))) {
					left = mid;
				} else {
					right = mid;
				}
			}
			ans = 0;
			if ((left == 0) || (left_evaluate(point, lowerStorage[p][left], lowerStorage[p][left - 1], lowerStorage[i], lowerSizes[i], LOWER_MODE))) {
				ans = MAX(ans, left);
			}
			if (left_evaluate(point, lowerStorage[p][left + 1], lowerStorage[p][left], lowerStorage[i], lowerSizes[i], LOWER_MODE)) {
				ans = MAX(ans, left + 1);
			}
			if ((right == 0) || (left_evaluate(point, lowerStorage[p][right], lowerStorage[p][right - 1], lowerStorage[i], lowerSizes[i], LOWER_MODE))) {
				ans = MAX(ans, right);
			}
			lp = MAX(lp, ans);
		} else if (i > p) {
			int left, right, mid, ans;
			left = 0;
			right = lowerCount - 1;
			while (left + 1 < right) {
				mid = (left + right) / 2;
				if ((mid > lowerCount - 2) || (right_evaluate(point, lowerStorage[p][mid], lowerStorage[p][mid + 1], lowerStorage[i], lowerSizes[i], LOWER_MODE))) {
					right = mid;
				} else {
					left = mid;
				}
			}
			ans = lowerSizes[i] - 1;
			if ((left > lowerCount - 2) || (right_evaluate(point, lowerStorage[p][left], lowerStorage[p][left + 1], lowerStorage[i], lowerSizes[i], LOWER_MODE))) {
				ans = MIN(ans, left);
			}
			if ((left > lowerCount - 3) || (right_evaluate(point, lowerStorage[p][left + 1], lowerStorage[p][left + 2], lowerStorage[i], lowerSizes[i], LOWER_MODE))) {
				ans = MIN(ans, left + 1);
			}
			if ((right > lowerCount - 2) || (right_evaluate(point, lowerStorage[p][right], lowerStorage[p][right + 1], lowerStorage[i], lowerSizes[i], LOWER_MODE))) {
				ans = MIN(ans, right);
			}
			rp = MIN(rp, ans);
		}
	}

	int *lowerLP;
	int *lowerRP;
	if (p == 0) {
		lowerLP = (int *) malloc(P * sizeof(int));
		lowerRP = (int *) malloc(P * sizeof(int));
	}
	MPI_Gather(&lp, 1, MPI_INT, lowerLP, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&rp, 1, MPI_INT, lowerRP, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int finalLowerCount = 0;
	int *finalLowerHull;

	if (p == 0) {
		// Final check for the computed convex hull
		int lowerNodeCount = 1;
		finalLowerHull = (int *) malloc(sizeof(int));
		for (i = 0; i < P; i++) {
			lowerNodeCount = finalLowerCount + lowerRP[i] - lowerLP[i] + 1;
			finalLowerHull = (int *) realloc(finalLowerHull, lowerNodeCount * sizeof(int));
			for (j = lowerLP[i]; j <= lowerRP[i]; j++) {
				while ((finalLowerCount > 1) && (cross_product_tri_index(point, lowerStorage[i][j], finalLowerHull[finalLowerCount - 1], finalLowerHull[finalLowerCount - 2]) > 0)) {
					finalLowerCount--;
				}
				finalLowerHull[finalLowerCount++] = lowerStorage[i][j];
			}
		}
	}
	
	// Algorithm Ends here
	MPI_Barrier(MPI_COMM_WORLD);
	if (p == 0) {
		double endTime = MPI_Wtime();
		printf("Execution Time: %e\n", endTime - startTime);

		if (N < 3000) {
			printf("%d, %d\n", 0, finalUpperCount - 1);
			for (j = 0; j < finalUpperCount; j++) {
				printf("%d\n", finalUpperHull[j]);
			}
			printf("%d, %d\n", 0, finalLowerCount - 1);
			for (j = 0; j < finalLowerCount; j++) {
				printf("%d\n", finalLowerHull[j]);
			}
		}
	}

	// Avoid Memory Leak
	
	if (p == 0) {
		free(finalUpperHull);
		free(finalLowerHull);
		free(upperLP);
		free(upperRP);
		free(lowerLP);
		free(lowerRP);
	}
	
	for (i = 0; i < P; i++) {
		free(lowerStorage[i]);
		free(upperStorage[i]);
	}
	
	free(q);
	free(point);
	free(requests);
	free(lowerSizes);
	free(upperSizes);
	free(localPoints);
	
	rc = MPI_Finalize();
	return 0;
}
