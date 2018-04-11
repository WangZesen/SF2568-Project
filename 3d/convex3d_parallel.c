#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)
#define RAND_MAX 2000000000
#define PI 3.14159265359
#define EPS 1e-15
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

// Struct Point and Its Auxiliary Functions
struct Point {
	double x, y, z;
};

// Copy a point
void copy_point(struct Point *x, const struct Point *other) {
	x->x = other->x;
	x->y = other->y;
	x->z = other->z;
}

// Return a unit vector
struct Point unit_vector() {
	struct Point x;
	x.x = x.y = x.z = 1;
	return x;
}

// Return Vector from point x to point y
struct Point subtract(const struct Point *x, const struct Point *y) {
	struct Point z;
	z.x = y->x - x->x;
	z.y = y->y - x->y;
	z.z = y->z - x->z;
	return z;
}

// return square of a number
double sqr(double x) {
	return x * x;
}

// return the length of a vector
double dist(const struct Point point) {
	return sqrt(sqr(point.x) + sqr(point.y) + sqr(point.z));
}

// return the length of a vector in 2d
double dist_2d(const struct Point point) {
	return sqrt(sqr(point.x) + sqr(point.y));
}

// return dot product of two vectors
double dot_prod(const struct Point *a, const struct Point *b) {
	return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

// return a normalized vector
struct Point normalize(const struct Point *a) {
	struct Point b;
	double len = dist(*a);
	if (len > 0) {
		b.x = a->x / len;
		b.y = a->y / len;
		b.z = a->z / len;
	} else {
		copy_point(&b, a);
	}
	return b;
}

// return the dot product of two vector
struct Point cross_prod(const struct Point *a, const struct Point *b) {
	struct Point z;
	z.x = a->y * b->z - a->z * b->y;
	z.y = a->z * b->x - a->x * b->z;
	z.z = a->x * b->y - a->y * b->x;
	return z;
}

// return the dot product of two vectors, and the second vector is normalized
double normal_dot_prod(const struct Point *a, const struct Point *b) {
	struct Point normalB = normalize(b);
	return dot_prod(a, &normalB);
}

// return the dot product of two vectors in 2d, and the second vector is normalized
double normal_dot_prod_2d(const struct Point *a, const struct Point *b, const struct Point *c) {
	struct Point supp1 = subtract(a, b);
	struct Point supp2 = subtract(b, c);
	double len = dist_2d(supp1);
	supp1.x = supp1.x / len;
	supp1.y = supp1.y / len;
	len = dist_2d(supp2);
	supp2.x = supp2.x / len;
	supp2.y = supp2.y / len;
	
	return supp1.x * supp2.x + supp1.y * supp2.y;
}

// return the cross product of two vectors in 2d
double cross_prod_2d(const struct Point *a, const struct Point *b, const struct Point *c) {
	struct Point supp1 = subtract(a, b);
	struct Point supp2 = subtract(b, c);
	return supp1.x * supp2.y - supp1.y * supp2.x;
}

// return true if three points are collinear
int on_line(const struct Point *a, const struct Point *b, const struct Point *c) {
	struct Point supp1 = subtract(a, b);
	struct Point supp2 = subtract(a, c);
	struct Point normal = cross_prod(&supp1, &supp2);
	struct Point unit = unit_vector();
	
	double area = dot_prod(&normal, &unit);
	if ((area < EPS) && (area > -EPS)) {
		return 1;
	} else {
		return 0;
	}
}

// used for sorting the positions of points
int compare_point(const void *a, const void *b) {
	if ( ((struct Point *) a)->x > ((struct Point *) b)->x ) {
		return 1;
	} else if ( ((struct Point *) a)->x < ((struct Point *) b)->x ) {
		return 0;
	} else if ( ((struct Point *) a)->y > ((struct Point *) b)->y ) {
		return 1;
	} else if ( ((struct Point *) a)->y < ((struct Point *) b)->y ) {
		return 0;
	} else if ( ((struct Point *) a)->z > ((struct Point *) b)->z ) {
		return 1;
	} else {
		return 0;
	} 
}

// Struct for Edge between Two Points
struct Edge {
	struct Edge *next;
	int id, dst, bound;
};

// add an edge endding with <id> and linking to <dst>
void construct_edge(struct Edge *start, int dst, int id) {
	struct Edge *newEdge = (struct Edge *) malloc(sizeof(struct Edge));
	newEdge->next = start->next;
	newEdge->id = id;
	newEdge->dst = dst;
	newEdge->bound = 0;
	start->next = newEdge;
}

// return the id of the edge endding with <id>
int find_edge_id(const struct Edge *start, int dst) {
	struct Edge *pivot = start->next;
	while (pivot != NULL) {
		if (pivot->dst == dst) {
			return pivot->id;
		}
		pivot = pivot->next;
	}
	return -1;
}

// return whether an edge is a bound edge
int find_edge_bound(const struct Edge *start, int dst) {
	struct Edge *pivot = start->next;
	while (pivot != NULL) {
		if (pivot->dst == dst) {
			return pivot->bound;
		}
		pivot = pivot->next;
	}
	return -1;
}

// set bound property of the edge linking to <dst> as <bound>
int set_edge_bound(const struct Edge *start, int dst, int bound) {
	struct Edge *pivot = start->next;
	while (pivot != NULL) {
		if (pivot->dst == dst) {
			pivot->bound = bound;
			return pivot->id;
		}
		pivot = pivot->next;
	}
	printf("Set Edge Bound Failed\n");
	return -1;
}

// change id property of the edge linking to <dst> as <id>
void edit_edge_id(struct Edge *start, int dst, int id) {
	struct Edge *pivot = start->next;
	while (pivot != NULL) {
		if (pivot->dst == dst) {
			pivot->id = id;
			return;
		}
		pivot = pivot->next;
	}
	printf("Edit Failed\n");
}

// change id property and dst property of the edge linking to <dst> as <id> and <newDst>
void edit_edge_dst_id(struct Edge *start, int dst, int newDst, int id) {
	struct Edge *pivot = start->next;
	while (pivot != NULL) {
		if (pivot->dst == dst) {
			pivot->id = id;
			pivot->dst = newDst;
			return;
		}
		pivot = pivot->next;
	}
	printf("Edit Failed\n");
}

// Struct for Facet
struct Facet {
	int a, b, c; // 3 Nodes on the facet
	int id; // facet id
	int valid; // whether facet belongs to convex hull
};

// return a facet
struct Facet construct_facet(int a, int b, int c, int id, int valid) {
	struct Facet facet;
	facet.a = a;
	facet.b = b;
	facet.c = c;
	facet.id = id;
	facet.valid = valid;
	return facet;
}

// change the property of a facet
void edit_facet(struct Facet *facet, int a, int b, int c, int id, int valid) {
	facet->a = a;
	facet->b = b;
	facet->c = c;
	facet->id = id;
	facet->valid = valid;
}

// copy a facet
void copy_facet(struct Facet *facet, const struct Facet *other) {
	facet->a = other->a;
	facet->b = other->b;
	facet->c = other->c;
	facet->id = other->id;
	facet->valid = other->valid;
}

// return the normal vector of a facet
struct Point facet_normal(const struct Point *points, const struct Facet *facet) {
	struct Point supp1 = subtract(points + facet->b, points + facet->a);
	struct Point supp2 = subtract(points + facet->c, points + facet->a);
	struct Point normal = cross_prod(&supp1, &supp2);	
	return normal;
}

// return the normal vector of a temporary facet
struct Point tem_facet_normal(const struct Point *a, const struct Point *b, const struct Point *c) {
	struct Point supp1 = subtract(b, a);
	struct Point supp2 = subtract(c, a);
	struct Point normal = cross_prod(&supp1, &supp2);	
	return normal;
}

// return true if <points[index]> is on facet <facet>
int on_facet(const struct Point *points, const struct Facet *facet, int index) {

	struct Point supp3 = subtract(points + index, points + facet->a);
	struct Point normal = facet_normal(points, facet);
	double proj = dot_prod(&normal, &supp3);
	if ((proj < EPS) && (proj > -EPS)) {
		return 1;
	} else {
		return 0;
	}
}

// return true if <points[index]> is above facet <facet>
int above_facet(const struct Point *points, const struct Facet *facet, int index) {
	struct Point supp3 = subtract(points + index, points + facet->a);
	struct Point normal = facet_normal(points, facet);
	double proj = dot_prod(&normal, &supp3);
	if (proj > EPS) {
		return 1;
	} else {
		return 0;
	}
}

// construct 3 edges for facet <facet>
void construct_edges(const struct Point *points, struct Edge *edges, struct Facet *facet) {
	construct_edge(edges + facet->a, facet->b, facet->id);
	construct_edge(edges + facet->b, facet->c, facet->id);
	construct_edge(edges + facet->c, facet->a, facet->id);
}

// change the valid property of some facets that can be reached by <points[index]>
void change_facet(const struct Point *points, struct Edge *edges, struct Facet *facets, struct Facet *facet, int *cnt, int index) {
	// delete current facet	
	facet->valid = 0; 
	int id;
	id = find_edge_id(edges + facet->b, facet->a);
	if (facets[id].valid) {
		if (above_facet(points, facets + id, index)) {
			change_facet(points, edges, facets, facets + id, cnt, index);
		} else {
			edit_facet(facets + *cnt, facet->a, facet->b, index, *cnt, 1);
			construct_edge(edges + facet->b, index, *cnt);
			construct_edge(edges + index, facet->a, *cnt);
			edit_edge_id(edges + facet->a, facet->b, *cnt);
			(*cnt)++;
		}
	}
	id = find_edge_id(edges + facet->c, facet->b);
	if (facets[id].valid) {
		if (above_facet(points, facets + id, index)) {
			change_facet(points, edges, facets, facets + id, cnt, index);
		} else {
			edit_facet(facets + *cnt, facet->b, facet->c, index, *cnt, 1);
			construct_edge(edges + facet->c, index, *cnt);
			construct_edge(edges + index, facet->b, *cnt);
			edit_edge_id(edges + facet->b, facet->c, *cnt);
			(*cnt)++;
		}
	}
	id = find_edge_id(edges + facet->a, facet->c);
	if (facets[id].valid) {
		if (above_facet(points, facets + id, index)) {
			change_facet(points, edges, facets, facets + id, cnt, index);
		} else {
			edit_facet(facets + *cnt, facet->c, facet->a, index, *cnt, 1);
			construct_edge(edges + facet->a, index, *cnt);
			construct_edge(edges + index, facet->c, *cnt);
			edit_edge_id(edges + facet->c, facet->a, *cnt);
			(*cnt)++;
		}
	}		
}

// Generate Random Positions
void generate(struct Point *points, int N, int shape) {
	int i;
	double scale = 1000.0;
	switch (shape) {
		case 0: { // Square
			for (i = 0; i < N; i++) {
				points[i].x = (rand() % RAND_MAX) / scale;
				points[i].y = (rand() % RAND_MAX) / scale;
				points[i].z = (rand() % RAND_MAX) / scale;
			}
			break;
		}
		case 1: { // Ball
			double r = RAND_MAX / scale / 2.0;
			struct Point center;
			center.x = r;
			center.y = r;
			center.z = r;
			for (i = 0; i < N; i++) {
				struct Point tmp;
				tmp.x = (rand() % RAND_MAX) / scale;
				tmp.y = (rand() % RAND_MAX) / scale;
				tmp.z = (rand() % RAND_MAX) / scale;	
				while (dist(subtract(&tmp, &center)) > r) {
					tmp.x = (rand() % RAND_MAX) / scale;
					tmp.y = (rand() % RAND_MAX) / scale;
					tmp.z = (rand() % RAND_MAX) / scale;
				}
				copy_point(&points[i], &tmp);
			}
			break;
		}
		case 2: { // Sphere
			double r = RAND_MAX / scale / 2.0;
			for (i = 0; i < N; i++) {
				double u = (rand() % RAND_MAX) / (double) RAND_MAX;
				double v = (rand() % RAND_MAX) / (double) RAND_MAX;
				double theta = 2 * PI * u;
				double phi = acos(2 * v - 1);
				points[i].x = sin(theta) * sin(phi) * r;
				points[i].y = cos(theta) * sin(phi) * r;
				points[i].z = cos(phi) * r;
			}
			break;
		}
	}
}

// Print the position of <points[index]>
void print_a_point(FILE *f, const struct Point *points, int index) {
	fprintf(f, "%.2f, %.2f, %.2f, ", points[index].x, points[index].y, points[index].z);
}

// Write Positions to File
void print_positions(const struct Point *points, int N) {
	int i;
	FILE *f;
	f = fopen("positions.txt", "w");
	for (i = 0; i < N; i++) {
		print_a_point(f, points, i);
		fprintf(f, "\n");
	}
	fclose(f);
}

// Print the positions of facets
void print_facets_write(struct Point *points, struct Facet *facets, int startCnt, int cnt) {
	int i;
	FILE *f;
	f = fopen("output.txt", "w");
	fprintf(f, "----------\n");
	for (i = startCnt; i < cnt; i++) {
		if (facets[i].valid) {
			print_a_point(f, points, facets[i].a);
			print_a_point(f, points, facets[i].b);
			print_a_point(f, points, facets[i].c);
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

// Print the positions of facets
void print_facets_append(struct Point *points, struct Facet *facets, int startCnt, int cnt) {
	int i;
	FILE *f;
	f = fopen("output.txt", "a");
	fprintf(f, "----------\n");
	for (i = startCnt; i < cnt; i++) {
		if (facets[i].valid) {
			print_a_point(f, points, facets[i].a);
			print_a_point(f, points, facets[i].b);
			print_a_point(f, points, facets[i].c);
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

// Print the positions of 2d convex hull (for test)
void print_chs_write() {
	FILE *f;
	f = fopen("chs.txt", "w");
	fclose(f);
}

// Print the positions of 2d convex hull (for test)
void print_chs_append(const struct Point *points, int *ch) {
	int i;
	FILE *f;
	f = fopen("chs.txt", "a");
	fprintf(f, "----------\n");
	i = 0;
	do {
		print_a_point(f, points, ch[i++]);
		fprintf(f, "\n");
	} while (ch[i] != ch[0]);
	fclose(f);
}

// Construct a convex hull from scratch
void construct(struct Point *points, struct Edge *edges, struct Facet *facets, int l, int r, int *cnt) {
	int i, j, startCnt;
	startCnt = *cnt;
	
	int *selected = (int *) calloc(r - l + 1, sizeof(int));
	int s[4];
	
	s[0] = l;
	s[1] = l + 1;
	
	selected[0] = 1;
	selected[1] = 1;
	
	// Find Three Non-colinear Points
	for (i = 2 + l; i < r; i++) {
		if (!on_line(points + l, points + l + 1, points + i)) {
			selected[i - l] = 1;
			s[2] = i;
			break;
		}
	}
	
	// Find Four Non-coplanar Points
	struct Facet facet = construct_facet(s[0], s[1], s[2], 0, 0);
	for (i = 2 + l; i < r; i++) {	
		if (!selected[i - l] && !on_facet(points, &facet, i)) {
			selected[i - l] = 1;
			s[3] = i;
			break;
		}
	}
	
	// Construct Initial Tetrahedron
	for (i = 0; i < 4; i++) {
		int a = s[i];
		int b = s[(i + 1) % 4];
		int c = s[(i + 2) % 4];
		int d = s[(i + 3) % 4];
		edit_facet(facets + *cnt, a, b, c, *cnt, 1);
		if (above_facet(points, facets + *cnt, d)) {
			edit_facet(facets + *cnt, a, c, b, *cnt, 1);
		}
		construct_edges(points, edges, facets + *cnt);
		(*cnt)++;
	}
	
	// Add New Nodes
	for (i = l + 2; i < r; i++) {
		if (!selected[i - l]) {
			for (j = startCnt; j < *cnt; j++) {
				if ((facets[j].valid) && (above_facet(points, facets + j, i))) {
					change_facet(points, edges, facets, facets + j, cnt, i);
					break;
				}
			}
		}
	}
	
	// Post Cleanup
	int index = startCnt;
	int *proj = (int *) malloc(sizeof(int) * (*cnt - startCnt + 1));
	
	for (i = startCnt; i < *cnt; i++) {
		while (facets[index].valid) {
			index++;
		}
		if ((facets[i].valid) && (i > index)) {
			copy_facet(facets + index, facets + i);
			facets[index].id = i;
			facets[i].valid = 0;
			proj[i - startCnt] = index;
		} else {
			proj[i - startCnt] = i;
			if (!facets[i].valid) {
				proj[i - startCnt] = -1;
			}
		}
	}

	for (i = l; i < r; i++) {
		struct Edge *edge = edges[i].next;
		struct Edge *pre = edges + i;
		while (edge != NULL) {
			if (proj[edge->id - startCnt] == -1) {
				pre->next = edge->next;
				free(edge);
				edge = pre->next;
			} else {
				edge->id = proj[edge->id - startCnt];
				pre = edge;
				edge = edge->next;
			}
		}
	}

	free(proj);
	free(selected);
	
	if (facets[index].valid) {
		*cnt = index + 1;
	} else {
		*cnt = index;
	}
}


// Construct the 2d convex hull according to the projection
int* convex_hull_2d(const struct Point *points, const struct Edge *edges, const struct Facet *facets, int fl, int fr) {
	int i;
	// Find Point with Smallest Y
	int minNode = -1;

	for (i = fl; i < fr; i++) {
		if (facets[i].valid) {
			if ((minNode == -1) || (points[facets[i].a].y < points[minNode].y)) {
				minNode = facets[i].a;
			}
			if (points[facets[i].b].y < points[minNode].y) {
				minNode = facets[i].b;
			}
			if (points[facets[i].c].y < points[minNode].y) {
				minNode = facets[i].c;
			}
		}
	}

	// Find the Second Point on CH
	int *ch = malloc(sizeof(int) * (fr - fl + 200));
	struct Point pseudo;
	pseudo.x = points[minNode].x - 1;
	pseudo.y = points[minNode].y;
	double maxDotProd = -10.0;
	int maxNode = -1;
	const struct Edge *edge = edges[minNode].next;
	
	while (edge != NULL) {
		double dotProd = normal_dot_prod_2d(&pseudo, points + minNode, points + edge->dst);
		if (dotProd > maxDotProd) {
			maxDotProd = dotProd;
			maxNode = edge->dst;
		}
		edge = edge->next;
	}
	
	// Construct 2D Convex Hull
	ch[0] = minNode;
	ch[1] = maxNode;
	int top = 1;
	while (ch[top] != minNode) {
		edge = edges[ch[top]].next;
		maxNode = -1;
		maxDotProd = -10.0;
		while (edge != NULL) {
			double dotProd = normal_dot_prod_2d(points + ch[top - 1], points + ch[top], points + edge->dst);
			if (dotProd > maxDotProd + EPS) {
				maxDotProd = dotProd;
				maxNode = edge->dst;
			}
			edge = edge->next;
		}
		ch[++top] = maxNode;
	}
	ch[++top] = minNode;
	return ch;
}

// Find the most left point of right convex hull
int find_left_most(const struct Point *points, int *ch) {
	int i = 0;
	int mostNode = 0;
	do {
		if (points[ch[i]].x < points[ch[mostNode]].x) {
			mostNode = i;
		}
		i++;
	} while (ch[i] != ch[0]);
	return mostNode;
}

// Find the most right point of left convex hull
int find_right_most(const struct Point *points, int *ch) {
	int i = 0;
	int mostNode = 0;
	do {
		if (points[ch[i]].x > points[ch[mostNode]].x) {
			mostNode = i;
		}
		i++;
	} while (ch[i] != ch[0]);
	return mostNode;
}

// Try to rotate the pivot on the left convex hull clockwisely
int clockwise_rotate(const struct Point *points, int *ch, int *pivot, int index) {
	int next = *pivot - 1;
	if (next < 0) {
		next = 0;
		while (ch[next + 1] != ch[0]) {
			next++;
		}
	}
	if (cross_prod_2d(points + ch[next], points + ch[*pivot], points + index) <= 0) {
		*pivot = next;
		return 1;
	}
	return 0;
}

// Try to rotate the pivot on the right convex hull counter-clockwisely
int counter_clockwise_rotate(const struct Point *points, int *ch, int *pivot, int index) {
	int next = *pivot + 1;
	if (ch[next] == ch[0]) {
		next = 0;
	}
	if (cross_prod_2d(points + ch[next], points + ch[*pivot], points + index) >= 0) {
		*pivot = next;
		return 1;
	}
	return 0;
}

// Find the next point of the common fillet
int find_best_fillet(const struct Point *point, const struct Edge *edges, const struct Facet *facets, const struct Point *normal, int index1, int index2) {
	double maxDotProd = -10.0;
	double minLength = -10.0;
	int maxIndex = -1e8;
	const struct Edge *edge = edges[index1].next;
	
	while (edge != NULL) {
		struct Point tem_normal = tem_facet_normal(point + index2, point + index1, point + edge->dst);
		double diff = normal_dot_prod(normal, &tem_normal) - maxDotProd;
		if (diff > EPS) {
			maxDotProd = normal_dot_prod(normal, &tem_normal);
			minLength = dist(tem_normal);
			maxIndex = edge->dst;
		} else if ((diff > - EPS) && (dist(tem_normal) < minLength)) {
			minLength = dist(tem_normal);
			maxIndex = edge->dst;
		}
		edge = edge->next;
	}
	edge = edges[index2].next;
	while (edge != NULL) {
		struct Point tem_normal = tem_facet_normal(point + index2, point + index1, point + edge->dst);
		double diff = normal_dot_prod(normal, &tem_normal) - maxDotProd;
		if (diff > EPS) {
			maxDotProd = normal_dot_prod(normal, &tem_normal);
			minLength = dist(tem_normal);
			maxIndex = - edge->dst;
		} else if ((diff > - EPS) && (dist(tem_normal) < minLength)) {
			minLength = dist(tem_normal);
			maxIndex = - edge->dst;
		}
		edge = edge->next;
	}
	return maxIndex;
}

// Use Flood-Fill to Delete Internal Facets
void delete_facets(const struct Point *points, const struct Edge *edges, struct Facet *facets, int startID, int size) {
	size = size * 2;
	int *queue = (int *) malloc(sizeof(int) * size);
	int l = 0, r = 1;
	queue[0] = startID;
	while (l != r) {
		int curID = queue[l];
		l = (l + 1) % size;
		if (facets[curID].valid) {
			facets[curID].valid = 0;
			if (!find_edge_bound(edges + facets[curID].a, facets[curID].b)) {
				int id = find_edge_id(edges + facets[curID].b, facets[curID].a);
				if (facets[id].valid) {
					queue[r] = id;
					r = (r + 1) % size;
				}
			} else {
				set_edge_bound(edges + facets[curID].a, facets[curID].b, 0);
			}
			if (!find_edge_bound(edges + facets[curID].b, facets[curID].c)) {
				int id = find_edge_id(edges + facets[curID].c, facets[curID].b);
				if (facets[id].valid) {
					queue[r] = id;
					r = (r + 1) % size;
				}
			} else {
				set_edge_bound(edges + facets[curID].b, facets[curID].c, 0);
			}
			if (!find_edge_bound(edges + facets[curID].c, facets[curID].a)) {
				int id = find_edge_id(edges + facets[curID].a, facets[curID].c);
				if (facets[id].valid) {
					queue[r] = id;
					r = (r + 1) % size;
				}
			} else {
				set_edge_bound(edges + facets[curID].c, facets[curID].a, 0);
			}
		} 
	}
	free(queue);
}

// Merge Two Convex Hull
int merge(const struct Point *points, struct Edge *edges, struct Facet *facets, int fl1, int fr1, int fl2, int fr2) {
	int i, j, cnt;
	
	// Pre Cleanup
	int minNodeIndex = 1e8;
	int maxNodeIndex = -1;
	for (i = fl2; i < fr2; i++) {
		minNodeIndex = MIN(minNodeIndex, facets[i].a);
		minNodeIndex = MIN(minNodeIndex, facets[i].b);
		minNodeIndex = MIN(minNodeIndex, facets[i].c);
		maxNodeIndex = MAX(maxNodeIndex, facets[i].a);
		maxNodeIndex = MAX(maxNodeIndex, facets[i].b);
		maxNodeIndex = MAX(maxNodeIndex, facets[i].c);
		copy_facet(facets + (i - fl2 + fr1), facets + i);
		facets[i - fl2 + fr1].id = i - fl2 + fr1;
		if (i != i - fl2 + fr1) {
			facets[i].valid = 0;
		}
	}
	
	for (i = minNodeIndex; i <= maxNodeIndex; i++) {
		struct Edge *edge = edges[i].next;
		while (edge != NULL) {
			edge->id -= fl2 - fr1;
			edge = edge->next;
		}
	}
	cnt = fr1 + (fr2 - fl2);
	fl2 = fr1;
	fr2 = cnt;
	
	
	for (i = fl1; i < fr1; i++) {
		minNodeIndex = MIN(minNodeIndex, facets[i].a);
		minNodeIndex = MIN(minNodeIndex, facets[i].b);
		minNodeIndex = MIN(minNodeIndex, facets[i].c);
	}
	
	// Find New Edge (Project Points to XY Plane)
	
	int *ch1 = convex_hull_2d(points, edges, facets, fl1, fr1);
	int *ch2 = convex_hull_2d(points, edges, facets, fl2, fr2);
	
	int pivot1 = find_right_most(points, ch1);
	int pivot2 = find_left_most(points, ch2);
	
	while (clockwise_rotate(points, ch1, &pivot1, ch2[pivot2]) || counter_clockwise_rotate(points, ch2, &pivot2, ch1[pivot1])) {
	}
	
	// Find fillets connecting two convex hull
	int suppNode1 = ch1[pivot1];
	int suppNode2 = ch2[pivot2];
	pivot1 = ch1[pivot1];
	pivot2 = ch2[pivot2];
	
	struct Point pseudo;
	pseudo.x = (points[pivot1].x + points[pivot2].x) / 2.0;
	pseudo.y = (points[pivot1].y + points[pivot2].y) / 2.0;
	pseudo.z = points[pivot1].z + points[pivot2].z;

	int start1 = -1;
	int start2 = -1;
	
	int *filletList = malloc(sizeof(int) * cnt * 2);
	int filletCount = 0;
	
	do {
		struct Point normal = tem_facet_normal(points + pivot1, points + pivot2, &pseudo);
		normal = normalize(&normal);
		int newIndex = find_best_fillet(points, edges, facets, &normal, pivot1, pivot2);
		filletList[filletCount++] = newIndex;
		
		if (newIndex >= 0) {
			copy_point(&pseudo, points + pivot1);
			start1 = set_edge_bound(edges + newIndex, pivot1, 1);
			pivot1 = newIndex;
		} else {
			copy_point(&pseudo, points + pivot2);
			start2 = set_edge_bound(edges + pivot2, - newIndex, 1);
			pivot2 = - newIndex;
		}

	} while ((pivot1 != suppNode1) || (pivot2 != suppNode2));
	
	// Use flood-fill to delete useless facets in between
	delete_facets(points, edges, facets, start1, fr1 - fl1);
	delete_facets(points, edges, facets, start2, fr2 - fl2);
	
	// Construct new convex hull
	for (i = 0; i < filletCount; i++) {
		int newIndex = filletList[i];
		if (newIndex >= 0) {
			edit_edge_id(edges + newIndex, pivot1, cnt);
			construct_edge(edges + pivot1, pivot2, cnt);
			construct_edge(edges + pivot2, newIndex, cnt);
			facets[cnt] = construct_facet(pivot1, pivot2, newIndex, cnt, 1);
			cnt++;
			pivot1 = newIndex;
		} else {
			newIndex = - newIndex;
			edit_edge_id(edges + pivot2, newIndex, cnt);
			construct_edge(edges + pivot1, pivot2, cnt);
			construct_edge(edges + newIndex, pivot1, cnt);
			facets[cnt] = construct_facet(pivot1, pivot2, newIndex, cnt, 1);
			cnt++;		
			pivot2 = newIndex;
		}
	}	
	
	free(filletList);
	free(ch1);
	free(ch2);
	
	// Post Cleanup
	
	int startCnt = fl1;
	int index = fl1;
	int *proj = (int *) malloc(sizeof(int) * (cnt - fl1 + 1));
	
	for (i = startCnt; i < cnt; i++) {
		while (facets[index].valid) {
			index++;
		}
		if ((facets[i].valid) && (i > index)) {
			copy_facet(facets + index, facets + i);
			facets[index].id = i;
			facets[i].valid = 0;
			proj[i - startCnt] = index;
		} else {
			proj[i - startCnt] = i;
			if (!facets[i].valid) {
				proj[i - startCnt] = -1;
			}
		}
	}

	for (i = minNodeIndex; i <= maxNodeIndex; i++) {
		struct Edge *edge = edges[i].next;
		struct Edge *pre = edges + i;
		while (edge != NULL) {
			if (proj[edge->id - startCnt] == -1) {
				pre->next = edge->next;
				free(edge);
				edge = pre->next;
			} else {
				edge->id = proj[edge->id - startCnt];
				pre = edge;
				edge = edge->next;
			}
		}
	}

	free(proj);
	if (facets[index].valid) {
		cnt = index + 1;
	} else {
		cnt = index;
	}
	
	return cnt;
}

// Construct the convex hull and return the maximum id of valid facets
int split_and_conquer(struct Point *points, struct Edge *edges, struct Facet *facets, int l, int r) {
	if (r - l < 200) {
		int facetCount = l * 6;
		construct(points, edges, facets, l, r, &facetCount);
		return facetCount;
	}
	int mid = (l + r) / 2;
	int facetCount1 = split_and_conquer(points, edges, facets, l, mid);
	int facetCount2 = split_and_conquer(points, edges, facets, mid, r);
	return merge(points, edges, facets, l * 6, facetCount1, mid * 6, facetCount2);
}

int main(int argc, char *argv[]) {

	MPI_Status status;
	int i, j, N, shape, P, p, rc, I, R, L, n;
	
	// Set Random Seed
	srand(1000);
	
	// Get Parameters
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
	
	// Initial the Positions of Points
	struct Point *points = (struct Point *) malloc(sizeof(struct Point) * (N + 1));
	generate(points, N, shape);
	qsort(points, N, sizeof(struct Point), compare_point);
	
	if ((p == 0) && (N < 3000)) {
		print_positions(points, N);
	}
	
	struct Edge *edges = (struct Edge *) malloc(sizeof(struct Edge) * (N + 1));
	for (i = 0; i < N; i++) {
		edges[i].next = NULL;
		edges[i].bound = 0;
	}
	struct Facet *facets = (struct Facet *) malloc(sizeof(struct Facet) * (N * 6 + 1));
	for (i = 0; i < 6 * N + 1; i++) {
		facets[i].valid = 0;
	}
	
	// Linear Load-balanc Data Distribution
	I = (N + P - p - 1) / P;
	R = N % P;
	L = N / P;
	n = p * L + MIN(p, R);
	
	// Algorithm Starts Here
	MPI_Barrier(MPI_COMM_WORLD);
	double startTime = MPI_Wtime();
	
	int cnt = split_and_conquer(points, edges, facets, n, n + I);
	
	int round = ceil(log(P) / log(2));
	int flip = 1;
	for (i = 0; i < round; i++) {
		if (((p ^ flip) < P) && ((p % flip) == 0)) {
			int dst = (p ^ flip);
			if (dst > p) {
				int start, end;
				MPI_Recv(&start, 1, MPI_INT, dst, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(&end, 1, MPI_INT, dst, 2, MPI_COMM_WORLD, &status);
				int *nodeA = (int *) malloc(sizeof(int) * (end - start));
				int *nodeB = (int *) malloc(sizeof(int) * (end - start));
				int *nodeC = (int *) malloc(sizeof(int) * (end - start));
				int *valid = (int *) malloc(sizeof(int) * (end - start));
				
				MPI_Recv(nodeA, end - start, MPI_INT, dst, 3, MPI_COMM_WORLD, &status);
				MPI_Recv(nodeB, end - start, MPI_INT, dst, 4, MPI_COMM_WORLD, &status);
				MPI_Recv(nodeC, end - start, MPI_INT, dst, 5, MPI_COMM_WORLD, &status);
				MPI_Recv(valid, end - start, MPI_INT, dst, 6, MPI_COMM_WORLD, &status);
				
				// Update edges received from other processors
				for (j = start; j < end; j++) {
					edit_facet(facets + j, nodeA[j - start], nodeB[j - start], nodeC[j - start], j, valid[j - start]);
					construct_edges(points, edges, facets + j);
				}
				
				// Merge the receivde convex hull with the current one
				cnt = merge(points, edges, facets, n * 6, cnt, start, end);

				free(nodeA);
				free(nodeB);
				free(nodeC);
				free(valid);
			} else { // send
				int start = 6 * n;
				MPI_Send(&start, 1, MPI_INT, dst, 1, MPI_COMM_WORLD);
				MPI_Send(&cnt, 1, MPI_INT, dst, 2, MPI_COMM_WORLD);

				int *nodeA = (int *) malloc(sizeof(int) * (cnt - start));
				int *nodeB = (int *) malloc(sizeof(int) * (cnt - start));
				int *nodeC = (int *) malloc(sizeof(int) * (cnt - start));
				int *valid = (int *) malloc(sizeof(int) * (cnt - start));
				for (j = start; j < cnt; j++) {
					nodeA[j - start] = facets[j].a;
					nodeB[j - start] = facets[j].b;
					nodeC[j - start] = facets[j].c;
					valid[j - start] = facets[j].valid;
				}
				// Send information of current convex hull to other processor
				MPI_Send(nodeA, cnt - start, MPI_INT, dst, 3, MPI_COMM_WORLD);
				MPI_Send(nodeB, cnt - start, MPI_INT, dst, 4, MPI_COMM_WORLD);
				MPI_Send(nodeC, cnt - start, MPI_INT, dst, 5, MPI_COMM_WORLD);
				MPI_Send(valid, cnt - start, MPI_INT, dst, 6, MPI_COMM_WORLD);

				free(nodeA);
				free(nodeB);
				free(nodeC);
				free(valid);
			}
		}
		flip += flip;
	}
	
	// Algorithm ends here
	MPI_Barrier(MPI_COMM_WORLD);
	if (p == 0) {
		double endTime = MPI_Wtime();
		printf("Execution Time: %e\n", endTime - startTime);
		if ((p == 0) && (N < 3000)) {
			print_facets_write(points, facets, 0, cnt);
		}
	}
	
	// Avoid Memory Leak
	for (i = 0; i < N; i++) {
		struct Edge *edge = edges[i].next;
		while (edge != NULL) {
			edges[i].next = edge->next;
			free(edge);
			edge = edges[i].next;
		}
	}

	free(points);
	free(edges);
	free(facets);
	
	rc = MPI_Finalize();
	return 0;
}

