
/*
 *
 */
void lin_kernighan(int size, double* x_of_cities, double* y_of_cities, int* opt_tour, double* opt_tour_len);

/*
 * swap two int
 */
void swap(int* a, int* b);


/*
 * generate a random array of 0..size-1
 */
void gen_rand_tour(int size, int* tour);

/*
 * compute the dist of two cities
 */
inline double dist(int x, int y, double* x_of_cities, double* y_of_cities);

/*
 * compute the length of a tour
 */
double tour_len(int size, int* tour, double* x_of_cities, double* y_of_cities);

/*
 *
 */
void prepare(int size, double* x_of_cities, double* y_of_cities, int** near);

/*
 * shift x positions
 */
void shift(int size, int x, int** tour);

/*
 * reverse a array from point a to point b - 1
 */
void reverse(int* a, int* b);

/*
 * generate a array of city y, in an increasing order of the dist between x and y
 */
int* gen_inc_array(int x, double* x_of_cities, double* y_of_cities);
