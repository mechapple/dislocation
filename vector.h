typedef struct{ double i,j,k; } Vector;
Vector crossProduct(Vector a,Vector b) { Vector c = {a.j*b.k - a.k*b.j, a.k*b.i - a.i*b.k, a.i*b.j - a.j*b.i};	return c; }
double dotProduct(Vector a, Vector b) {	return a.i*b.i+a.j*b.j+a.k*b.k; }
void printVector(Vector a) {std::cout << a.i << " " << a.j << " " << a.k << std::endl;}
