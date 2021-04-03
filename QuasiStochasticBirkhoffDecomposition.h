#include <algorithm>
#include <iostream>
#include <limits>
#include <math.h>
#include <set>
#include <string>
#include <vector>

struct Graph {
    std::vector<int> vertices;
    std::vector<std::set<int>> edges;
};

typedef std::vector<double> column;
typedef std::vector<column> matrix;


std::vector<double> get_minimal_negativity_birkhoff_decomposition(matrix m);