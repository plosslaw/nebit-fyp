#include <algorithm>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <set>
#include <string>
#include <vector>

struct Graph {
    std::vector<int> vertices;
    std::vector<std::set<int>> edges;
};

typedef std::vector<double> row;
typedef std::vector<row> matrix;

int parseProgramMode();
void printProgressBar();
matrix parseUserInputMatrix();
int parseInteger();
double parseDouble();
std::vector<double> parseMatrixRows(int row_number, int size);

bool verifyQuasiBistochasticMatrix(matrix m);

void printNegativityResults(double negativity, std::vector<std::pair<double, std::string>> min_negativity_decomposition);
void printMostNegativeSingleQubitRotation();

std::vector<double> get_minimal_negativity_birkhoff_decomposition(matrix m);