#include <algorithm>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <map>
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
std::vector<int> convertCSVToVector(std::string csv_string);

void printNegativityResults(double negativity, std::vector<std::pair<double, std::vector<int>>> min_negativity_decomposition);
void printMostNegativeSingleQubitRotation();

std::pair<double, std::vector<std::pair<double, std::vector<int>>>> get_minimal_negativity_birkhoff_decomposition(matrix m);
void DFS_all_birkhoff_decomposition(std::vector<std::vector<std::pair<double, std::vector<int>>>> &results, matrix &m, 
    std::vector<int> permutation, std::set<int> permutation_set, std::vector<std::pair<double, std::vector<int>>> &intermediate_results, int remaining_entries);