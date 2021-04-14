#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

class BipartiteGraph {
    // m is the number of vertices on the left side (number of columns)
    // n is the number of vertices on the right side (number of rows)
    int m, n;

    std::vector<std::vector<int>> edges;
    std::vector<int> pairM, pairN, distance;

public:
    BipartiteGraph(int m, int n);

    // Add edges from left side to right side
    // u is vertex from left side, v is vertex from right side
    void addEdge(int u, int v);

    bool bfs();
    bool dfs(int u);

    int hopcroftKarp();

    std::vector<int> getPermutationMatrix();
};

typedef std::vector<double> row;
typedef std::vector<row> matrix;

int parseProgramMode();
void printProgressBar(double progress);
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
    std::vector<std::vector<int>> &all_permutation_matrices, std::vector<std::pair<double, std::vector<int>>> &intermediate_results, int remaining_entries);


long long factorial(int n);

std::pair<double, std::vector<std::pair<double, std::vector<int>>>> performRandomTrialsOfBirkhoffDecompositions(int num_of_trials, matrix user_matrix);
void DFS_random_birkhoff_decomposition(std::vector<std::pair<double, std::vector<int>>> &results, matrix &m, int remaining_entries, std::vector<std::vector<int>> &preferred_decomposition, int pref_counter, bool &failed);
