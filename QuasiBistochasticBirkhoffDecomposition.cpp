#include "QuasiBistochasticBirkhoffDecomposition.h"
#include <algorithm>
#include <sstream>
#include <string>

using namespace std;


const double EPS = 1e-6;
const double INF = 1e7;

int main() {
    // printProgressBar();
    matrix user_matrix;

    int option = parseProgramMode();
    pair<double, vector<pair<double, vector<int>>>> min_negativity_decomposition;
    switch (option) {
    case 1:
        cout<<"\nOption 1 chosen\n\n";

        while (true) {
            user_matrix = parseUserInputMatrix();
            cout<<"Checking if input matrix is quasi-bistochastic...\n\n";
            if (verifyQuasiBistochasticMatrix(user_matrix)) {
                cout<<"Input matrix is quasi-bistochastic\n\n";
                break;
            }
            cout<<"Input matrix is not quasi-bistochastic, the sum of each column and sum of each row has to be 1.\nPlease enter a new quasi-bistochastic matrix\n\n";
        }

        cout<<"Computing negativity of minimal negativity Birkhoff decomposition...\n\n";

        min_negativity_decomposition = get_minimal_negativity_birkhoff_decomposition(user_matrix);

        printNegativityResults(min_negativity_decomposition.first, min_negativity_decomposition.second);
        break;
    case 2:
        cout<<"\nOption 2 chosen\n\n";
        cout<<"This feature will be added in future updates when we find a better way to reduce \nthe time complexity of finding the minimal negativity Birkhoff decomposition\n\n";
        break;
    case 3:
        cout<<"\nProgram will now end. Have a nice day.\n\n";
        break;
    default:
        throw "Unexpected path\n";
    }
    


    return 0;
}







int parseProgramMode() {
    string welcome_message = "\nWhat would you like to do?\n";
    string first_option_message = 
        "\t1. Compute the negativity of quasi-bistochastic matrix using \n\tminimal negativity Birkhoff-von Neumann Decomposition\n\n";
    string second_option_message = 
        "\t2. Compute the rotation matrix for the most negative single qubit unitary \n\toperation using minimal negativity Birkhoff-von Neumann Decomposition\n\n";
    string third_option_message = "\t3. Exit program\n\n";
    string start_message = welcome_message + first_option_message + second_option_message + third_option_message;
    cout<<start_message;

    int input;
    while (true) {
        input = parseInteger();
        if (input == 1 || input == 2 || input == 3) {
            break;
        } else {
            cout<<"Invalid option, please choose one of the options listed\n";
        }
    }

    return input;
}

void printProgressBar() {
    float progress = 0.0;
    int barWidth = 70;
    while (progress <= 1.0) {

        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();

        progress += 0.01; // for demonstration only
    }
    std::cout << "[";
    int pos = barWidth;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << "DONE\r";
    std::cout.flush();

    std::cout << std::endl;
}

matrix parseUserInputMatrix() {
    matrix m;
    cout<<"The quasi-bistochastic matrix has to be a square matrix.\nPlease enter the size of quasi-bistochastic matrix\n";
    cout<<"e.g. Enter 4 for 4x4 matrix\n\n";
    int size;

    while (true) {
        size = parseInteger();
        if (size > 0) {
            break;
        } else {
            // cout<<"Please enter an integer greater than 0\n";
        }
    }

    for (int i = 0; i < size; i++) {
        vector<double> v = parseMatrixRows(i, size);
        m.push_back(v);
    }
    cout<<"\n==========Input Matrix==========\n\n";
    for (row ele : m) {
        for (double ele2 : ele) {
            cout<<fixed<<setprecision(3)<<ele2<<"\t";
        }
        cout<<"\n\n";
    }

    return m;
}

int parseInteger() {
    int input;
    while (true) {
        string input_string;
        getline(cin,input_string);
        input_string.erase(remove_if(input_string.begin(), input_string.end(), ::isspace), input_string.end());
        char *p;
        input = strtol(input_string.c_str(), &p, 10);
        if (*p) {
            cout<<"Invalid input, please enter an appropriate integer\n";
        } else {
            break;
        }
    }
    return input;
}

vector<double> parseMatrixRows(int row_number, int size) {
    vector<double> input_row;
    cout<<"Enter row entries (double values) for row "<<row_number<<'\n';

    bool continue_loop = true;
    while (continue_loop) {
        vector<string> temp_vector;
        vector<double> temp_input_row;
        string temp_string;
        string input_row_string;
        getline(cin, input_row_string);
        stringstream ss(input_row_string);
        while (getline(ss, temp_string, ' ')) {
            temp_vector.push_back(temp_string);
        }
        if (size != temp_vector.size()) {
            cout<<"Incorrect number of entries detected, re-enter row entries (double values) for row "<<row_number<<'\n';
            continue;
        }

        for (int i = 0; i < temp_vector.size(); i++) {
            string ele = temp_vector[i];
            char *p;
            double temp_entry_value = strtod(ele.c_str(), &p);
            if (*p) {
                cout<<"Invalid entry value detected, re-enter row entries (double values) for row "<<row_number<<'\n';
                break;
            } else {
                temp_input_row.push_back(temp_entry_value);
            }
            if (i == temp_vector.size() - 1) {
                // last element
                input_row = temp_input_row;
                continue_loop = false;
            }
        }

    }
    return input_row;
}

double parseDouble() {
    double input;
    while (true) {
        string input_string;
        getline(cin,input_string);
        input_string.erase(remove_if(input_string.begin(), input_string.end(), ::isspace), input_string.end());
        char *p;
        input = strtod(input_string.c_str(), &p);
        if (*p) {
            cout<<"Invalid input, please enter an appropriate double value\n";
        } else {
            break;
        }
    }
    return input;
}

bool verifyQuasiBistochasticMatrix(matrix m) {
    for (row ele : m) {
        double sum = 0;
        for (double ele2 : ele) {
            sum += ele2;
        }
        if (sum - 1 >= EPS || 1 - sum >= EPS) {
            return false;
        }
    }

    for (int i = 0; i < m.size(); i++) {
        double sum = 0;
        for (row ele : m) {
            sum += ele[i];
        }

        if (sum - 1 >= EPS || 1 - sum >= EPS) {
            return false;
        }
    }
    return true;
}

void printNegativityResults(double negativity, vector<pair<double, vector<int>>> min_negativity_decomposition) {
    if (min_negativity_decomposition.size() == 0) {
        return;
    }

    vector<pair<double, vector<int>>> positive_part;
    vector<pair<double, vector<int>>> negative_part;
    int size = min_negativity_decomposition[0].second.size();
    for (pair<double, vector<int>> ele : min_negativity_decomposition) {
        if (ele.first < 0) {
            negative_part.push_back(ele);
        }else if (ele.first > 0) {
            positive_part.push_back(ele);
        }
    }
    cout<<"======POSITIVE PART DECOMPOSITION======\n";
    for (pair<double, vector<int>> ele : positive_part) {
        cout<<ele.first<<" * [ ";
        for (int perm_ele : ele.second) {
            cout<<perm_ele<<" ";
        }
        cout<<"]\n";
    }
    cout<<"\n======NEGATIVE PART DECOMPOSITION======\n";
    for (pair<double, vector<int>> ele : negative_part) {
        cout<<ele.first<<" * [ ";
        for (int perm_ele : ele.second) {
            cout<<perm_ele<<" ";
        }
        cout<<"]\n";
    }

    matrix positive_matrix;
    for (int i = 0; i < size; i++) {
        row zero_row(size, 0);
        positive_matrix.push_back(zero_row);
    }
    for (pair<double, vector<int>> ele : positive_part) {
        for (int j = 0; j < ele.second.size(); j ++) {
            int perm_ele = ele.second[j];
            positive_matrix[j][perm_ele] += ele.first;
        }
    }

    matrix negative_matrix;
    for (int i = 0; i < size; i++) {
        row zero_row(size, 0);
        negative_matrix.push_back(zero_row);
    }
    for (pair<double, vector<int>> ele : negative_part) {
        for (int j = 0; j < ele.second.size(); j ++) {
            int perm_ele = ele.second[j];
            negative_matrix[j][perm_ele] += ele.first;
        }
    }



    cout<<"\n======POSITIVE MATRIX PART======\n";
    for (row ele : positive_matrix) {
        for (double ele2 : ele) {
            cout<<fixed<<setprecision(3)<<ele2<<"\t";
        }
        cout<<"\n\n";
    }

    cout<<"\n======NEGATIVE MATRIX PART======\n";
    for (row ele : negative_matrix) {
        for (double ele2 : ele) {
            cout<<fixed<<setprecision(3)<<ele2<<"\t";
        }
        cout<<"\n\n";
    }

    cout<<"The negativity of the minimal negativity Birkhoff Decomposition is \n\n";
    cout<<fixed<<setprecision(5)<<negativity<<"\n\n";


}







pair<double, vector<pair<double, vector<int>>>> get_minimal_negativity_birkhoff_decomposition(matrix m) {
    vector<pair<double, vector<int>>> results;
    double delta = 0;
    for (row ele : m) {
        for (double ele2 : ele) {
            delta = min(delta, ele2);
        }
    }
    delta = abs(delta); //make delta positive
    // cout<<delta<<'\n';

    /**
     * S_ = (1+d*delta) A - d*delta B
     */

    matrix positive_part;
    matrix all_one_matrix_normalized;
    
    int num_of_nonzero_entries_in_positive = 0;

    for (row ele : m) {
        vector<double> temp_row;
        double normalized_ones = 1.0/m.size(); // 1/d
        double positive_norm_constant = 1.0/(1.0 + m.size()*delta); // 1/(1+d*delta)
        vector<double> temp_ones_row(m.size(), normalized_ones);
        for (double ele2 : ele) {
            if ((ele2 + delta)*positive_norm_constant < EPS) {
                temp_row.push_back(0);
                // cout<<"0\t";
            } else {
                temp_row.push_back((ele2 + delta)*positive_norm_constant);
                num_of_nonzero_entries_in_positive++;
                // cout<<(ele2 + delta)*positive_norm_constant<<"\t";
            }
        }
        // cout<<'\n';
        positive_part.push_back(temp_row);
        all_one_matrix_normalized.push_back(temp_ones_row);
    }
    vector<vector<int>> all_permutation_matrices;
    vector<int> base_vector;
    for (int i = 0; i < m.size(); i++) {
        base_vector.push_back(i);
    }
    long long num_of_permutation_matrices = factorial(m.size());
    for (long long i = 0; i<num_of_permutation_matrices; i++) {
        vector<int> v = base_vector;
        all_permutation_matrices.push_back(v);
        next_permutation(base_vector.begin(), base_vector.end());
    }

    // for (vector<int> ele_vec : all_permutation_matrices) {
    //     for (int ele : ele_vec) {
    //         cout<<ele<<" ";
    //     }
    //     cout<<'\n';
    // }

    vector<vector<pair<double, vector<int>>>> positive_part_normalized_decompositions;
    vector<pair<double, vector<int>>> positive_part_intermediate_results;
    DFS_all_birkhoff_decomposition(positive_part_normalized_decompositions, positive_part, all_permutation_matrices, positive_part_intermediate_results, num_of_nonzero_entries_in_positive);
    vector<vector<pair<double, vector<int>>>> all_ones_normalized_decompositions;
    // vector<pair<double, vector<int>>> sample_decomposition;
    // sample_decomposition.push_back({0.25, vector<int>{0,1,2,3}});
    // sample_decomposition.push_back({0.25, vector<int>{1,0,3,2}});
    // sample_decomposition.push_back({0.25, vector<int>{3,2,1,0}});
    // sample_decomposition.push_back({0.25, vector<int>{2,3,0,1}});
    // all_ones_normalized_decompositions.push_back(sample_decomposition);
    vector<pair<double, vector<int>>> all_ones_intermediate_results;
    DFS_all_birkhoff_decomposition(all_ones_normalized_decompositions, all_one_matrix_normalized, all_permutation_matrices, all_ones_intermediate_results, m.size() * m.size());


    vector<vector<pair<double, vector<int>>>> positive_part_normalized_decompositions_condensed;
    vector<vector<pair<double, vector<int>>>> all_ones_normalized_decompositions_condensed;

    set<string> unique_set_ones;
    for (vector<pair<double, vector<int>>> all_ones_decomposition : all_ones_normalized_decompositions) {
        vector<string> temp_vector;
        for (pair<double, vector<int>> all_ones_term : all_ones_decomposition) {
            string key = to_string(all_ones_term.first); // add coefficient of permutation matrix to key
            if (key.size() > 7) {
                key = key.substr(0,7);
            }
            for (int ele : all_ones_term.second) {
                key += to_string(ele) + ",";
            }
            temp_vector.push_back(key);
        }
        sort(temp_vector.begin(), temp_vector.end());
        string set_key = "";
        for (string ele : temp_vector) {
            set_key += ele + ",";
        }
        if (unique_set_ones.find(set_key) == unique_set_ones.end()) {
            unique_set_ones.insert(set_key);
            all_ones_normalized_decompositions_condensed.push_back(all_ones_decomposition);
        }
    }


    set<string> unique_set_positive_part;
    for (vector<pair<double, vector<int>>> positive_part_decomposition : positive_part_normalized_decompositions) {
        vector<string> temp_vector;
        for (pair<double, vector<int>> positive_part_term : positive_part_decomposition) {
            string key = to_string(positive_part_term.first); // add coefficient of permutation matrix to key
            if (key.size() > 7) {
                key = key.substr(0,7);
            }
            for (int ele : positive_part_term.second) {
                key += to_string(ele) + ",";
            }
            temp_vector.push_back(key);
        }
        sort(temp_vector.begin(), temp_vector.end());
        string set_key = "";
        for (string ele : temp_vector) {
            set_key += ele + ",";
        }
        if (unique_set_positive_part.find(set_key) == unique_set_positive_part.end()) {
            unique_set_positive_part.insert(set_key);
            positive_part_normalized_decompositions_condensed.push_back(positive_part_decomposition);
        }
    }


    // for (vector<pair<double, vector<int>>> positive_part_decomposition : positive_part_normalized_decompositions) {
    //     cout<<"===POSSIBLE DECOMPOSITION===\n\n";
    //     for (pair<double, vector<int>> positive_part_term : positive_part_decomposition) {
    //         cout<<positive_part_term.first<<" [ ";
    //         for (int perm_matrix_ele : positive_part_term.second) {
    //             cout<<perm_matrix_ele<<" ";
    //         }
    //         cout<<"] ";
    //     }
    //     cout<<'\n';
    // }

    long long chosen_positive_decomposition;
    long long chosen_ones_decomposition;
    double minimal_negativity = INF;
    for (int i = 0; i < positive_part_normalized_decompositions_condensed.size(); i++) {
        for (int j = 0; j < all_ones_normalized_decompositions_condensed.size(); j++) {
            vector<pair<double, vector<int>>> positive_part_decomposition = positive_part_normalized_decompositions_condensed[i];
            vector<pair<double, vector<int>>> all_ones_decomposition = all_ones_normalized_decompositions_condensed[j];
            map<string, double> summing_up_coefficients;

            double temp_negativity = 0;
            for (pair<double, vector<int>> positive_part_term : positive_part_decomposition) {
                string map_key = "";
                for (int ele : positive_part_term.second) {
                    map_key += to_string(ele)+",";
                }
                summing_up_coefficients.insert({map_key, (1 + m.size() * delta) * positive_part_term.first});
                // summing_up_coefficients[map_key] = (1 + m.size() * delta) * positive_part_term.first;
            }
            for (pair<double, vector<int>> all_ones_term : all_ones_decomposition) {
                string map_key = "";
                for (int ele : all_ones_term.second) {
                    map_key += to_string(ele)+",";
                }
                if (summing_up_coefficients.find(map_key) == summing_up_coefficients.end()) {
                    temp_negativity += (m.size() * delta) * all_ones_term.first;
                } else {
                    double temp_coefficient = summing_up_coefficients[map_key] - (m.size() * delta) * all_ones_term.first;
                    if (temp_coefficient < 0) {
                        temp_negativity += abs(temp_coefficient);
                    }
                }
            }
            if (minimal_negativity > temp_negativity) {
                minimal_negativity = temp_negativity;
                chosen_positive_decomposition = i;
                chosen_ones_decomposition = j;
            }
        }
    }
    // cout<<minimal_negativity<<'\n';
    map<string, double> summing_up_coefficients;

    for (pair<double, vector<int>> positive_part_term : positive_part_normalized_decompositions_condensed[chosen_positive_decomposition]) {
        string map_key = "";
        for (int ele : positive_part_term.second) {
            map_key += to_string(ele) + ",";
        }
        summing_up_coefficients.insert({map_key, (1 + m.size() * delta) * positive_part_term.first});
        // summing_up_coefficients[map_key] = (1 + m.size() * delta) * positive_part_term.first;
    }
    for (pair<double, vector<int>> all_ones_term : all_ones_normalized_decompositions_condensed[chosen_ones_decomposition]) {
        string map_key = "";
        for (int ele : all_ones_term.second) {
            map_key += to_string(ele) + ",";
        }
        if (summing_up_coefficients.find(map_key) == summing_up_coefficients.end()) {
            results.push_back({(m.size() * delta * -1) * all_ones_term.first, all_ones_term.second});
        } else {
            summing_up_coefficients[map_key] -= (m.size() * delta) * all_ones_term.first;

        }
    }

    for (pair<string, double> ele : summing_up_coefficients) {
        results.push_back({ele.second, convertCSVToVector(ele.first)});
    }

    return {minimal_negativity, results};
}

void DFS_all_birkhoff_decomposition(vector<vector<pair<double, vector<int>>>> &results, matrix &m, 
    vector<vector<int>> &all_permutation_matrices, vector<pair<double, vector<int>>> &intermediate_results, int remaining_entries) {
    if (remaining_entries == 0) {
        vector<pair<double, vector<int>>>copy_intermediate_results = intermediate_results;
        results.push_back(copy_intermediate_results);
        return;
    }

    vector<vector<int>>all_permutation_matrices_copy = all_permutation_matrices;
    for (int i = 0; i < all_permutation_matrices.size(); i++) {
        vector<int> permutation = all_permutation_matrices[i];
        double min_ele = m[0][permutation[0]];
        int new_zero_entries = 0;
        matrix copy_m;
        vector<pair<double, vector<int>>> copy_intermediate_results;
        bool skip = false;
        for (int j = 0; j < permutation.size(); j ++) {
            if (m[j][permutation[j]] < EPS) {
                m[j][permutation[j]] = 0;
                skip = true;
                all_permutation_matrices_copy.erase(all_permutation_matrices_copy.begin()+i);
                break;
            }
            min_ele = min(min_ele, m[j][permutation[j]]);
        }
        if (skip) continue; 

        copy_m = m;
        copy_intermediate_results = intermediate_results;
        for (int j = 0; j < permutation.size(); j++) {
            copy_m[j][permutation[j]] -= min_ele;
            if (copy_m[j][permutation[j]] < EPS) {
                copy_m[j][permutation[j]] = 0;
                new_zero_entries++;
            }
        }
        copy_intermediate_results.push_back({min_ele, permutation});
        DFS_all_birkhoff_decomposition(results, copy_m, all_permutation_matrices_copy, copy_intermediate_results, remaining_entries - new_zero_entries);
    }
}

vector<int> convertCSVToVector(string csv_string) {
    vector<int> results_vector;
    string temp_string;
    string csv_string_copy = csv_string;
    if (csv_string.at(csv_string.size()-1) == ',') {
        // csv_string_trimmed = csv_string.substr(0, csv_string.size()-1);
        csv_string_copy.pop_back();
    }
    stringstream ss(csv_string_copy);
    while (getline(ss, temp_string, ',')) {
        results_vector.push_back(stoi(temp_string));
    }
    return results_vector;
}


long long factorial(int n) {
    if (n == 1 || n == 0) {
        return 1;
    }
    return n * factorial(n-1);
}