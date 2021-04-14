#include "QuasiBistochasticBirkhoffDecomposition.h"

using namespace std;


const double EPS = 1e-7;
const double INF = 1e7;
const double PI = 3.1415926535897932384626433832;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
#define uid uniform_int_distribution<int>

int main() {
    // printProgressBar();
    matrix user_matrix;

    int option = parseProgramMode();
    int num_of_trials;
    double step_interval;
    pair<double, vector<pair<double, vector<int>>>> min_negativity_decomposition;
    pair<double, pair<matrix, matrix>> positive_negative_decomposition;
    pair<double, pair<matrix, matrix>> most_negative_rotation_results;
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
        cout<<"Input the step interval in radians (at least 0.0001) for the search for \nthe most negative single qubit unitary\n\n";
        // cout<<"This feature will be added in future updates when we find a better way to reduce \nthe time complexity of finding the minimal negativity Birkhoff decomposition\n\n";
        while (true) {
            step_interval = parseDouble();
            if (step_interval >= 0.0001) {
                break;
            }
            cout<<"Please enter a positive step interval larger than 0.0001 radians\n";
        }
        cout<<'\n';
        
        most_negative_rotation_results = calculateMostNegativeSingleQubitUnitary(step_interval);

        printMostNegativeSingleQubitRotation(most_negative_rotation_results.first, most_negative_rotation_results.second.first, most_negative_rotation_results.second.second);
        break;
    case 3:
        cout<<"\nOption 3 chosen\n\n";
        cout<<"How many times do you want to perform the randomized search for minimal decomposition \nfor the given quasi-stochastic matrix ?\n\n";
        while (true) {
            num_of_trials = parseInteger();
            if (num_of_trials > 0) {
                break;
            }
            cout<<"Please enter a positive integer for the number of random trials to perform\n";
        }
        cout<<'\n';

        while (true) {
            user_matrix = parseUserInputMatrix();
            cout<<"Checking if input matrix is quasi-bistochastic...\n\n";
            if (verifyQuasiBistochasticMatrix(user_matrix)) {
                cout<<"Input matrix is quasi-bistochastic\n\n";
                break;
            }
            cout<<"Input matrix is not quasi-bistochastic, the sum of each column and sum of each row has to be 1.\nPlease enter a new quasi-bistochastic matrix\n\n";
        }
        cout<<"Computing the minimal negativity of random Birkhoff decompositions out of "<<num_of_trials<<" random trial(s)...\n\n";

        min_negativity_decomposition = performRandomTrialsOfBirkhoffDecompositions(num_of_trials, user_matrix);

        printNegativityResults(min_negativity_decomposition.first, min_negativity_decomposition.second);
        break;

    case 4:
        cout<<"\nOption 4 chosen\n\n";

        while (true) {
            user_matrix = parseUserInputMatrix();
            cout<<"Checking if input matrix is quasi-bistochastic...\n\n";
            if (verifyQuasiBistochasticMatrix(user_matrix)) {
                cout<<"Input matrix is quasi-bistochastic\n\n";
                break;
            }
            cout<<"Input matrix is not quasi-bistochastic, the sum of each column and sum of each row has to be 1.\nPlease enter a new quasi-bistochastic matrix\n\n";
        }
        cout<<"Computing the minimal negativity of quasi-bistochastic matrix using \ngreedy minimal negativity decomposition...\n\n";

        positive_negative_decomposition = performGreedyMinimalNegativityDecomposition(user_matrix);
        printPositiveNegativeResults(positive_negative_decomposition.first, positive_negative_decomposition.second);
        cout<<"\nVerifying decomposition...\n\n";
        if (verifyDecomposition(user_matrix, positive_negative_decomposition.second.first, positive_negative_decomposition.second.second)) {
            cout<<"Decomposition is valid!\n\n";
        } else {
            cout<<"Decomposition is invalid!\n\n";
        }
        
        break;

    case 5:
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
    string third_option_message = 
        "\t3. Approximate the minimal negativity of quasi-bistochastic matrix by performing \n\tsome number of random trials\n\n";
    string fourth_option_message = 
        "\t4. Compute minimal negativity of quasi-bistochastic matrix using greedy \n\tminimal negativity decomposition\n\n";

    string exit_option_message = "\t5. Exit program\n\n";
    string start_message = welcome_message + first_option_message + second_option_message +
        third_option_message + fourth_option_message + exit_option_message;
    cout<<start_message;

    int input;
    while (true) {
        input = parseInteger();
        if (input == 1 || input == 2 || input == 3 || input == 4 || input == 5) {
            break;
        } else {
            cout<<"Invalid option, please choose one of the options listed\n";
        }
    }

    return input;
}

void printProgressBar(double progress) {
    // double progress = 0.0;
    int barWidth = 70;
    // while (progress <= 1.0) {

    //     cout << "[";
    //     int pos = barWidth * progress;
    //     for (int i = 0; i < barWidth; ++i) {
    //         if (i < pos) std::cout << "=";
    //         else if (i == pos) std::cout << ">";
    //         else std::cout << " ";
    //     }
    //     cout << "] " << int(progress * 100.0) << " %\r";
    //     cout.flush();

    //     progress += 0.01; // for demonstration only
    // }
    if (progress < 1.0) {

        cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        cout << "] " << int(progress * 100.0) << " %\r";
        cout.flush();

        progress += 0.01; // for demonstration only
    } else {
        cout << "[";
        int pos = barWidth;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) cout << "=";
            else if (i == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << "DONE\r";
        cout.flush();

        cout <<'\n'<< endl;
    }
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
    cout<<'\n';

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
    cout<<"Enter row entries (double values, ensure entries are accurate up to 10 dp) for row "<<row_number<<'\n';

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

pair<double, vector<pair<double, vector<int>>>> performRandomTrialsOfBirkhoffDecompositions(int num_of_trials, matrix user_matrix) {
    double minimal_negativity = INF;
    vector<pair<double, vector<int>>> results;
    vector<pair<double, vector<int>>> chosen_positive_decomposition;
    vector<pair<double, vector<int>>> chosen_ones_decomposition;

    double delta = 0;
    for (row ele : user_matrix) {
        for (double ele2 : ele) {
            delta = min(delta, ele2);
        }
    }
    delta = abs(delta); //make delta positive

    /**
     * S_ = (1+d*delta) A - d*delta B
     */

    matrix positive_part;
    matrix all_one_matrix_normalized;
    
    int num_of_nonzero_entries_in_positive = 0;

    for (row ele : user_matrix) {
        vector<double> temp_row;
        double normalized_ones = 1.0/user_matrix.size(); // 1/d
        double positive_norm_constant = 1.0/(1.0 + user_matrix.size()*delta); // 1/(1+d*delta)
        vector<double> temp_ones_row(user_matrix.size(), normalized_ones);
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

    /**
     * Perform random trials
     */
    double progress = 0.0;
    double step_interval = 1.0/num_of_trials;
    vector<double> all_negativities;

    for (int i = 0; i < num_of_trials; i++) {
        printProgressBar(progress);
        double temp_minimal_negativity = 0.0;
        vector<pair<double, vector<int>>> temp_ones_decomposition;
        vector<pair<double, vector<int>>> temp_positive_decomposition;
        vector<vector<int>> no_preference;
        bool failed = false;

        DFS_random_birkhoff_decomposition(temp_ones_decomposition, all_one_matrix_normalized, user_matrix.size() * user_matrix.size(), no_preference, 0, failed);
        vector<vector<int>> preferred_decomposition;

        for (pair<double, vector<int>> ele : temp_ones_decomposition) {
            preferred_decomposition.push_back(ele.second);
        }

        DFS_random_birkhoff_decomposition(temp_positive_decomposition, positive_part, num_of_nonzero_entries_in_positive, preferred_decomposition, 0, failed);
        if (failed) {
            continue;
        }
        map<string, double> temp_summing_up_coefficients;

        for (pair<double, vector<int>> positive_part_term : temp_positive_decomposition) {
            string map_key = "";
            for (int ele : positive_part_term.second) {
                map_key += to_string(ele)+",";
            }
            temp_summing_up_coefficients.insert({map_key, (1 + user_matrix.size() * delta) * positive_part_term.first});
            // summing_up_coefficients[map_key] = (1 + m.size() * delta) * positive_part_term.first;
        }
        for (pair<double, vector<int>> all_ones_term : temp_ones_decomposition) {
            string map_key = "";
            for (int ele : all_ones_term.second) {
                map_key += to_string(ele)+",";
            }
            if (temp_summing_up_coefficients.find(map_key) == temp_summing_up_coefficients.end()) {
                temp_minimal_negativity += (user_matrix.size() * delta) * all_ones_term.first;
            } else {
                double temp_coefficient = temp_summing_up_coefficients[map_key] - (user_matrix.size() * delta) * all_ones_term.first;
                if (temp_coefficient < 0) {
                    temp_minimal_negativity += abs(temp_coefficient);
                }
            }
        }

        // for debugging to show negativities of all tried decompositions
        // all_negativities.push_back(temp_minimal_negativity); 

        if (temp_minimal_negativity < minimal_negativity) {
            minimal_negativity = temp_minimal_negativity;
            chosen_ones_decomposition = temp_ones_decomposition;
            chosen_positive_decomposition = temp_positive_decomposition;
        }
        progress += step_interval;
    }
    printProgressBar(progress);

    for (double ele : all_negativities) {
        cout<<"tried negativity: "<<ele<<'\n';
    }

    cout<<endl;


    map<string, double> summing_up_coefficients;

    for (pair<double, vector<int>> positive_part_term : chosen_positive_decomposition) {
        string map_key = "";
        for (int ele : positive_part_term.second) {
            map_key += to_string(ele) + ",";
        }
        summing_up_coefficients.insert({map_key, (1 + user_matrix.size() * delta) * positive_part_term.first});
        // summing_up_coefficients[map_key] = (1 + m.size() * delta) * positive_part_term.first;
    }
    for (pair<double, vector<int>> all_ones_term : chosen_ones_decomposition) {
        string map_key = "";
        for (int ele : all_ones_term.second) {
            map_key += to_string(ele) + ",";
        }
        if (summing_up_coefficients.find(map_key) == summing_up_coefficients.end()) {
            results.push_back({(user_matrix.size() * delta * -1) * all_ones_term.first, all_ones_term.second});
        } else {
            summing_up_coefficients[map_key] -= (user_matrix.size() * delta) * all_ones_term.first;

        }
    }

    for (pair<string, double> ele : summing_up_coefficients) {
        results.push_back({ele.second, convertCSVToVector(ele.first)});
    }

    return {minimal_negativity, results};
}


void DFS_random_birkhoff_decomposition(vector<pair<double, vector<int>>> &results, matrix &m, int remaining_entries, 
    vector<vector<int>> &preferred_decomposition, int pref_counter, bool &failed) {
    if (remaining_entries == 0) {
        return;
    }

    int new_zero_entries = 0;
    matrix copy_m = m;
    int counter = 0;

    if (pref_counter < preferred_decomposition.size()) {
        vector<int> preferred_perm = preferred_decomposition[pref_counter];
        double min_ele = m[0][preferred_perm[0]];
        for (int i = 0; i < preferred_perm.size(); i++) {
            if (m[i][preferred_perm[i]] < EPS) {
                if (m[i][preferred_perm[i]] != 0) {
                    counter++;
                }
                copy_m[i][preferred_perm[i]] = 0;
            }
            min_ele = min(min_ele, copy_m[i][preferred_perm[i]]);
        }

        if (min_ele == 0) {
            return DFS_random_birkhoff_decomposition(results, copy_m, remaining_entries - counter, preferred_decomposition, pref_counter + 1, failed);
        }

        for (int i = 0; i < preferred_perm.size(); i++) {
            copy_m[i][preferred_perm[i]] -= min_ele;
            if (copy_m[i][preferred_perm[i]] < EPS) {
                copy_m[i][preferred_perm[i]] = 0;
                new_zero_entries++;
            }
        }
        results.push_back({min_ele, preferred_perm});

        return DFS_random_birkhoff_decomposition(results, copy_m, remaining_entries - new_zero_entries, preferred_decomposition, pref_counter + 1, failed);

    }
    // Use Hopcroft-Karp algorithm to find perfect matching
    BipartiteGraph bipartiteGraph(m.size(), m.size());
    for (int i = 0; i < m.size(); i++) {
        row ele_row = m[i];
        for (int j = 0; j < ele_row.size(); j++) {
            double ele = ele_row[j];
            if (ele > EPS) {
                bipartiteGraph.addEdge(j+1, i+1);
            }
        }
    }
    int num_of_matchings = bipartiteGraph.hopcroftKarp();
    
    if (num_of_matchings < m.size()) {
        cout<<"Hopcroft-Karp algorithm failed to find perfect matching, retrying...\n";
        // results.clear();
        failed = true;
        return;
        // for (int i = 0; i < m.size(); i++) {
        //     for (int j = 0; j < m.size(); j++) {
        //         if (m[i][j] < EPS) {
        //             if (m[i][j] != 0) {
        //                 counter++;
        //             }
        //             copy_m[i][j] = 0;
        //         }
        //     }
        // }
        // return DFS_random_birkhoff_decomposition(results, copy_m, remaining_entries - counter, preferred_decomposition, pref_counter);
    }
    
    vector<int> current_perm = bipartiteGraph.getPermutationMatrix();

    double min_ele = m[0][current_perm[0]];
    for (int i = 0; i < current_perm.size(); i++) {
        if (m[i][current_perm[i]] < EPS) {
            if (m[i][current_perm[i]] != 0) {
                counter++;
            }
            copy_m[i][current_perm[i]] = 0;
        }
        min_ele = min(min_ele, copy_m[i][current_perm[i]]);
    }

    if (min_ele == 0) {
        return DFS_random_birkhoff_decomposition(results, copy_m, remaining_entries - counter, preferred_decomposition, pref_counter, failed);
    }

    for (int i = 0; i < current_perm.size(); i++) {
        copy_m[i][current_perm[i]] -= min_ele;
        if (copy_m[i][current_perm[i]] < EPS) {
            copy_m[i][current_perm[i]] = 0;
            new_zero_entries++;
        }
    }
    results.push_back({min_ele, current_perm});

    return DFS_random_birkhoff_decomposition(results, copy_m, remaining_entries - new_zero_entries, preferred_decomposition, pref_counter, failed);

}

BipartiteGraph::BipartiteGraph(int m, int n) {
    this->m = m;
    this->n = n;
    for (int i = 0; i <= m; i++) {
        this->edges.push_back(vector<int>{});
    }
}

void BipartiteGraph::addEdge(int u, int v) {
    this->edges[u].push_back(v);
}

bool BipartiteGraph::bfs() {
    queue<int> Q;

    for (int i = 1; i <=this->m; i++) {
        if (this->pairM[i] == 0) {
            // unmatched vertex
            this->distance[i] = 0;
            Q.push(i);

        } else {
            this->distance[i] = INT_MAX;
        }
    }

    this->distance[0] = INT_MAX;

    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();

        if (this->distance[u] < this->distance[0]) {
            for (int ele : this->edges[u]) {
                if (this->distance[this->pairN[ele]] == INT_MAX) {
                    // unexplored edge
                    this->distance[this->pairN[ele]] = this->distance[u] + 1;
                    Q.push(this->pairN[ele]);
                }
            }

        }
    }

    return (this->distance[0] != INT_MAX);

}

bool BipartiteGraph::dfs(int u) {
    if (u != 0) {
        vector<int> copy_edges = this->edges[u];
        uid randInt(0, copy_edges.size() -1);
        for (int i = 0; i < copy_edges.size(); i++) {
            int random_int_1 = randInt(rng);
            int random_int_2 = randInt(rng);
            int buffer = copy_edges[random_int_1];
            copy_edges[random_int_1] = copy_edges[random_int_2];
            copy_edges[random_int_2] = buffer;
        }

        for (int ele : copy_edges) {
            if (this->distance[this->pairN[ele]] == this->distance[u] + 1) {
                if (dfs(this->pairN[ele]) == true) {
                    this->pairM[u] = ele;
                    this->pairN[ele] = u;
                    return true;
                }
            }
        }

        this->distance[u] = INT_MAX;
        return false;
    }
    return true;

}

int BipartiteGraph::hopcroftKarp() {
    this->pairM = vector<int>(this->m + 1, 0);
    this->pairN = vector<int>(this->n + 1, 0);
    this->distance = vector<int>(this->m + 1, 0);
    int result = 0;

    while (bfs()) {
        for (int i = 1; i <= this->m; i++) {
            if (this->pairM[i] == 0 && dfs(i)) {
                result++;
            }

        }
    }

    return result;
}

vector<int> BipartiteGraph::getPermutationMatrix() {
    vector<int> permutation_matrix;
    for (int i = 1; i < this->pairN.size(); i++) {
        // Get the corresponding column index for a given row index
        permutation_matrix.push_back(pairN[i] - 1);
    }
    return permutation_matrix;
}


bool less_than_pair_first (const pair<double, pair<int, int>> &p1, const pair<double, pair<int,int>> &p2) {
    return (p1.first < p2.first);
}


pair<double, pair<matrix, matrix>>  performGreedyMinimalNegativityDecomposition(matrix m) {
    double minimal_negativity = 0.0;
    pair<matrix, matrix> positive_negative_decomposition;
    matrix positive_part = m;
    matrix negative_part;

    for (int i = 0; i < m.size(); i++) {
        negative_part.push_back(vector<double>(m.size(), 0));
    }

    vector<pair<double, pair<int, int>>>negative_element_ranking;

    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m.size(); j++) {
            if (m[i][j] < 0) {
                negative_element_ranking.push_back({m[i][j] ,{i,j}});
            }
        }
    }

    sort(negative_element_ranking.begin(), negative_element_ranking.end(), less_than_pair_first);

    // for debugging and viewing negative elements of quasi-bistochastic matrix
    // cout<<"ranking list of negative elements: ";
    // for (pair<double, pair<int, int>> ele : negative_element_ranking) {
    //     cout<<ele.first<<" ";
    // }

    // cout<<endl;

    while (!negative_element_ranking.empty()) {
        vector<pair<double, pair<int, int>>> negative_element_ranking_temp = negative_element_ranking;

        vector<int> current_perm(m.size(), -1);
        set<int> perm_set_row; // row index
        set<int> perm_set_column; // column index

        pair<double, pair<int, int>> most_negative_element = negative_element_ranking.front();

        perm_set_column.insert(most_negative_element.second.first);
        perm_set_row.insert(most_negative_element.second.second);
        current_perm[most_negative_element.second.first] = most_negative_element.second.second;
        // remove most_negative_element from ranking
        negative_element_ranking.erase(negative_element_ranking.begin());

        for (int i = 0; i < negative_element_ranking_temp.size(); i++) {
            pair<double, pair<int, int>> ele = negative_element_ranking_temp[i];
            if (perm_set_row.find(ele.second.second) != perm_set_row.end() || 
                perm_set_column.find(ele.second.first) != perm_set_column.end()) {
                continue;
            }
            current_perm[ele.second.first] = ele.second.second;
            perm_set_column.insert(ele.second.first);
            perm_set_row.insert(ele.second.second);
            // remove ele from ranking
            vector<pair<double, pair<int, int>>>::iterator pos = find(negative_element_ranking.begin(), negative_element_ranking.end(), ele);
            negative_element_ranking.erase(pos);
        }

        if (perm_set_row.size() < m.size()) {
            for (int i = 0; i < m.size(); i++) {
                if (current_perm[i] == -1) {
                    for (int j = 0; j < m.size(); j++) {
                        if (perm_set_row.find(j) == perm_set_row.end()) {
                            current_perm[i] = j;
                            perm_set_row.insert(j);
                            perm_set_column.insert(i);
                            break;
                        }
                    }
                }
            }
        }

        // successfully created perm matrix
        for (int i = 0; i < m.size(); i++) {
            // update positive part and negative part
            positive_part[i][current_perm[i]] += abs(most_negative_element.first);
            negative_part[i][current_perm[i]] -= abs(most_negative_element.first);
        }

        minimal_negativity += abs(most_negative_element.first);

    }
    positive_negative_decomposition.first = positive_part;
    positive_negative_decomposition.second = negative_part;

    return {minimal_negativity, positive_negative_decomposition};
}

void printPositiveNegativeResults(double negativity, pair<matrix, matrix> positive_negative_decomposition) {
    cout<<"\n======POSITIVE MATRIX PART======\n";
    for (row ele : positive_negative_decomposition.first) {
        for (double ele2 : ele) {
            cout<<fixed<<setprecision(3)<<ele2<<"\t";
        }
        cout<<"\n\n";
    }

    cout<<"\n======NEGATIVE MATRIX PART======\n";
    for (row ele : positive_negative_decomposition.second) {
        for (double ele2 : ele) {
            cout<<fixed<<setprecision(3)<<ele2<<"\t";
        }
        cout<<"\n\n";
    }

    cout<<"The negativity of the greedy minimal negativity decomposition is \n\n";
    cout<<fixed<<setprecision(5)<<negativity<<"\n\n";

}


bool verifyDecomposition(matrix m, matrix positive, matrix negative) {
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m.size(); j++) {
            if (abs(m[i][j] - positive[i][j] - negative[i][j]) > EPS) {
                cout<<"leftover difference: "<<abs(m[i][j] - positive[i][j] - negative[i][j])<<endl;
                return false;
            }
        }
    }
    return true;
}



pair<double, pair<matrix, matrix>> calculateMostNegativeSingleQubitUnitary(double step_interval) {
    double max_negativity = 0.0;
    matrix chosen_rotation;
    matrix chosen_quasi_bistochastic;

    double chosen_alpha = 0.0;
    double chosen_beta = 0.0;
    double chosen_gamma = 0.0;

    matrix rotation;
    matrix quasi_bistochastic;

    for (int i = 0; i < 3; i++) {
        rotation.push_back(vector<double>(3, 0));
    }

    for (int i = 0; i < 4; i++) {
        quasi_bistochastic.push_back(vector<double>(4, 0));
    }
    vector<double> n_0 = {0, sqrt(2.0/3), 1/sqrt(3)};
    vector<double> n_1 = {0, -1*sqrt(2.0/3), 1/sqrt(3)};
    vector<double> n_2 = {sqrt(2.0/3), 0, -1/sqrt(3)};
    vector<double> n_3 = {-1*sqrt(2.0/3), 0, -1/sqrt(3)};

    vector<vector<double>> n_set;
    n_set.push_back(n_0);
    n_set.push_back(n_1);
    n_set.push_back(n_2);
    n_set.push_back(n_3);

    double progress = 0.0;
    double progress_interval = step_interval/2/PI;

    // construct rotation matrix using Tait-Bryan angles
    for (double alpha = 0; alpha < 2*PI; alpha += step_interval) {
        printProgressBar(progress);
        for (double beta = 0; beta < 2*PI; beta += step_interval) {
            for (double gamma = 0; gamma < 2*PI; gamma += step_interval) {
                rotation[0][0] = cos(alpha) * cos(beta);
                rotation[0][1] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma);
                rotation[0][2] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma);
                rotation[1][0] = sin(alpha) * cos(beta);
                rotation[1][1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma);
                rotation[1][2] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma);
                rotation[2][0] = -1 * sin(beta);
                rotation[2][1] = cos(beta) * sin(gamma);
                rotation[2][2] = cos(beta) * cos(gamma);

                vector<double> n_a;
                vector<double> n_b;
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        quasi_bistochastic[i][j] =  1.0/4 + 3.0/4 * dot_product(matrix_vector_multiply(rotation, n_set[i]), n_set[j]);
                    }
                }

                pair<double, pair<matrix, matrix>>temp_greedy_decomposition = performGreedyMinimalNegativityDecomposition(quasi_bistochastic);
                if (temp_greedy_decomposition.first > max_negativity) {
                    max_negativity = temp_greedy_decomposition.first;
                    chosen_quasi_bistochastic = quasi_bistochastic;
                    chosen_rotation = rotation;

                    chosen_alpha = alpha;
                    chosen_beta = beta;
                    chosen_gamma = gamma;
                }
            }
        }

        progress += progress_interval;
    }
    printProgressBar(progress);

    cout<<"\n Tait-Bryan Angles (alpha, beta, gamma): " <<chosen_alpha<<" "<<chosen_beta<<" "<<chosen_gamma<<"\n\n";

    return {max_negativity, {chosen_rotation, chosen_quasi_bistochastic}};
}

void printMostNegativeSingleQubitRotation(double negativity, matrix rotation, matrix quasi_bistochastic) {
    cout<<"\n======ROTATION MATRIX======\n";
    for (row ele : rotation) {
        for (double ele2 : ele) {
            cout<<fixed<<setprecision(4)<<ele2<<"\t";
        }
        cout<<"\n\n";
    }

    cout<<"\n======QUASI-BISTOCHASTIC MATRIX======\n";
    for (row ele : quasi_bistochastic) {
        for (double ele2 : ele) {
            cout<<fixed<<setprecision(4)<<ele2<<"\t";
        }
        cout<<"\n\n";
    }

    cout<<"The negativity of the most negative single qubit unitary using \ngreedy minimal negativity decomposition is \n\n";
    cout<<fixed<<setprecision(5)<<negativity<<"\n\n";
}



double dot_product(vector<double> v1, vector<double> v2) {
    double result = 0.0;
    for (int i = 0; i < v1.size(); i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

vector<double> matrix_vector_multiply(matrix m, vector<double> v) {
    vector<double> result(m.size(), 0.0);
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m.size(); j++) {
            result[i] += m[i][j] * v[j];
        }
    }

    return result;
}