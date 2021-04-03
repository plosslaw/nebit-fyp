#include "QuasiBistochasticBirkhoffDecomposition.h"
#include <sstream>
#include <string>

using namespace std;


const double EPS = 1e-7;

int main() {
    printProgressBar();
    matrix user_matrix;

    int option = parseProgramMode();
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

        break;
    case 2:
        cout<<"\nOption 2 chosen\n\n";
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
            cout<<"Please enter an integer greater than 0\n";
        }
    }

    for (int i = 0; i < size; i++) {
        vector<double> v = parseMatrixRows(i, size);
        m.push_back(v);
    }
    cout<<"\n==========Input Matrix==========\n\n";
    for (row ele : m) {
        for (double ele2 : ele) {
            cout<<ele2<<"\t";
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

void printNegativityResults(double negativity, vector<pair<double, string>> min_negativity_decomposition) {
    cout<<fixed<<setprecision(5)<<negativity<<endl;
}







vector<pair<double, string>> get_minimal_negativity_birkhoff_decomposition(matrix m) {
    vector<pair<double, string>> result;
    double delta = 0;
    for (row ele : m) {
        for (double ele2 : ele) {
            delta = min(delta, ele2);
        }
    }

    /**
     * S_ = (1+d*delta) A - d*delta B
     */

    matrix positive_part;
    matrix all_one_matrix_normalized;

    for (row ele : m) {
        vector<double> temp_row;
        double normalized_ones = 1.0/m.size(); // 1/d
        double positive_norm_constant = 1.0/(1.0 + m.size()*delta); // 1/(1+d*delta)
        vector<double> temp_ones_row(m.size(), normalized_ones);
        for (double ele2 : ele) {
            if ((ele2 + delta)*positive_norm_constant < EPS) {
                temp_row.push_back(0);
            } else {
                temp_row.push_back((ele2 + delta)*positive_norm_constant);
            }
        }
        positive_part.push_back(temp_row);
        all_one_matrix_normalized.push_back(temp_ones_row);
    }

    vector<vector<pair<double, string>>> positive_part_normalized_decomposition =
        DFS_all_birkhoff_decomposition(positive_part);
    
    vector<vector<pair<double, string>>> all_ones_normalized_decomposition = 
        DFS_all_birkhoff_decomposition(all_one_matrix_normalized);


    for (vector<pair<double, string>> positive_part_ele : positive_part_normalized_decomposition) {
        for (vector<pair<double, string>> all_ones_ele : all_ones_normalized_decomposition) {

        }
    }

    return result;
}

vector<vector<pair<double, string>>> DFS_all_birkhoff_decomposition(matrix &m) {

}
