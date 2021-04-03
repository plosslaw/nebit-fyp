#include "QuasiStochasticBirkhoffDecomposition.h"

using namespace std;

int main() {
    string welcome_message = "\nWhat would you like to do?\n";
    string first_option_message = 
        "\t1. Compute the negativity of quasi-bistochastic matrix using \n\tminimal negativity Birkhoff-von Neumann Decomposition\n";
    string second_option_message = 
        "\t2. Compute the most negative single qubit unitary operation \n\tusing minimal negativity Birkhoff-von Neumann Decomposition\n";
    string third_option_message = "\t3. Exit program\n\n";
    string start_message = welcome_message + first_option_message + second_option_message + third_option_message;
    cout<<start_message;

    
    while (true) {
        string input_string;
        getline(cin,input_string);
        char *p;
        long input = strtol(input_string.c_str(), &p, 10);
        while (*p) {
            cout<<"Invalid input, please enter the index number of the options listed\n";
            getline(cin,input_string);
            input = strtol(input_string.c_str(), &p, 10);
        }
        if (input == 1) {
            cout<<"option 1 chosen\n";
            break;
        } else if (input == 2) {
            cout<<"option 2 chosen\n";
            break;
        } else if (input == 3) {
            cout<<"\nProgram will now end. Have a nice day.\n\n";
            break;
        } else {
            cout<<"Invalid option, please choose one of the options listed\n";
        }
        
    }


    return 0;
}