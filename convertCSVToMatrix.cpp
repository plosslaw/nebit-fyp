#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

int main() {
    int size_matrix;
    cout<<"Input size of matrix\n\n";
    cin>>size_matrix;
    
    string input_csv;
    cout<<"\nInput CSV data\n\n";
    cin>>input_csv;

    freopen("output.txt", "w", stdout);
    vector<vector<string>> results;
    for (int i = 0; i < size_matrix; i++) {
        results.push_back(vector<string>{});

    }
    string temp_string;
    string csv_copy = input_csv;
    if (input_csv.at(input_csv.size()-1) == ',') {
        csv_copy.pop_back();
    }
    stringstream ss(csv_copy);
    int row_number = 0;
    int counter = 0;
    while (getline(ss >> std::ws, temp_string, ',')) {
        results[row_number].push_back(temp_string);
        counter++;
        if (counter >= size_matrix) {
            row_number++;
            counter = 0;
        }
    }

    for (vector<string> ele_vector : results) {
        for (string ele : ele_vector) {
            cout<<ele<<" ";
        }
        cout<<'\n';
    }


    fclose(stdout);

    return 0;
}
