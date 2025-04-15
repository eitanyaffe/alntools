#include <iostream>
#include <fstream>
using namespace std;

void save_number(const string& filename, size_t number) {
    ofstream file(filename, ios::binary);
    if (!file.is_open()) {
        cerr << "Error opening file for writing: " << filename << endl;
        abort();
    }

    cout << "Saving number: " << number << endl;
    file.write(reinterpret_cast<const char*>(&number), sizeof(number));
    file.close();
}

size_t load_number(const string& filename) {
    ifstream file(filename, ios::binary);
    if (!file.is_open()) {
        cerr << "Error opening file for reading: " << filename << endl;
        abort();
    }

    size_t number;
    file.read(reinterpret_cast<char*>(&number), sizeof(number));
    file.close();

    cout << "Loaded number: " << number << endl;
    return number;
}

int xmain() {
    string filename = "test_number.bin";
    size_t original_number = 42;

    save_number(filename, original_number);
    size_t loaded_number = load_number(filename);

    if (original_number == loaded_number) {
        cout << "Success! Numbers match." << endl;
    } else {
        cout << "Error: Numbers do not match." << endl;
    }

    return 0;
}
