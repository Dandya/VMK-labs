
#include "labs.h"

// extern std::vector<std::vector<int>> Lab1();
// extern void Lab2(std::vector<std::vector<int>> f_vec);
// extern std::pair<int, std::vector<int>> Lab3(std::vector<std::vector<int>> f_vec);
// extern void Lab4();

int main(){
	// Test data
	std::vector<int> s = {
			11, 16, 12, 21, 47, 31, 36, 39, 26, 0,  40, 25, 1,  22, 55, 6,
			57, 4,  62, 27, 15, 7,  61, 51, 34, 63, 50, 38, 20, 33, 49, 5,
			13, 41, 53, 2,  54, 43, 10, 8,  42, 35, 60, 48, 9,  14, 30, 44,
			23, 52, 3,  29, 58, 19, 28, 46, 18, 32, 56, 37, 59, 24, 17, 45};
	Labs labs(s);
	// Labs labs;

	labs.Lab1();
	labs.Lab2();
	labs.Lab3();

        // std::vector<std::vector<int>> f_vec = Lab1();

	// cout << "\n\n\tЛабораторная работа №2\n\n" << endl;
	// Lab2(f_vec);

	// cout << "\n\n\tЛабораторная работа №3\n\n" << endl;
	// pair<int, std::vector<int>> non_corr_immun_func = Lab3(f_vec);

	// cout << "\n\n\tЛабораторная работа №4\n\n" << endl;
	/*f_vec[non_corr_immun_func.first]={0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0,
																			0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0,
																			0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0,
																			0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0};*/
	// Lab4(/*f_vec[non_corr_immun_func.first], non_corr_immun_func.second*/);
}