#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <ctime>

#define N_NUMBER 6

struct Labs {
	std::vector<int> s_vec;
	std::vector<std::vector<int>> f_vec;

	Labs() {}

	Labs(const std::vector<int> s) : s_vec(s) {}

	void Lab1();

	void Lab2();

	void Lab3();

	void Lab4();

 private:
	std::vector<int> CreateVecRand(int two_pow_n);

	std::vector<int> GetFuncF(int numb_bit, const std::vector<int>& s_vec);

	int GetWeight(const std::vector<int>& f_vec);

	// lab2.cpp
	bool IsHighlyEquiprobable(const std::vector<int>& vec);

	// lab3.cpp
	std::vector<double> CreateFurieKoef(const std::vector<int>& f_vec);

	std::pair<bool, int> IsCorrelativeImmunity(const std::vector<double>& walsh);

	std::pair<bool, int> Elasticity (const std::vector<int>& f_vec,
			std::pair<bool, int> corr_immun, const std::vector<double>& walsh);

	void PrintSpectrum(const std::vector<std::vector<int>>& delta);

	std::vector<std::string> BestLinearApproximations(const std::vector<int>& stat_coeff);

	template<typename T>
	void PrintVec(const std::vector<T>& vec, const std::string& desc) {
		std::cout << desc <<": ";
		for(int i = 0; i < vec.size(); i++)
				std::cout << vec[i] << " ";
		std::cout << "\n";
	}

	std::string GetBits(unsigned int number, unsigned int len) {
		std::string bits = "";
		for (size_t i = len; i > 0; i--) {
			bits += ((number & (1 << i)) > 0) ? '1' : '0';
		}
		return bits;
	}

	std::vector<int> GetBitsVec(int number, int bit_len){
		std::vector<int> bits(bit_len);

		for(int i = 0; i < bit_len; ++i){
			bool bit = (number & (1 << i)) > 0;
			bits[bit_len-i-1] = bit ? 1 : 0;
		}
		return bits;
	}

	double GetPredominanceZero(const std::vector<int>& f_vec) {
		// d(f) = 1 - 2P{f=1} = 1 - ||f||/2^{N_NUMBER-1}
		return 1.0 - static_cast<double>(GetWeight(f_vec))/pow(2, N_NUMBER-1);
	}
};