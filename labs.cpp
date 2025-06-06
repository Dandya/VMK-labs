#include "labs.h"

std::vector<int>
Labs::CreateVecRand(int two_pow_n){
	std::vector<int> s_vec;
	srand(time(nullptr));
	for(int i=0; i<two_pow_n; i++){
		int buf = rand()%two_pow_n;
		while (std::find(s_vec.begin(), s_vec.end(), buf) != s_vec.end())
			buf = rand()%two_pow_n;
		s_vec.push_back(buf);
	}
	return s_vec;
}

std::vector<int>
Labs::GetFuncF(int numb_bit, const std::vector<int>& s_vec){
	std::vector<int> f;
	int bit = 1 << numb_bit;
	for(int i = 0; i < s_vec.size(); i++) {
		int buf = bit & s_vec[i];
		buf >>= numb_bit;
		f.push_back(buf);
	}
	return f;
}

int
Labs::GetWeight(const std::vector<int>& vec) {
	int counter=0;
	for (int i = 0; i < vec.size(); i++)
		if(vec[i])
			counter++;
	return counter;
}