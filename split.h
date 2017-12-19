std::vector<std::string> split_into_words(std::string line) {
	std::stringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);
	return vstrings;
}

std::vector<double> split_into_doubles(std::string line) {
	std::stringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);
	std::vector<double> vnumbers;
	
	for(std::vector<std::string>::iterator it = vstrings.begin(); it != vstrings.end(); ++it)
	{
		std::stringstream geek(*it); double num; geek >> num; vnumbers.push_back(num);
	}
	return vnumbers;
}

std::vector<int> split_into_integers(std::string line) {
	std::stringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);
	std::vector<int> vnumbers;
	
	for(std::vector<std::string>::iterator it = vstrings.begin(); it != vstrings.end(); ++it)
	{
		std::stringstream geek(*it); int num; geek >> num; vnumbers.push_back(num);
	}
	return vnumbers;
}
