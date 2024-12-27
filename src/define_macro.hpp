#ifndef DEFINE_MACRO_HPP_
#define DEFINE_MACRO_HPP_

#include <iostream>
#include <fstream>

#define FAILED_TO_OPEN(file_stream, file_name) if(!file_stream) {std::cerr << "##ERROR : Failed to open " << file_name << std::endl; exit(1); }

#define SQR(x) ( (x)*(x) )
#define CUB(x) ( (x)*(x)*(x) )
#define SIGN(x) ( (x >= 0.0 ? 1.0 : -1.0))

#endif /* DEFINE_MACRO_HPP_ */