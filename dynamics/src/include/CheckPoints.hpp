/*!
 * @file CheckPoints.hpp
 * @date 1 Jun 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __CHECKPOINTS_HPP
#define __CHECKPOINTS_HPP

#include "dgVector.hpp"
#include "cgVector.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>



namespace Nextsim {

/*!
 * This namespace collects the auxiliary routines
 */
namespace CheckPoints {


template<int DG, template<int> class VectorType>
void saveData(const std::string& fileName, int n, VectorType<DG>& matrix)
{

	const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
    
    std::ostringstream ss;
    ss << fileName << "." << std::setw(5) << std::setfill('0') << n << ".txt";

	std::ofstream file(ss.str());
	assert(file.is_open());
    
		file << matrix.format(CSVFormat);
		file.close();
	
}

template<int DG>
void loadData (const std::string & path, Nextsim::CellVector<DG>& matrix) {
    
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }

#pragma omp parallel for
    for (size_t i = 0; i < rows; ++i) 
      for (size_t j = 0; j < DG; ++j) 
          matrix(i,j) = values[i*DG+j];

}


template<int CG>
void loadData (const std::string & path, Nextsim::CGVector<CG>& matrix) {
    
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }

#pragma omp parallel for
    for (size_t i = 0; i < rows; ++i) 
          matrix(i,0) = values[i];
}


} /* namespace CheckPoints */

} /* namespace Nextsim */

#endif /* __CHECKPOINTS_HPP */
