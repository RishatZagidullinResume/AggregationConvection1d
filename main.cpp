#include <algorithm>
#include <memory>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <utility>
#include "solver/solver.h"
#include <cstdlib>
#include <limits>
#include <iterator>
#include <sstream>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

template<typename T>
std::vector<T> split(const std::string& line)
{
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is), 
                          std::istream_iterator<T>());
}

int main(int argc, char ** argv)
{
    const auto stream_max = 
                   std::numeric_limits<std::streamsize>::max();
    std::ofstream output;
    std::ifstream input;
    if (argc !=4)
    {
        std::cout << "need 2 output file names" << std::endl;
        return 2;
    }
    solver::AdvDifCoag1d solver = solver::AdvDifCoag1d();

    std::ofstream check;
    check.open(argv[3]);

    int append = atoi(argv[2]);
    if (append != 0){
        input.open(argv[1]);
        output.open(argv[1], std::ios_base::app);
        for (int i = 0; i < append*solver.max_x; i++)
            input.ignore(stream_max, '\n');
        std::string line;
        int cc = 0;
        while (std::getline(input, line)){
            std::vector<double> vec = split<double>(line);
            for (int m = 0; m < solver.size; m++)
                solver.init_layer[m+solver.size*cc] = vec[m];
            cc++;
        }
        input.close();
    } else
        output.open(argv[1]);

    for (int t = 0; t < solver.TIME; t++)
    {
        printProgress((double)t/solver.TIME);
        solver.iteration(output, check, t);
    }

    output.close();
    check.close();
    return 0;
}

