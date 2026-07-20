// Copyright (C) 2017 by SIB Swiss Institute of Bioinformatics, Julien Dorier and Dimos Goundaroulis.
// 
// This file is part of project Knoto-ID.
// 
// Knoto-ID is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// Knoto-ID is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Knoto-ID.  If not, see <http://www.gnu.org/licenses/>.

#ifndef KNOTTEDCORE_HH
#define KNOTTEDCORE_HH

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <map>
#include <ctime>
#include <Polygon.hh>
#include <PlanarDiagram.hh>
#include <PolynomialInvariant.hh>
// #include <Version.h>
#include <Random.hh>
// #include <KnotoidNames.hh>
#include <boost/regex.hpp>
#include <vector>
#include <string>
// #include <fstream>
#include <sstream>
// #include <stdexcept>

class KnottedCore {
public:
    KnottedCore(
        const std::vector<double>& x_coords,
        const std::vector<double>& y_coords,
        const std::vector<double>& z_coords,
        const std::string& formattedData,
        const std::vector<std::vector<double> > projectionlist_projections,
        const std::vector<double> projectionlist_weights,
        const bool& all_chains = false,
        const std::string& closure_method = "open",
        const bool& planar = true,
        const bool& arrow_polynomial = true,
        const bool& cyclic = false,
        const bool& cyclic_input = false,
        const bool& kc_search_path = false,
        const bool& reduction_3d=true,
        const bool& simplify_diagram=true,
        const time_t timeout = 1,
        const std::string& jones_method = "recursive",
        const long& max_nb_random_moves_III = 100000,
        const long& max_nb_unsuccessfull_random_moves_III = 2000,
        const bool& debug = false
    );

std::tuple< std::map<std::string, std::vector<std::string>>, std::map<std::string, std::vector<std::string>>, std::map<std::string, std::vector<std::string>>> Run();

private:
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    std::vector<double> z_coords;
    std::string formattedData;
    std::vector<std::vector<double> > projectionlist_projections;
    std::vector<double> projectionlist_weights;
    bool all_chains;
    std::string closure_method;
    bool planar;
    bool arrow_polynomial;
    bool cyclic;
    bool cyclic_input;
    const bool kc_search_path;
    const bool reduction_3d;
    const bool simplify_diagram;
    const time_t timeout;
    std::string jones_method;
    long max_nb_random_moves_III;
    long max_nb_unsuccessfull_random_moves_III;
    bool debug;

    // std::string readFileContents(const std::string& filename);
    std::string get_jones(Polygon& polygon, double& frequ, unsigned long nb_projections, std::string closure_method);
    std::vector<std::string> split_string(std::string str,std::string sep_list);
    std::vector<std::string> split_input(std::string str);
    // int mod(int x,int N){if(x>=0)return x%N;else return (x+N*(2-x/N))%N;}
    // Other private methods and members as needed for processing
};

#endif // KNOTTEDCORE_H
