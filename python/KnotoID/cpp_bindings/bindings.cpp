#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <pybind11/iostream.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "Random.hh" 
#include "Polygon.hh"
#include "PolynomialInvariant.hh"
#include "PlanarDiagram.hh"
#include "KnottedCore.hh"
#include "DoubleBranchedCover.hh"

namespace py = pybind11;

// Route C++ std::cout/std::cerr into Python's sys.stdout/sys.stderr for the
// duration of a call, so the Python side can capture or silence them (see
// KnotoID._suppress_cpp_output). Applied as a call_guard on the methods below.
using ostream_redirect = py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>;

// Setup random number generator
boost::random::mt19937 rng; // Mersenne Twister RNG

void set_seed(unsigned int seed) {
    rng.seed(seed);
};


// Pybind11 module
PYBIND11_MODULE(knotoID_cpp, m) {

    m.def("set_seed", &set_seed, "Set the seed for the random number generator");

    py::class_<PolynomialInvariant>(m, "PolynomialInvariant")
        .def(py::init<PlanarDiagram&, bool, bool, bool>(),
             py::arg("diagram"), py::arg("flag_planar"),
             py::arg("flag_arrow_polynomial") = false, py::arg("flag_debug") = false,
             ostream_redirect())
        .def("set_timeout", &PolynomialInvariant::set_timeout)
        .def("get_timeout", &PolynomialInvariant::get_timeout)
        .def("get_polynomial_simple", &PolynomialInvariant::get_polynomial_simple, ostream_redirect())
        .def("get_polynomial_recursive", &PolynomialInvariant::get_polynomial_recursive,
             py::arg("method") = "default", py::arg("flag_silent") = false,
             ostream_redirect())
        ;


    py::class_<Polygon>(m, "Polygon")
        .def(py::init([](py::array_t<double> points, bool flag_cyclic, bool flag_debug) {
            py::buffer_info info = points.request();
            if (info.ndim != 2 || info.shape[1] != 3) {
                throw std::runtime_error("Input should be an nx3 array");
            }

            auto ptr = static_cast<double *>(info.ptr);
            std::vector<double> x, y, z;
            for (size_t i = 0; i < info.shape[0]; i++) {
                x.push_back(ptr[i * 3]);
                y.push_back(ptr[i * 3 + 1]);
                z.push_back(ptr[i * 3 + 2]);
            }
            return Polygon(x, y, z, flag_cyclic, flag_debug);
        }), py::arg("points"), py::arg("flag_cyclic"), py::arg("flag_debug") = false, ostream_redirect())
        .def("set_debug", &Polygon::set_debug)
        .def("set_closure", &Polygon::set_closure, py::arg("dx"), py::arg("dy"), py::arg("dz"), py::arg("method") = "direct", ostream_redirect())
        .def("get_nb_points", &Polygon::get_nb_points)
        .def("simplify_polygon", &Polygon::simplify_polygon, py::arg("dx") = 0, py::arg("dy") = 0, py::arg("dz") = 1, ostream_redirect())
        .def("get_planar_diagram", &Polygon::get_planar_diagram, py::arg("dx"), py::arg("dy"), py::arg("dz"), py::arg("flag_planar"), ostream_redirect())
        .def("project", &Polygon::project, py::arg("dx") = 0, py::arg("dy") = 0, py::arg("dz") = 1, ostream_redirect())
        // .def("get_polygon", static_cast<Polygon (Polygon::*)(int, int)>(&Polygon::get_polygon), py::arg("start_index"), py::arg("end_index"))
        .def("get_polygon", static_cast<Polygon (Polygon::*)(int, int, bool)>(&Polygon::get_polygon), py::arg("start_index"), py::arg("end_index"), py::arg("flag_cyclic"), ostream_redirect())
        ;

    py::class_<PlanarDiagram>(m, "PlanarDiagram")
    .def(py::init<bool, bool>(), py::arg("flag_planar"), py::arg("flag_debug") = false, ostream_redirect())
    .def(py::init<std::vector<Crossing>&, bool, bool>(), py::arg("crossings"), py::arg("flag_planar"), py::arg("flag_debug") = false, ostream_redirect())
    .def("set_debug", &PlanarDiagram::set_debug)
    .def("get_nb_crossings", &PlanarDiagram::get_nb_crossings, py::arg("flag_ignore_endpoints") = false)
    .def("get_nb_arcs", &PlanarDiagram::get_nb_arcs)
    .def("get_crossings", &PlanarDiagram::get_crossings)
    .def("load_from_string_extended_gauss_code", &PlanarDiagram::load_from_string_extended_gauss_code, py::arg("str"), py::arg("flag_cyclic"), ostream_redirect())
    .def("load_from_pd_code",
         [](PlanarDiagram &self,
            const std::vector<std::vector<int>> &crossings,
            const std::vector<int> &exterior_arcs) {
             // Build a PD[X[...]] / KnotTheory string from a nested list of
             // crossings (each = 4 arc labels) and reuse the C++ parser.
             // Arcs listed in exterior_arcs are wrapped in r[] (needed only for
             // planar open knotoids). PD codes carry no signs: the over/under
             // strand is encoded by the order of the four arc labels.
             std::set<int> ext(exterior_arcs.begin(), exterior_arcs.end());
             std::ostringstream oss;
             oss << "PD[";
             for (size_t c = 0; c < crossings.size(); ++c) {
                 if (crossings[c].size() != 4)
                     throw std::invalid_argument("each crossing must be a list of exactly 4 arc labels");
                 if (c) oss << ",";
                 oss << "X[";
                 for (int k = 0; k < 4; ++k) {
                     if (k) oss << ",";
                     int arc = crossings[c][k];
                     if (arc < 0)
                         throw std::invalid_argument("arc labels must be non-negative integers");
                     if (ext.count(arc)) oss << "r[" << arc << "]";
                     else oss << arc;
                 }
                 oss << "]";
             }
             oss << "]";
             self.load_from_string_KnotTheory(oss.str());
         },
         py::arg("crossings"), py::arg("exterior_arcs") = std::vector<int>(),
         ostream_redirect())
    .def("to_pd_code",
         [](PlanarDiagram &self) {
             // Export regular crossings as [a0,a1,a2,a3] rows (0-based, same
             // KnotTheory column order Knoodle uses). Endpoints/empty crossings
             // are skipped. No sign column: the order encodes the over/under.
             std::vector<std::vector<int>> out;
             for (auto &c : self.get_crossings()) {
                 int a0 = c.get_arc(0), a1 = c.get_arc(1), a2 = c.get_arc(2), a3 = c.get_arc(3);
                 if (a0 >= 0 && a1 >= 0 && a2 >= 0 && a3 >= 0)
                     out.push_back({a0, a1, a2, a3});
             }
             return out;
         },
         ostream_redirect())
    .def("get_planar_graph", &PlanarDiagram::get_planar_graph, py::arg("nb_additional_nodes_per_arc_noloop") = 0, py::arg("nb_additional_nodes_per_arc_empty_loop") = 1, py::arg("nb_additional_nodes_per_arc_enclosing_loop") = 3, ostream_redirect())
    .def("is_cyclic", &PlanarDiagram::is_cyclic)
    .def("change_orientation", &PlanarDiagram::change_orientation, ostream_redirect())
    .def("reorder_crossings", &PlanarDiagram::reorder_crossings, ostream_redirect())
    .def("close", &PlanarDiagram::close, py::arg("flag_overpass") = false, ostream_redirect())
    .def("simplify", &PlanarDiagram::simplify, ostream_redirect())
    .def("simplify_with_random_reidemeister_moves_III", &PlanarDiagram::simplify_with_random_reidemeister_moves_III, py::arg("maxiterations"), py::arg("maxiterations_unsuccessfull"), ostream_redirect())
    .def("create_regions", &PlanarDiagram::create_regions, ostream_redirect())
    .def("get_writhe", &PlanarDiagram::get_writhe, ostream_redirect())
    .def("check", &PlanarDiagram::check, ostream_redirect())
    ;

    py::class_<Polynomial>(m, "Polynomial")
    .def(py::init<>())
    .def(py::init<const std::vector<std::string>&>(), py::arg("var_names"))
    .def(py::init<const std::string&, const std::vector<std::string>&>(),
            py::arg("string_polynomial"), py::arg("var_names") = std::vector<std::string>(),
            ostream_redirect())
    .def("add", (void (Polynomial::*)(const Polynomial&, bool)) &Polynomial::add, 
            py::arg("p"), py::arg("flag_force") = false)
    .def("add", (void (Polynomial::*)(double, std::string, long)) &Polynomial::add, 
            py::arg("coefficient"), py::arg("var") = "", py::arg("exponent") = 1)
    .def("multiply", &Polynomial::multiply, 
            py::arg("p"), py::arg("flag_force") = false)
    .def("get_variable_names", &Polynomial::get_variable_names)
    .def("load_from_string", &Polynomial::load_from_string,
            py::arg("input"), py::arg("flag_keep_var_names") = false,
            ostream_redirect())
    .def("to_string", &Polynomial::to_string)
    .def("__eq__", &Polynomial::operator==, py::is_operator())
    .def("__repr__", [](const Polynomial& p) {
        return "<Polynomial: " + p.to_string() + ">";
    })
    ;

py::class_<KnottedCore>(m, "KnottedCore")
    .def(py::init([](py::array_t<double> input, 
                     const std::string& formattedData, 
                     const std::vector<std::vector<double>>& projectionlist_projections, 
                     const std::vector<double>& projectionlist_weights, 
                     bool all_chains, 
                     const std::string& closure_method, 
                     bool planar, 
                     bool arrow_polynomial, 
                     bool cyclic, 
                     bool cyclic_input,
                     bool kc_search_path, 
                     bool reduction_3d, 
                     bool simplify_diagram, 
                     time_t timeout, 
                     const std::string& jones_method, 
                     long max_nb_random_moves_III, 
                     long max_nb_unsuccessfull_random_moves_III, 
                     bool debug) {
        py::buffer_info info = input.request();
        if (info.ndim != 2 || info.shape[1] != 3) {
            throw std::runtime_error("Input should be a 2D array with 3 columns");
        }

        std::vector<double> x_coords(info.shape[0]), y_coords(info.shape[0]), z_coords(info.shape[0]);
        double *ptr = static_cast<double *>(info.ptr);
        for (size_t i = 0; i < info.shape[0]; ++i) {
            x_coords[i] = ptr[i * 3];
            y_coords[i] = ptr[i * 3 + 1];
            z_coords[i] = ptr[i * 3 + 2];
        }

        return new KnottedCore(x_coords, y_coords, z_coords, formattedData, projectionlist_projections, projectionlist_weights, all_chains, closure_method, planar, arrow_polynomial, cyclic, cyclic_input, kc_search_path, reduction_3d, simplify_diagram, timeout, jones_method, max_nb_random_moves_III, max_nb_unsuccessfull_random_moves_III, debug);
    }), 
    py::arg("input"),
    py::arg("formattedData"),
    py::arg("projectionlist_projections"),
    py::arg("projectionlist_weights"),
    py::arg("all_chains") = false,
    py::arg("closure_method") = "open",
    py::arg("planar") = true,
    py::arg("arrow_polynomial") = true,
    py::arg("cyclic") = false,
    py::arg("cyclic_input") = false,
    py::arg("kc_search_path") = false,
    py::arg("reduction_3d") = true,
    py::arg("simplify_diagram") = true,
    py::arg("timeout") = 1,
    py::arg("jones_method") = "recursive",
    py::arg("max_nb_random_moves_III") = 100000,
    py::arg("max_nb_unsuccessfull_random_moves_III") = 2000,
    py::arg("debug") = false)
    .def("run", &KnottedCore::Run, ostream_redirect());
    // py::bind_map<std::map<std::string, std::vector<std::string>>>(m, "StringVectorMap");


    // ----- double branched cover invariant for planar knotoids -----
    py::class_<DBCLift>(m, "DBCLift")
        .def_readonly("gauss", &DBCLift::gauss)
        .def_readonly("signs", &DBCLift::signs)
        .def_readonly("outside", &DBCLift::outside)
        .def_readonly("inside", &DBCLift::inside)
        .def("__repr__", [](const DBCLift &l) {
            return "<DBCLift gauss=" + std::to_string(l.gauss.size()) + " crossings>";
        });

    m.def("double_branched_cover_lift", &double_branched_cover_lift,
          py::arg("gauss"), py::arg("signs"), py::arg("outer"), ostream_redirect());

    m.def("double_branched_cover_polynomial",
          (Polynomial (*)(const std::vector<int> &, const std::vector<int> &, const std::vector<int> &, time_t))
              &double_branched_cover_polynomial,
          py::arg("gauss"), py::arg("signs"), py::arg("outer"), py::arg("timeout") = 0, ostream_redirect());

    m.def("double_branched_cover_polynomial",
          (Polynomial (*)(PlanarDiagram &, time_t)) &double_branched_cover_polynomial,
          py::arg("diagram"), py::arg("timeout") = 0, ostream_redirect());
}