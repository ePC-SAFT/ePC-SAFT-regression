#ifndef EPCSAFT_REGRESSION_PURE_SATURATION_FIT_INTERNAL_HPP
#define EPCSAFT_REGRESSION_PURE_SATURATION_FIT_INTERNAL_HPP

#include <Python.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace epcsaft_regression::internal {

constexpr std::size_t residuals_per_row = 4;
constexpr std::size_t baygi_equilibrium_max_iterations = 60;
constexpr std::size_t baygi_equilibrium_max_backtracks = 16;
constexpr double baygi_equilibrium_relative_rank_threshold = 1.0e-12;
constexpr double baygi_smooth_absolute_delta = 1.0e-4;

struct PyObjectDeleter final {
    void operator()(PyObject* object) const noexcept { Py_XDECREF(object); }
};
using OwnedPyObject = std::unique_ptr<PyObject, PyObjectDeleter>;

struct Row final {
    std::string row_id;
    std::string component_id;
    double temperature;
    double pressure;
    double liquid_density;
    std::string source_id;
};

struct Payload final {
    std::vector<std::string> identity;
    std::vector<Row> rows;
    std::vector<double> start;
    std::vector<double> confirmation_start;
    std::vector<double> lower;
    std::vector<double> upper;
    std::vector<double> parameter_scale;
    double amount;
    double molar_mass;
    std::vector<double> weights;
    std::array<double, 2> liquid_volume_bounds;
    std::array<double, 2> vapor_volume_bounds;
    double topology_separation;
    int max_iterations;
    double function_tolerance;
    double gradient_tolerance;
    double parameter_tolerance;
    std::array<double, 2> reporting_pressure_bounds;
    double confirmation_liquid_start_multiplier;
    double confirmation_vapor_start_multiplier;
    double confirmation_parameter_delta;
    double confirmation_cost_delta;
    double reporting_pressure_closure;
    double reporting_mu_closure;
    int num_threads;
};

inline std::size_t parameter_count(const Payload& payload) {
    return payload.start.size();
}
inline std::size_t row_count(const Payload& payload) {
    return payload.rows.size();
}
inline std::size_t variable_count(const Payload& payload) {
    if (payload.identity[1] == "monoethanolamine") {
        return parameter_count(payload);
    }
    return parameter_count(payload) + 2 * row_count(payload);
}
inline std::size_t residual_count(const Payload& payload) {
    if (payload.identity[1] == "monoethanolamine") {
        return 2 * row_count(payload);
    }
    return residuals_per_row * row_count(payload);
}

bool positive_finite(double value);
std::vector<double> doubles(PyObject* object, std::size_t expected, const char* label);
Payload parse_payload(PyObject* object);
Row parse_row(PyObject* object, const std::string& component_id, std::size_t source_index);
std::size_t reporting_row_count(const Payload& payload);

}  // namespace epcsaft_regression::internal

#endif
