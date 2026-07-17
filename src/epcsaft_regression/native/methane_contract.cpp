#include "methane_fit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_regression::internal {
namespace {

std::string text(PyObject* object, const char* label) {
    if (!PyUnicode_Check(object)) {
        throw std::invalid_argument(std::string(label) + " must be text");
    }
    Py_ssize_t size = 0;
    const char* value = PyUnicode_AsUTF8AndSize(object, &size);
    if (value == nullptr) {
        throw std::invalid_argument(std::string(label) + " must be UTF-8 text");
    }
    return std::string(value, static_cast<std::size_t>(size));
}

std::vector<std::string> texts(PyObject* object, std::size_t expected, const char* label) {
    PyObject* sequence = PySequence_Fast(object, label);
    if (sequence == nullptr) throw std::invalid_argument(label);
    if (PySequence_Fast_GET_SIZE(sequence) != static_cast<Py_ssize_t>(expected)) {
        Py_DECREF(sequence);
        throw std::invalid_argument(std::string(label) + " has the wrong length");
    }
    std::vector<std::string> values;
    values.reserve(expected);
    for (std::size_t index = 0; index < expected; ++index) {
        values.push_back(text(
            PySequence_Fast_GET_ITEM(sequence, static_cast<Py_ssize_t>(index)), label
        ));
    }
    Py_DECREF(sequence);
    return values;
}

template <std::size_t Size>
std::array<double, Size> fixed_doubles(PyObject* object, const char* label) {
    const std::vector<double> parsed = doubles(object, Size, label);
    std::array<double, Size> values{};
    std::copy(parsed.begin(), parsed.end(), values.begin());
    return values;
}

constexpr std::array<const char*, 9> expected_row_ids{
    "nist-methane-sat-100-k", "nist-methane-sat-110-k", "nist-methane-sat-120-k",
    "nist-methane-sat-130-k", "nist-methane-sat-140-k", "nist-methane-sat-150-k",
    "nist-methane-sat-160-k", "nist-methane-sat-170-k", "nist-methane-sat-180-k",
};
constexpr std::array<double, 9> expected_pressures{
    34375.892, 88130.038, 191430.08, 367319.94, 641181.43,
    1039961.3, 1592078.0, 2328348.8, 3285180.7,
};
constexpr std::array<double, 9> expected_liquid_densities{
    438.88524, 424.77725, 409.90234, 394.03734, 376.86505,
    357.89846, 336.31495, 310.50203, 276.22850,
};

std::vector<std::string> expected_identity() {
    return {
        "nist-webbook-methane-saturation-100-180-k-v1",
        "methane", "K", "Pa", "kg/m3",
        "nist-webbook-srd69-methane-saturation",
        "NIST Chemistry WebBook, SRD 69, methane (CAS 74-82-8) fluid properties",
        "Saturation properties query, 100 K through 180 K in 10 K increments",
        "https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74828&Type=SatP&Digits=8&THigh=180&TLow=100&TInc=10&RefState=DEF&TUnit=K&PUnit=Pa&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm",
        "2026-07-17",
        "NIST Standard Reference Data retained as compact source-backed candidate evidence; redistribution and use remain subject to the NIST SRD terms",
        "Exact retained CSV fields and decimal strings; CRLF line endings normalized to LF. Retained source SHA-256: a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227; packaged SHA-256: dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c.",
        "temperature", "K", "pressure", "Pa", "saturated_liquid_mass_density", "kg/m3",
        "a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227",
        "dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c",
        "pure-methane-saturation-lifted-volumes-v1",
        "segment_count", "segment_diameter_angstrom", "dispersion_energy_over_k_kelvin",
        "1", "angstrom", "K",
        "liquid_pressure", "vapor_pressure", "chemical_potential_equality", "liquid_density",
        "mol", "m3/mol", "epcsaft.native_sdk.v1",
        "4b10cb899c94687cae734980285badb224dc95e6",
        "f92f79c8d6f614660e5c201b7061c9b02b5cd1a25a4ed8c8fee0b59adaabf2bf",
        "sha256:5f836aa84935df70be2e5cffae51b178a7b797c2cee036e9ff47d8097ca94bbf",
        "p_j = start_j + parameter_scale_j * z_j",
        "V_liquid = (molar_mass / observed_liquid_density) * exp(u_liquid)",
        "V_vapor = (R * T / observed_pressure) * exp(u_vapor)",
        "P_report = observed_pressure * exp(u_pressure)",
        "DENSE_QR", "SILENT",
    };
}

}  // namespace

bool positive_finite(double value) {
    return std::isfinite(value) && value > 0.0;
}

std::vector<double> doubles(PyObject* object, std::size_t expected, const char* label) {
    PyObject* sequence = PySequence_Fast(object, label);
    if (sequence == nullptr) throw std::invalid_argument(label);
    const Py_ssize_t size = PySequence_Fast_GET_SIZE(sequence);
    if (size != static_cast<Py_ssize_t>(expected)) {
        Py_DECREF(sequence);
        throw std::invalid_argument(std::string(label) + " has the wrong length");
    }
    std::vector<double> values;
    values.reserve(expected);
    for (Py_ssize_t index = 0; index < size; ++index) {
        const double value = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(sequence, index));
        if (PyErr_Occurred() != nullptr) {
            Py_DECREF(sequence);
            throw std::invalid_argument(std::string(label) + " must contain numbers");
        }
        values.push_back(value);
    }
    Py_DECREF(sequence);
    return values;
}

Row parse_row(PyObject* object, std::size_t source_index) {
    PyObject* sequence = PySequence_Fast(object, "source row must be a sequence");
    if (sequence == nullptr || PySequence_Fast_GET_SIZE(sequence) != 6) {
        Py_XDECREF(sequence);
        throw std::invalid_argument("source row must contain exactly six fields");
    }
    PyObject** items = PySequence_Fast_ITEMS(sequence);
    Row row{
        text(items[0], "row_id"),
        text(items[1], "species"),
        PyFloat_AsDouble(items[2]),
        PyFloat_AsDouble(items[3]),
        PyFloat_AsDouble(items[4]),
        text(items[5], "source_id"),
    };
    Py_DECREF(sequence);
    const double expected_temperature = 100.0 + 10.0 * static_cast<double>(source_index);
    if (PyErr_Occurred() != nullptr
        || row.row_id != expected_row_ids[source_index]
        || row.species != "methane"
        || row.temperature != expected_temperature
        || row.pressure != expected_pressures[source_index]
        || row.liquid_density != expected_liquid_densities[source_index]
        || row.source_id != "nist-webbook-srd69-methane-saturation") {
        throw std::invalid_argument("source row violates the exact retained identity and values");
    }
    return row;
}

Payload parse_payload(PyObject* object) {
    PyObject* sequence = PySequence_Fast(object, "methane fit payload must be a sequence");
    if (sequence == nullptr) throw std::invalid_argument("methane fit payload must be a sequence");
    if (PySequence_Fast_GET_SIZE(sequence) != 24) {
        Py_DECREF(sequence);
        throw std::invalid_argument("methane fit payload must contain exactly 24 fields");
    }
    PyObject** items = PySequence_Fast_ITEMS(sequence);
    Payload payload{};
    payload.identity = texts(items[0], 43, "compiled problem identity");
    if (payload.identity != expected_identity()) {
        Py_DECREF(sequence);
        throw std::invalid_argument("compiled problem identity does not match the first slice");
    }
    PyObject* rows = PySequence_Fast(items[1], "training rows must be a sequence");
    if (rows == nullptr || PySequence_Fast_GET_SIZE(rows) != static_cast<Py_ssize_t>(row_count)) {
        Py_XDECREF(rows);
        Py_DECREF(sequence);
        throw std::invalid_argument("training rows must contain exactly four rows");
    }
    constexpr std::array<std::size_t, row_count> source_indices{1, 3, 5, 7};
    for (std::size_t index = 0; index < row_count; ++index) {
        payload.rows[index] = parse_row(
            PySequence_Fast_GET_ITEM(rows, static_cast<Py_ssize_t>(index)), source_indices[index]
        );
    }
    Py_DECREF(rows);
    payload.start = fixed_doubles<3>(items[2], "parameter start");
    payload.lower = fixed_doubles<3>(items[3], "parameter lower bounds");
    payload.upper = fixed_doubles<3>(items[4], "parameter upper bounds");
    payload.parameter_scale = fixed_doubles<3>(items[5], "parameter scales");
    payload.amount = PyFloat_AsDouble(items[6]);
    payload.molar_mass = PyFloat_AsDouble(items[7]);
    payload.weights = fixed_doubles<4>(items[8], "residual weights");
    payload.liquid_volume_bounds = fixed_doubles<2>(items[9], "liquid volume bounds");
    payload.vapor_volume_bounds = fixed_doubles<2>(items[10], "vapor volume bounds");
    payload.topology_separation = PyFloat_AsDouble(items[11]);
    const long max_iterations = PyLong_AsLong(items[12]);
    payload.function_tolerance = PyFloat_AsDouble(items[13]);
    payload.gradient_tolerance = PyFloat_AsDouble(items[14]);
    payload.parameter_tolerance = PyFloat_AsDouble(items[15]);
    payload.reporting_pressure_bounds = fixed_doubles<2>(items[16], "reporting pressure bounds");
    payload.confirmation_liquid_start_multiplier = PyFloat_AsDouble(items[17]);
    payload.confirmation_vapor_start_multiplier = PyFloat_AsDouble(items[18]);
    payload.confirmation_parameter_delta = PyFloat_AsDouble(items[19]);
    payload.confirmation_cost_delta = PyFloat_AsDouble(items[20]);
    payload.reporting_pressure_closure = PyFloat_AsDouble(items[21]);
    payload.reporting_mu_closure = PyFloat_AsDouble(items[22]);
    const long num_threads = PyLong_AsLong(items[23]);
    Py_DECREF(sequence);
    if (PyErr_Occurred() != nullptr) {
        throw std::invalid_argument("methane fit payload contains a nonnumeric scalar");
    }
    payload.max_iterations = static_cast<int>(max_iterations);
    payload.num_threads = static_cast<int>(num_threads);
    if (payload.start != std::array<double, 3>{1.08, 3.555744, 157.5315}
        || payload.lower != std::array<double, 3>{0.5, 2.0, 50.0}
        || payload.upper != std::array<double, 3>{3.5, 5.0, 400.0}
        || payload.parameter_scale != std::array<double, 3>{0.1, 0.1, 10.0}
        || payload.amount != 1.0 || payload.molar_mass != 0.016043
        || payload.weights != std::array<double, 4>{0.25, 0.25, 0.25, 0.25}
        || payload.liquid_volume_bounds != std::array<double, 2>{2.0e-5, 1.0e-4}
        || payload.vapor_volume_bounds != std::array<double, 2>{1.5e-4, 0.1}
        || payload.topology_separation != 1.0e-3
        || payload.max_iterations != 500
        || payload.function_tolerance != 1.0e-10
        || payload.gradient_tolerance != 1.0e-10
        || payload.parameter_tolerance != 1.0e-10
        || payload.reporting_pressure_bounds != std::array<double, 2>{1.0e3, 1.0e7}
        || payload.confirmation_liquid_start_multiplier != 1.01
        || payload.confirmation_vapor_start_multiplier != 0.98
        || payload.confirmation_parameter_delta != 1.0e-5
        || payload.confirmation_cost_delta != 1.0e-8
        || payload.reporting_pressure_closure != 1.0e-8
        || payload.reporting_mu_closure != 1.0e-8
        || payload.num_threads != 1) {
        throw std::invalid_argument("methane fit payload does not match the first slice contract");
    }
    return payload;
}

}  // namespace epcsaft_regression::internal
