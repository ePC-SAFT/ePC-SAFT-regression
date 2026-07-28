#ifndef EPCSAFT_REGRESSION_GENERAL_FIT_HPP
#define EPCSAFT_REGRESSION_GENERAL_FIT_HPP

#include <Python.h>

namespace epcsaft_regression {

PyObject* parameter_capabilities_python(PyObject* capsule);
PyObject* evaluate_general_pair_python(
    PyObject* capsule, PyObject* payload, PyObject* variables
);
PyObject* solve_general_pair_python(PyObject* capsule, PyObject* payload);

}  // namespace epcsaft_regression

#endif
