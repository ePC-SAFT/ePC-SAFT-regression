#ifndef EPCSAFT_REGRESSION_BORN_DIAMETER_FIT_HPP
#define EPCSAFT_REGRESSION_BORN_DIAMETER_FIT_HPP

#include <Python.h>

namespace epcsaft_regression {

PyObject* evaluate_born_python(
    PyObject* capsules,
    PyObject* payload,
    PyObject* diameters
);
PyObject* solve_born_python(PyObject* capsules, PyObject* payload);

}  // namespace epcsaft_regression

#endif
