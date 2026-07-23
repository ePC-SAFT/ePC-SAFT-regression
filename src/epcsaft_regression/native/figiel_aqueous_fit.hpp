#pragma once

#include <Python.h>

namespace epcsaft_regression {

PyObject* evaluate_figiel_aqueous_python(
    PyObject* capsules, PyObject* payload, PyObject* parameters
);
PyObject* solve_figiel_aqueous_python(PyObject* capsules, PyObject* payload);

}  // namespace epcsaft_regression
