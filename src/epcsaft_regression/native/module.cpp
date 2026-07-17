#include "pure_saturation_fit.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace epcsaft_regression {
namespace {

struct NativeSdkPrefix final {
    std::uint32_t abi_version;
    std::size_t table_size;
    std::size_t result_size;
    const void* model_context;
    epcsaft_evaluate_pure_phase_v1 evaluate_pure_phase;
};

}  // namespace

const epcsaft_native_sdk_v1* checked_provider_table(PyObject* capsule) {
    if (!PyCapsule_CheckExact(capsule)) {
        PyErr_SetString(PyExc_TypeError, "provider transport must be an exact CPython capsule");
        return nullptr;
    }
    void* pointer = PyCapsule_GetPointer(capsule, EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME);
    if (pointer == nullptr) {
        return nullptr;
    }
    NativeSdkPrefix prefix{};
    std::memcpy(&prefix, pointer, sizeof(prefix));
    if (prefix.abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION) {
        PyErr_SetString(PyExc_RuntimeError, "provider native SDK ABI version mismatch");
        return nullptr;
    }
    constexpr std::size_t minimum_size =
        offsetof(epcsaft_native_sdk_v1, evaluate_pure_phase_parameters)
        + sizeof(epcsaft_evaluate_pure_phase_parameters_v1);
    if (prefix.table_size < minimum_size) {
        PyErr_SetString(PyExc_RuntimeError, "provider native SDK table lacks the parameterized tail");
        return nullptr;
    }
    const auto* table = static_cast<const epcsaft_native_sdk_v1*>(pointer);
    if (table->parameterized_result_size != sizeof(epcsaft_parameterized_phase_block_result_v1)) {
        PyErr_SetString(PyExc_RuntimeError, "provider parameterized result size mismatch");
        return nullptr;
    }
    if (table->model_context == nullptr || table->evaluate_pure_phase_parameters == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "provider parameterized evaluator is unavailable");
        return nullptr;
    }
    return table;
}

Py_ssize_t minimum_parameterized_table_size() noexcept {
    return static_cast<Py_ssize_t>(
        offsetof(epcsaft_native_sdk_v1, evaluate_pure_phase_parameters)
        + sizeof(epcsaft_evaluate_pure_phase_parameters_v1)
    );
}

Py_ssize_t parameterized_result_size() noexcept {
    return static_cast<Py_ssize_t>(sizeof(epcsaft_parameterized_phase_block_result_v1));
}

PyObject* transport_info(PyObject* capsule) {
    const epcsaft_native_sdk_v1* table = checked_provider_table(capsule);
    if (table == nullptr) {
        return nullptr;
    }
    epcsaft_parameterized_phase_block_result_v1 result{};
    result.struct_size = sizeof(result);
    const int status = table->evaluate_pure_phase_parameters(
        table->model_context,
        130.0,
        1.0,
        1.0e-3,
        1.0,
        3.7039,
        150.03,
        &result
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1 || result.status != status) {
        PyErr_Format(
            PyExc_RuntimeError,
            "provider parameterized evaluator rejected the transport probe: %s",
            result.error
        );
        return nullptr;
    }
    const std::size_t fingerprint_length = strnlen(
        result.parameter_fingerprint, EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE
    );
    return Py_BuildValue(
        "(knns#)",
        static_cast<unsigned long>(table->abi_version),
        static_cast<Py_ssize_t>(table->table_size),
        static_cast<Py_ssize_t>(table->parameterized_result_size),
        result.parameter_fingerprint,
        static_cast<Py_ssize_t>(fingerprint_length)
    );
}

}  // namespace epcsaft_regression

namespace {

PyObject* py_transport_info(PyObject*, PyObject* capsule) {
    return epcsaft_regression::transport_info(capsule);
}

PyObject* py_minimum_parameterized_table_size(PyObject*, PyObject*) {
    return PyLong_FromSsize_t(epcsaft_regression::minimum_parameterized_table_size());
}

PyObject* py_parameterized_result_size(PyObject*, PyObject*) {
    return PyLong_FromSsize_t(epcsaft_regression::parameterized_result_size());
}

PyObject* py_evaluate(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    PyObject* variables = nullptr;
    if (!PyArg_ParseTuple(args, "OOO:evaluate", &capsule, &payload, &variables)) {
        return nullptr;
    }
    return epcsaft_regression::evaluate_python(capsule, payload, variables);
}

PyObject* py_solve(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    PyObject* reporting_rows = nullptr;
    if (!PyArg_ParseTuple(args, "OOO:solve", &capsule, &payload, &reporting_rows)) {
        return nullptr;
    }
    return epcsaft_regression::solve_python(capsule, payload, reporting_rows);
}

PyMethodDef methods[] = {
    {"transport_info", py_transport_info, METH_O, "Validate the installed provider capsule."},
    {
        "minimum_parameterized_table_size",
        py_minimum_parameterized_table_size,
        METH_NOARGS,
        "Return the minimum v1 parameterized table size."
    },
    {
        "parameterized_result_size",
        py_parameterized_result_size,
        METH_NOARGS,
        "Return the required parameterized result size."
    },
    {"evaluate", py_evaluate, METH_VARARGS, "Evaluate exact pure-saturation residuals and Jacobian."},
    {"solve", py_solve, METH_VARARGS, "Fit pure-saturation parameters and evaluate reporting rows."},
    {nullptr, nullptr, 0, nullptr},
};

PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_native",
    "Native Ceres pure-saturation regression transport.",
    -1,
    methods,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
};

}  // namespace

PyMODINIT_FUNC PyInit__native() {
    return PyModule_Create(&module);
}
