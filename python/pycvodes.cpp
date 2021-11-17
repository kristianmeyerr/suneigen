#include "pycvodes.h"
#include "suneigen.h"
#include "py_model_ode.h"
#include <utility>
#include <string>
#include <typeinfo>

namespace py = pybind11;
using namespace pybind11::literals;

py::tuple cvodes(const std::function<py::array_t<double>(double, py::array_t<double>)> &ode,
                 const std::function<py::array_t<double>(double, py::array_t<double>)> &jac,
                 const py::array_t<double>& x0, const py::array_t<double>& times,
                 const py::dict& options) {

    // Process inputs
    py::buffer_info times_buf = times.request();
    if (times_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one for the times input.");
    auto *times_ptr = static_cast<double *>(times_buf.ptr);

    int nt = times.size();
    std::vector<realtype> ts(nt);
    for (int j = 0; j < nt; ++j) {
        ts[j] = times_ptr[j];
    }

    // Create a model with python callbacks.
    auto model = std::unique_ptr<suneigen::Model>(
            new suneigen::pymodel::PyModel_ODE(ode, jac, x0));

    // Set desired output timepoints
    model->setTimepoints(ts);

    // Create a solver instance
    auto solver = model->getSolver();

    // Set options
    for (auto item : options){

        auto option = py::cast<std::string>(item.first);

        // Transform to lower case
        std::transform(option.begin(), option.end(), option.begin(),
                       [](unsigned char c){ return std::tolower(c); });

        if (option == "abstol"){
            auto atol = py::cast<double>(item.second);
            solver->setAbsoluteTolerance(atol);
        } else if (option == "reltol") {
            auto rtol = py::cast<double>(item.second);
            solver->setRelativeTolerance(rtol);
        } else {
            throw std::runtime_error("Option: {} is not supported."_s.format(option));
        }
    }

    // Optionally set integration tolerance



    // Create an application instance
    auto app = suneigen::SunApplication();

    // Run the simulation
    auto rdata = app.runSimulation(*solver,*model);

    // Get a numpy result vector
    auto result = py::array_t<double>(rdata->x.size());
    py::buffer_info buffer_result = result.request();
    auto *ptr_result = static_cast<double *>(buffer_result.ptr);
    for (unsigned int i = 0; i < rdata->x.size(); ++i) {
        ptr_result[i] = rdata->x[i];
    }

    // Reshape into a 2d array (nt x nx)
    int nx = model->nx;
    result.resize({nt, nx});

    return py::make_tuple(times, result);

}