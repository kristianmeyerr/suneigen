#include "suneigen.h"
#include "py_model_ode.h"
#include <utility>
#include <string>
#include <typeinfo>

namespace suneigen {

    namespace py = pybind11;
    using namespace pybind11::literals;

    py::tuple cvodes(const pyfunction_cb &ode,
                     const pyjac_cb &jac,
                     const py::array_t<double>& x0,
                     const py::array_t<double>& times,
                     const int& nnz,
                     const pybind11::array_t<double>& p,
                     const pybind11::array_t<double>& k,
                     const py::dict& options) {

        py::buffer_info times_buf = times.request();
        if (times_buf.ndim != 1)
            throw std::runtime_error("Number of dimensions must be one for the times input.");
        auto times_ptr = static_cast<double*>(times_buf.ptr);
        std::vector<double> times_(times_ptr, times_ptr + times.size());

        py::buffer_info p_buf = p.request();
        if (p_buf.ndim != 1)
            throw std::runtime_error("Number of dimensions must be one for the p input.");
        auto p_ptr = static_cast<double*>(p_buf.ptr);
        std::vector<double> p_(p_ptr, p_ptr + p.size());

        py::buffer_info k_buf = k.request();
        if (k_buf.ndim != 1)
            throw std::runtime_error("Number of dimensions must be one for the k input.");
        auto k_ptr = static_cast<double*>(k_buf.ptr);
        std::vector<double> k_(k_ptr, k_ptr + k.size());

        // Size of problem
        size_t nx = x0.size();
        size_t np = p_.size();
        size_t nk = k_.size();
        size_t ny = 0;
        size_t nz = 0;
        size_t ne = 0;
        size_t ubw = 0;
        size_t lbw = 0;

        // Create a model
        auto model = std::unique_ptr<suneigen::Model>(
                new pymodel::PyModel_ODE(
                        nx,
                        np,
                        nk,
                        ny,
                        nz,
                        ne,
                        nnz,
                        ubw,
                        lbw,
                        ode,
                        jac,
                        x0));

        // Set desired output timepoints
        model->setTimepoints(times_);

        // Set parameters
        model->setParameters(p_);

        // Set fixed parameters
        model->setFixedParameters(k_);

        // Create a solver instance
        auto solver = model->getSolver();

        // Set options
        for (auto item : options) {

            auto option = py::cast<std::string>(item.first);

            // Transform to lower case
            std::transform(option.begin(), option.end(), option.begin(),
                           [](unsigned char c) { return std::tolower(c); });

            if (option == "abstol") {
                auto atol = py::cast<double>(item.second);
                solver->setAbsoluteTolerance(atol);
            } else if (option == "reltol") {
                auto rtol = py::cast<double>(item.second);
                solver->setRelativeTolerance(rtol);
            } else if (option == "linsol") {
                auto ls = py::cast<std::string>(item.second);
                // Transform to lower case
                std::transform(ls.begin(), ls.end(), ls.begin(),
                               [](unsigned char c) { return std::tolower(c); });
                if (ls == "dense")
                    solver->setLinearSolver(suneigen::LinearSolver::dense);
                else if (ls == "superlu")
                    solver->setLinearSolver(suneigen::LinearSolver::SuperLU);
                else {
                    throw std::runtime_error("Linear solver: {} is not supported. Use dense or superlu"_s.format(ls));
                }
            } else {
                throw std::runtime_error("Option: {} is not supported."_s.format(option));
            }
        }

        // Create an application instance
        auto app = suneigen::SunApplication();

        // Run the simulation
        auto rdata = app.runSimulation(*solver, *model);

        // Get a numpy result vector
        py::array_t<double> result = py::array_t<double>(static_cast<ssize_t>(rdata->x.size()), rdata->x.data());

        // Reshape into a 2d array (nt x nx)
        size_t nt = times_.size();
        result.resize({nt, nx});

        rdata->printStatistics();

        return py::make_tuple(times, result);

    }

}