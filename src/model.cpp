#include "model.h"
#include "suneigen.h"
#include "exception.h"
#include "defines.h"

#include <iostream>

namespace suneigen {

    Model::Model(ModelDimensions const & model_dimensions)
    : ModelDimensions(model_dimensions), derived_state_(model_dimensions), state_is_non_negative_(false) {
    }

    void Model::initialize(Vector &x, Vector& dx, VectorArray &sx,
                           VectorArray& sdx, bool computeSensitivities) {
        (void) dx; (void)sdx;
        initializeStates(x);
        if (computeSensitivities)
            initializeStateSensitivities(sx, x);
    }

    void Model::initializeStates(Vector& x) {
        std::fill(derived_state_.x_rdata_.begin(), derived_state_.x_rdata_.end(), 0.0);
        fx0(x);
    }

    void Model::fx0(Vector &x) {
        (void)x;
        fx0(derived_state_.x_rdata_.data(), tstart_);
        fx_solver(x.data(), derived_state_.x_rdata_.data());

    }

    void Model::fx_solver(realtype *x_solver, const realtype *x_rdata) {
        std::copy_n(x_rdata, nx, x_solver);
    }

    void Model::initializeStateSensitivities(VectorArray &sx,
                                             Vector const &x) {
        (void) sx;
        (void) x;
    }

    int Model::checkFinite(gsl::span<const realtype> array, const char* fun) const {

        auto result = app->checkFinite(array, fun);
        /*
        if (result != SUNEIGEN_SUCCESS) {
            app->checkFinite(state_.fixedParameters, "k");
            app->checkFinite(state_.unscaledParameters, "p");
            app->checkFinite(derived_state_.w_, "w");
            app->checkFinite(simulation_parameters_.ts_, "t");
        }
        */

        return result;
    }

    void Model::setTimepoints(const std::vector<realtype> &ts) {
        if (!std::is_sorted(ts.begin(), ts.end()))
            throw SunException("Encountered non-monotonic timepoints, please order"
                               " timepoints such that they are monotonically"
                               " increasing!");
        ts_ = ts;
    }

    void Model::fx_rdata(Vector& x_rdata, const Vector& x) {
        std::copy_n(x.data(), nx, x_rdata.data());
        if (always_check_finite_)
            checkFinite(x_rdata.getVector(), "x_rdata");
    }

    void Model::fsx_rdata(VectorArray &sx_rdata, const VectorArray &sx) {
        for (unsigned int ip = 0; ip < np(); ip++) {
            std::copy_n(sx.data(ip), nx, sx_rdata.data(ip));
        }
    }

    const_N_Vector Model::computeX_pos(const_N_Vector x) {
        if (any_state_non_negative_) {
            for (unsigned int ix = 0; ix < x_pos_tmp_.getLength(); ++ix) {
                x_pos_tmp_.at(ix) =
                        (state_is_non_negative_.at(ix) && NV_Ith_S(x, ix) < 0)
                        ? 0
                        : NV_Ith_S(x, ix);
            }
            return x_pos_tmp_.getNVector();
        }
        return x;
    }

    double Model::getTimepoint(const size_t it) const { return ts_.at(it); }

    std::vector<realtype> const &Model::getTimepoints() const { return ts_; }

    double Model::t0() const { return tstart_; }

    size_t Model::nt() const { return ts_.size(); }

    size_t Model::np() const { return static_cast<ModelDimensions const&>(*this).np; }

    void Model::setT0(double t0) { tstart_ = t0; }



}
