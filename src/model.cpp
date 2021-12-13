#include "model.h"
#include "suneigen.h"
#include "exception.h"
#include "defines.h"

#include <iostream>

namespace suneigen {

    Model::Model(ModelDimensions const & model_dimensions,
                 SimulationParameters simulation_parameters, std::vector<int> z2event)
    : ModelDimensions(model_dimensions),
        derived_state_(model_dimensions),
        z2event_(std::move(z2event)),
        state_is_non_negative_(nx,false),
        simulation_parameters_(std::move(simulation_parameters))
    {
        Expects(model_dimensions.np == simulation_parameters_.parameters.size());
        Expects(model_dimensions.nk == simulation_parameters_.fixedParameters.size());

        state_.h.resize(ne, 0.0);
        state_.x.resize(nx, 0.0);
        state_.parameters = simulation_parameters_.parameters;
        state_.fixedParameters = simulation_parameters_.fixedParameters;
    }

    void Model::initialize(Vector &x, Vector& dx, VectorArray &sx,
                           VectorArray& /*sdx*/, bool computeSensitivities) {

        initializeStates(x);
        if (computeSensitivities)
            initializeStateSensitivities(sx, x);

        if (ne){
            initHeaviside(x, dx);
        }
    }

    void Model::initializeStates(Vector& x) {

        if (x0data_.empty()) {
            fx0(x);
        } else {
            std::vector<realtype> x0_solver(nx, 0.0);
            fx_solver(x0_solver.data(), x0data_.data());
            std::copy(x0_solver.cbegin(), x0_solver.cend(), x.data());
        }
    }

    void Model::initializeStateSensitivities(VectorArray &sx,
                                             Vector const &x) {
        if (sx0data_.empty()) {
            fsx0(sx, x);
        } else {
            std::vector<realtype> sx0_solver_slice(nx, 0.0);
            for (unsigned int ip = 0; ip < np(); ip++) {
                fsx_solver(sx0_solver_slice.data(), &sx0data_.at(ip * nx));
                for (unsigned int ix = 0; ix < nx; ix++) {
                    sx.at(ix, ip) = sx0_solver_slice.at(ix);
                }
            }
        }
    }

    void Model::initHeaviside(Vector const &x, Vector const &dx) {
        std::vector<realtype> rootvals(ne, 0.0);
        froot(simulation_parameters_.tstart_, x, dx, rootvals);
        for (unsigned int ie = 0; ie < ne; ie++) {
            if (rootvals.at(ie) < 0) {
                state_.h.at(ie) = 0.0;
            } else {
                state_.h.at(ie) = 1.0;
            }
        }
    }

    size_t Model::np() const { return static_cast<ModelDimensions const&>(*this).np; }

    size_t Model::nk() const { return state_.fixedParameters.size(); }

    unsigned int Model::nMaxEvent() const { return nmaxevent_; }

    void Model::setNMaxEvent(unsigned int nmaxevent) { nmaxevent_ = nmaxevent; }

    size_t Model::nt() const { return simulation_parameters_.ts_.size(); }

    std::vector<realtype> const &Model::getParameters() const {
        return simulation_parameters_.parameters;
    }

    void Model::getEvent(gsl::span<realtype> z, const unsigned int ie, const realtype t,
                         const Vector &x) {
        fz(ie, t, x);
        writeSliceEvent(derived_state_.z_, z, static_cast<int>(ie));
    }

    void Model::fz(const unsigned int ie, const realtype t, const Vector &x) {

        derived_state_.z_.assign(nz, 0.0);

        fz(derived_state_.z_.data(), ie, t, x.data(), state_.parameters.data(),
                state_.fixedParameters.data(), state_.h.data());
    }

    void Model::setParameters(const std::vector<realtype> &p) {
        if (p.size() != np())
            throw SunException("Dimension mismatch. Size of parameters does not "
                               "match number of model parameters.");
        simulation_parameters_.parameters = p;
        state_.parameters.resize(simulation_parameters_.parameters.size());

        for (gsl::span<realtype>::index_type ip = 0;
             ip < simulation_parameters_.parameters.size(); ++ip) {
            state_.parameters[ip] = simulation_parameters_.parameters[ip];
        }
    }

    void Model::setFixedParameters(const std::vector<realtype> &k) {
        if (k.size() != nk())
            throw SunException("Dimension mismatch. Size of fixedParameters does "
                               "not match number of fixed model parameters.");
        state_.fixedParameters = k;
    }

    std::vector<realtype> const &Model::getTimepoints() const { return simulation_parameters_.ts_; }

    void Model::setTimepoints(const std::vector<realtype> &ts) {
        if (!std::is_sorted(ts.begin(), ts.end()))
            throw SunException("Encountered non-monotonic timepoints, please order"
                               " timepoints such that they are monotonically"
                               " increasing!");
        simulation_parameters_.ts_ = ts;
    }

    double Model::getTimepoint(const size_t it) const { return simulation_parameters_.ts_.at(it); }

    double Model::t0() const { return simulation_parameters_.tstart_; }

    void Model::setT0(double t0) { simulation_parameters_.tstart_ = t0; }

    std::vector<bool> const &Model::getStateIsNonNegative() const {
        return state_is_non_negative_;
    }

    void Model::setStateIsNonNegative(std::vector<bool> const &nonNegative) {

        if (state_is_non_negative_.size() != static_cast<unsigned int>(nx)) {
            throw SunException("Dimension of input stateIsNonNegative (%u) does "
                               "not agree with number of state variables (%d)",
                               state_is_non_negative_.size(), nx);
        }
        state_is_non_negative_ = nonNegative;
        any_state_non_negative_ =
                std::any_of(state_is_non_negative_.begin(), state_is_non_negative_.end(),
                            [](bool x) { return x; });
    }

    void Model::setAllStatesNonNegative() {
        setStateIsNonNegative(std::vector<bool>(nx, true));
    }

    std::vector<realtype> Model::getInitialStates() {
        if(!x0data_.empty()) {
            return x0data_;
        }

        /* Initial states have not been set explicitly on this instance, so we
         * compute it, but don't save it, as this would have to be invalidated upon
         * changing parameters etc.
         */
        std::vector<realtype> x0(nx, 0.0);
        fx0(x0.data(), simulation_parameters_.tstart_,
                state_.parameters.data(), state_.fixedParameters.data());
        return x0;
    }

    void Model::setInitialStates(const std::vector<realtype> &x0) {
        if (x0.size() != nx && !x0.empty())
            throw SunException("Dimension mismatch. Size of x0 does not match "
                               "number of model states.");

        if (x0.empty()) {
            x0data_.clear();
            return;
        }

        x0data_ = x0;
    }

    std::vector<realtype> Model::getInitialStateSensitivities() {
        if(!sx0data_.empty()) {
            return sx0data_;
        }

        /* Initial state sensitivities have not been set explicitly on this
         * instance, so we compute it, but don't save it, as this would have to be
         * invalidated upon changing parameters etc.
         */
        std::vector<realtype> sx0(nx * np(), 0.0);
        auto x0 = getInitialStates();
        for (unsigned int ip = 0; ip < np(); ip++) {
            fsx0(sx0.data(), simulation_parameters_.tstart_, x0.data(), state_.parameters.data(), ip);
        }
        return sx0;
    }

    void Model::setInitialStateSensitivities(const std::vector<realtype> &sx0) {

        if (sx0.size() != nx * np() && !sx0.empty())
            throw SunException("Dimension mismatch. Size of sx0 does not match "
                               "number of model states * number of parameter "
                               "selected for sensitivities.");

        if (sx0.empty()) {
            sx0data_.clear();
            return;
        }

        std::vector<realtype> sx0_rdata(nx * np(), 0.0);
        for (unsigned int ip = 0; ip < np(); ip++) {

            for (unsigned int ix = 0; ix < nx; ++ix) {
                sx0_rdata.at(ip * nx + ix) =
                        sx0.at(ip * nx + ix);
            }
        }
        setUnscaledInitialStateSensitivities(sx0_rdata);
    }

    void Model::setUnscaledInitialStateSensitivities(
            const std::vector<realtype> &sx0) {
        if (sx0.size() != nx * np() && !sx0.empty())
            throw SunException("Dimension mismatch. Size of sx0 does not match "
                               "number of model states * number of parameter "
                               "selected for sensitivities.");

        if (sx0.empty()) {
            sx0data_.clear();
            return;
        }

        sx0data_ = sx0;
    }


    void Model::updateHeaviside(const std::vector<int> &rootsfound) {
        for (unsigned int ie = 0; ie < ne; ie++) {
            state_.h.at(ie) += rootsfound.at(ie);
        }
    }

    void Model::fx0(Vector &x) {

        std::fill(derived_state_.x_rdata_.begin(), derived_state_.x_rdata_.end(), 0.0);

        fx0(derived_state_.x_rdata_.data(), simulation_parameters_.tstart_,
                state_.parameters.data(), state_.fixedParameters.data());
        fx_solver(x.data(), derived_state_.x_rdata_.data());

        if (always_check_finite_) {
            checkFinite(derived_state_.x_rdata_, "x0 x_rdata");
            checkFinite(x.getVector(), "x0 x");
        }
    }

    void Model::fsx0(VectorArray &sx, const Vector &x) {

        for (unsigned int ip = 0; ip < np(); ip++) {
            std::fill(derived_state_.sx_rdata_.begin(),
                      derived_state_.sx_rdata_.end(), 0.0);
            fsx0(derived_state_.sx_rdata_.data(), simulation_parameters_.tstart_,
                 x.data(), state_.parameters.data(), ip);
            fsx_solver(sx.data(ip), derived_state_.sx_rdata_.data());
        }
    }

    void Model::fx_solver(realtype *x_solver, const realtype *x_rdata) {
        std::copy_n(x_rdata, nx, x_solver);
    }

    void Model::fsx_solver(realtype *sx_solver, const realtype *sx_rdata) {
        /* for the moment we do not need an implementation of fsx_solver as
         * we can simply reuse fx_solver and replace states by their
         * sensitivities */
        fx_solver(sx_solver, sx_rdata);
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

    void Model::writeSliceEvent(gsl::span<const realtype> slice,
                                gsl::span<realtype> buffer, const int ie) {
        checkBufferSize(buffer, slice.size());
        checkBufferSize(buffer, z2event_.size());
        for (unsigned izt = 0; izt < z2event_.size(); ++izt)
            if (z2event_.at(izt) - 1 == ie)
                buffer[izt] = slice[izt];
    }

    void Model::writeSensitivitySliceEvent(gsl::span<const realtype> slice,
                                           gsl::span<realtype> buffer,
                                           const int ie) {
        checkBufferSize(buffer, slice.size());
        checkBufferSize(buffer, z2event_.size() * np());
        for (unsigned int ip = 0; ip < np(); ++ip)
            for (unsigned izt = 0; izt < z2event_.size(); ++izt)
                if (z2event_.at(izt) - 1 == ie)
                    buffer[ip * ne + izt] = slice[ip * ne + izt];
    }

    const_N_Vector Model::computeX_pos(const_N_Vector x) {
        if (any_state_non_negative_) {
            for (unsigned int ix = 0; ix < derived_state_.x_pos_tmp_.getLength(); ++ix) {
                derived_state_.x_pos_tmp_.at(ix) =
                        (state_is_non_negative_.at(ix) && NV_Ith_S(x, ix) < 0)
                        ? 0
                        : NV_Ith_S(x, ix);
            }
            return derived_state_.x_pos_tmp_.getNVector();
        }
        return x;
    }

    void Model::getEventTimeSensitivity(std::vector<realtype> &stau,
                                        const realtype t, const unsigned int ie,
                                        const Vector &x,
                                        const VectorArray &sx) {

        std::fill(stau.begin(), stau.end(), 0.0);
        (void)t; (void)ie; (void)x;(void)sx;

        for (unsigned int ip = 0; ip < np(); ip++) {
            fstau(&stau.at(ip), t, x.data(), state_.parameters.data(), state_.h.data(), sx.data(ip),
                  ip, ie);
        }
    }

    void Model::addStateEventUpdate(Vector &x, const unsigned int ie, const realtype t,
                                    const Vector &xdot,
                                    const Vector &xdot_old) {

        derived_state_.deltax_.assign(nx, 0.0);

        // compute update
        fdeltax(derived_state_.deltax_.data(), t, x.data(), state_.parameters.data(),
                state_.h.data(), ie, xdot.data(),
                xdot_old.data());

        if (always_check_finite_) {
            SunApplication app;
            app.checkFinite(derived_state_.deltax_, "deltax");
        }

        // update @todo Make this an eigen version
        // amici_daxpy(nx, 1.0, derived_state_.deltax_.data(), 1, x.data(), 1);
    }

    void Model::addStateSensitivityEventUpdate(VectorArray &sx, const unsigned int ie,
                                               const realtype t,
                                               const Vector &x_old,
                                               const Vector &xdot,
                                               const Vector &xdot_old,
                                               const std::vector<realtype> &stau) {
        (void)sx;
        (void)ie;
        (void)t;
        (void)x_old;
        (void)xdot;
        (void)xdot_old;
        (void)stau;
    }

    int Model::checkFinite(gsl::span<const realtype> array, const char* fun) const {
        SunApplication app;
        auto result = app.checkFinite(array, fun);
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

    void Model::setAlwaysCheckFinite(bool alwaysCheck) {
        always_check_finite_ = alwaysCheck;
    }

    bool Model::getAlwaysCheckFinite() const { return always_check_finite_; }



}
