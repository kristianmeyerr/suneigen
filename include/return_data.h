#ifndef SUNEIGEN_RETURN_DATA_H
#define SUNEIGEN_RETURN_DATA_H

#include "defines.h"
#include "vector.h"
#include "model.h"
#include "misc.h"
#include "forward_problem.h"

#include <vector>

namespace suneigen {

    class ReturnData;
    class Solver;
    class ForwardProblem;
    struct SimulationState;

    /**
     * @brief Stores all data to be returned by amici::runAmiciSimulation.
     *
     * NOTE: multi-dimensional arrays are stored in row-major order (C-style)
     */
    class ReturnData: public ModelDimensions {
    public:
        /**
         * @brief Default constructor
         */
        ReturnData() = default;

        /**
         * @brief constructor that uses information from model and solver to
         * appropriately initialize fields
         * @param solver solver instance
         * @param model model instance
         */
        ReturnData(Solver const& solver, const Model& model);

        ReturnData(std::vector<realtype> ts,
                   ModelDimensions const& model_dimensions,
                   unsigned int nmaxevent, size_t nt, SensitivityOrder sensi, SensitivityMethod sensi_meth
                   );

        ~ReturnData() = default;

        /**
         * @brief constructor that uses information from model and solver to
         * appropriately initialize fields
         * @param fwd simulated forward problem, pass `nullptr` to ignore
         * @param model matching model instance
         * @param solver matching solver instance
         */
        void processSimulationObjects(ForwardProblem const *fwd,
                                      Model &model, Solver const &solver);

        /**
         * @brief Print final simulation statistics
         */
        void printStatistics();

        /**
         * timepoints (shape `nt`)
         */
        std::vector<realtype> ts;

        /** time derivative (shape `nx`) */
        std::vector<realtype> xdot;

        /**
         * Jacobian of differential equation right hand side (shape `nx` x `nx`, row-major)
         */
        std::vector<realtype> J;

        /** event output (shape `nmaxevent` x `nz`, row-major) */
        std::vector<realtype> z;

        /**
         * parameter derivative of event output (shape `nmaxevent` x `nz`, row-major)
         */
        std::vector<realtype> sz;

        /** event trigger output (shape `nmaxevent` x `nz`, row-major)*/
        std::vector<realtype> rz;

        /**
         * parameter derivative of event trigger output (shape `nmaxevent` x `nz` x `nplist`, row-major)
         */
        std::vector<realtype> srz;

        /** state (shape `nt` x `nx`, row-major) */
        std::vector<realtype> x;

        /**
         * parameter derivative of state (shape `nt` x `nplist` x `nx`, row-major)
         */
        std::vector<realtype> sx;

        /** number of integration steps forward problem (shape `nt`) */
        std::vector<int> numsteps;

        /** number of right hand side evaluations forward problem (shape `nt`) */
        std::vector<int> numrhsevals;

        /** number of error test failures forward problem (shape `nt`) */
        std::vector<int> numerrtestfails;

        /**
         * number of non-linear solver convergence failures forward problem (shape `nt`)
         */
        std::vector<int> numnonlinsolvconvfails;

        /**
         * The number of non-linear solver iterations (shape `nt`)
         */
         std::vector<int> numnonlinsolviter;

        /**
        * The number of jacobian evaluations (shape `nt`)
        */
        std::vector<int> numjacevals;

        /** employed order forward problem (shape `nt`) */
        std::vector<int> order;

        /** computation time of forward solve [ms] */
        double cpu_time = 0.0;

        /** The initial simulation time point */
        double t0 = std::nan("");

        /** The absolute tolerance */
        realtype abstol = std::nan("");

        /** The relative tolerance */
        realtype reltol = std::nan("");

        /** initial state (shape `nx`) */
        std::vector<realtype> x0;

        /** initial sensitivities (shape `nplist` x `nx`, row-major) */
        std::vector<realtype> sx0;

        /** status code */
        int status = 0;

        /** The solver used */
        std::string solver_type;

        /** maximal number of occurring events (for every event type) */
        unsigned int nmaxevent{0};

        /** number of considered timepoints */
        size_t nt{0};

        /** sensitivity order */
        SensitivityOrder sensi{SensitivityOrder::none};

        /** maximal number of newton iterations for steady state calculation */
        size_t newton_maxsteps{0};

        /** sensitivity method */
        SensitivityMethod sensi_meth{SensitivityMethod::none};

        LinearMultistepMethod lmm_solver {LinearMultistepMethod::BDF};

        LinearSolver ls_solver {LinearSolver::dense};

    protected:

        /** timepoint for model evaluation*/
        realtype t_{0};

        /** partial state vector */
        Vector x_solver_;

        /** partial time derivative of state vector */
        Vector dx_solver_;

        /** partial sensitivity state vector array */
        VectorArray sx_solver_;

        /** The state vector */
        Vector x_rdata_;

        /** The sensitivity state vector array */
        VectorArray sx_rdata_;

        /** array of number of found roots for a certain event type
        * (shape `ne`) */
        std::vector<unsigned int> nroots_;

        void initiainitializeReporting();

        /**
         * @brief extracts results from forward problem
         * @param fwd forward problem
         * @param model model that was used for forward simulation
         */
        void processForwardProblem(ForwardProblem const &fwd,
                                   Model &model);

        /**
         * @brief extracts results from solver
         * @param solver solver that was used for forward/backward simulation
         */
        void processSolver(Solver const &solver);

        /**
         * @brief sets member variables and model state according to provided simulation state
         * @param state simulation state provided by Problem
         * @param model model that was used for forward/backward simulation
         */
        void readSimulationState(SimulationState const &state, Model &model);

        /**
         * @brief Set state variables, outputs and respective
         * sensitivities to NaN (typically after integration failure)
         * @param it_start time index at which to start invalidating
         */
        void invalidate(size_t it_start);

        /**
         * @brief Checks whether forward sensitivity analysis is performed
         * @return boolean indicator
         */
        bool computingFSA() const {
            return (sensi_meth == SensitivityMethod::forward &&
                    sensi >= SensitivityOrder::first);
        }

        /**
         * @brief Extracts output information for data-points, expects that
         * x_solver_ and sx_solver_ were set appropriately
         * @param it timepoint index
         * @param model model that was used in forward solve
         */
        void getDataOutput(size_t it, Model &model);

        /**
         * @brief Extracts output information for events, expects that x_solver_
         * and sx_solver_ were set appropriately
         * @param t event timepoint
         * @param rootidx information about which roots fired
         * (1 indicating fired, 0/-1 for not)
         * @param model model that was used in forward solve
         */
        void getEventOutput(realtype t, std::vector<int> rootidx, Model &model);

        /**
         * @brief Extracts data information for forward sensitivity analysis,
         * expects that x_solver_ and sx_solver_ were set appropriately
         * @param it index of current timepoint
         * @param model model that was used in forward solve
         */
        void getDataSensisFSA(size_t it, Model &model);

    };
}

#endif //SUNEIGEN_RETURN_DATA_H
