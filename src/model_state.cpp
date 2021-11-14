#include "model_state.h"

namespace suneigen {

    ModelStateDerived::ModelStateDerived(const ModelDimensions &dim)
            : J_(static_cast<sunindextype>(dim.nx), static_cast<sunindextype>(dim.nx),
                 static_cast<sunindextype>(dim.nnz), CSC_MAT),
              x_rdata_(dim.nx, 0.0),
              sx_rdata_(dim.nx, 0.0),
              x_pos_tmp_(dim.nx) {}

}
