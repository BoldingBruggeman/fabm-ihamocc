add_library(fabm_models_ihamocc OBJECT
            ihamocc_model_library.F90
			shared.F90
			oxygen.F90
            carbon.F90
            bromo.F90
            cfc.F90
            dms.F90
            alkalinization.F90
            iron.F90
            ndep.F90
            preftrc.F90
            phytoplankton.F90
            zooplankton.F90
            detritus.F90
            # List additional source code file
           )

add_dependencies(fabm_models_ihamocc fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2023 - BB
#-----------------------------------------------------------------------
