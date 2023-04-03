add_library(fabm_models_ihamocc OBJECT
            ihamocc_model_library.F90
			shared.F90
			oxygen.F90
            # List additional source code file
           )

add_dependencies(fabm_models_ihamocc fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2023 - BB
#-----------------------------------------------------------------------
