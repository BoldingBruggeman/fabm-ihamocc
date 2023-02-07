add_library(fabm_models_ihammoc OBJECT
            ihammoc_model_library.F90
            # List additional source code file
           )

add_dependencies(fabm_models_ihammoc fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2023 - BB
#-----------------------------------------------------------------------
