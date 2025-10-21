module PiezoresistiveDesign

#bababab

using Optim, LineSearches
using Mabla.FEINNs

using DrWatson, DataFrames, JLD2, CSV

using BlockArrays, LinearAlgebra
using Gridap, Gridap.MultiField, Gridap.Algebra, Gridap.Adaptivity, GridapTopOpt# , GridapSolvers

using ForwardDiff, ReverseDiff

using GridapEmbedded, GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: DifferentiableTriangulation

using Gridap.Arrays
using Gridap.Fields
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.TensorValues
using GridapGmsh

using Zygote 

include(srcdir("ADext.jl"))   

include(srcdir("forwardproblem/measures.jl"))
include(srcdir("forwardproblem/spaces.jl"))
include(srcdir("forwardproblem/functionals.jl"))
include(srcdir("forwardproblem/IPforwardpass.jl"))
include(srcdir("forwardproblem/TOforwardpass.jl"))
include(srcdir("forwardproblem/electrodepositions.jl"))
include(srcdir("forwardproblem/ramping.jl"))
#include(srcdir("forwardproblem/embeddedcollection.jl"))
# include(srcdir("optimisation/reconstruction.jl"))
include(srcdir("sweepconfig/sweepconfig.jl"))
include(srcdir("sweepconfig/saverun.jl"))
include(srcdir("sweepconfig/runcasefunc.jl"))
include(srcdir("driver.jl"))
include(srcdir("optimisation/shapes.jl"))
include(srcdir("optimisation/helpers.jl"))
include(srcdir("optimisation/optimisation.jl"))


#include(srcdir("chainrules.jl"))

export run_problem
export pbs_sweeper
export save_job_dicts
export run_case_function

export PDEConstrainedLoss2
export CustomPDEConstrainedFunctionals

import GridapTopOpt: StateParamIntegrandWithMeasure
import GridapTopOpt: StateParamMap
#import GridapTopOpt: IntegrandWithMeasure

println("done loading")


end