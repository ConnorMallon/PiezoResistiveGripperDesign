module postprocess

using DataFrames
using DrWatson; 

@quickactivate 
#include(srcdir("PreLoad.jl"))

import DrWatson.collect_results

using CSV, DataFrames

cols_to_filter = ["ts","Δt","w1","n","cost","pd","Vₘₐₓ","Rₚₑₗ","ϵ_macro","Rf","iterations","effective_vol","cbf_its","c_out","m_stored"]

model = "test" #"123_more_its"

results15  = DrWatson.collect_results(datadir(joinpath("results",model))) #,white_list=cols_to_filter)
@show names(results15)
results15.J1[1]

using Gridap
using Gridap.Arrays
import Gridap.Geometry: get_node_coordinates
using GridapEmbedded
using GridapEmbedded.LevelSetCutters

ϕf = ϕf.-0.55

geo1 = DiscreteGeometry(ϕf,point_to_coords,name="")
geo2 = DiscreteGeometry(-ϕf,point_to_coords,name="")

cutgeo1= cut(bgmodel, geo1)
cutgeo2= cut(bgmodel, geo2)

# Setup interpolation meshes
Ω1_act = Triangulation(cutgeo1)#,cutgeo1)
Ω2_act = Triangulation(cutgeo2)#,cutgeo2)

Γ1_act = EmbeddedBoundary(cutgeo1)
writevtk(Γ1_act,"Γ1_act_geo1")

writevtk(Ω2_act,"Ω1_act_geo2")

end

