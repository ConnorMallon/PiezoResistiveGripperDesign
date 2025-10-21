using Gridap, Gridap.FESpaces, Gridap.MultiField, Gridap.CellData, Gridap.Helpers, Gridap.Fields
using Gridap.FESpaces: _change_argument, _compute_cell_ids, autodiff_array_gradient
using Gridap.MultiField: blocks, MultiFieldFEFunction, mortar
using Gridap.CellData: is_change_possible
using Test

function concat_contribs(size,contribs...)
  ArrayBlock([contribs...],size)
end

GetIndex(k) = i->getindex(i,k)

function _check_trians(f::MultiFieldFEFunction)
  trians = map(get_triangulation,f.fe_space.spaces)
  trian = first(trians)
  all(t -> t===trian, trians)
end

function FESpaces._gradient(f,uh::MultiFieldFEFunction,fuh::DomainContribution)
  if _check_trians(uh)
    _mf_gradient_same_trian(f,uh,fuh)
  else
    _mf_gradient(f,uh,fuh)
  end
end

function _mf_gradient_same_trian(f,uh,fuh)
  terms = DomainContribution()
  for trian in get_domains(fuh)
    g = _change_argument(Gridap.gradient,f,trian,uh)
    cell_u = get_cell_dof_values(uh)
    cell_id = _compute_cell_ids(uh,trian)
    cell_grad = autodiff_array_gradient(g,cell_u,cell_id)
    add_contribution!(terms,trian,cell_grad)
  end
  terms
end

function _mf_gradient(f,uh,fuh)
  _uh = uh.single_fe_functions
  terms = [DomainContribution() for _ in 1:num_fields(uh)]
  for k in 1:num_fields(uh)
    cell_u = get_cell_dof_values(uh[k])
    for trian in get_domains(fuh)
      g = _change_argument(Gridap.gradient,uk->f((_uh[1:k-1]...,uk,_uh[k+1:end]...)),trian,uh[k])
      cell_id = _compute_cell_ids(uh[k],trian)
      cell_grad = autodiff_array_gradient(g,cell_u,cell_id)
      add_contribution!(terms[k],trian,cell_grad)
    end
  end

  contribs = DomainContribution()
  for trian in get_domains(fuh)
    trian_to_contrib = lazy_map(GetIndex(trian),terms)
    contrib_to_touched = fill([true for _ in 1:num_fields(uh)],length(first(trian_to_contrib)));
    mf_cell_grad = lazy_map(concat_contribs,contrib_to_touched,trian_to_contrib...);
    add_contribution!(contribs,trian,mf_cell_grad)
  end
  contribs
end