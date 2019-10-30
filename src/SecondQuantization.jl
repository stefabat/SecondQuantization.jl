
import Base.*
import Base.+
import Base.print
import Base.eval

#### Orbital index ####

"""Struct representing an orbital index"""
struct OrbitalIndex
  index::Symbol
  class::String
end

get_index(orbind::OrbitalIndex) = orbind.index
get_class(orbind::OrbitalIndex) = orbind.class

"""Returns `true` if the index represents an inactive orbital"""
function is_inactive(orbind::OrbitalIndex)
  return get_class(orbind) == "inactive"
end

"""Returns `true` if the index represents an active orbital"""
function is_active(orbind::OrbitalIndex)
  return get_class(orbind) == "active"
end

"""Returns `true` if the index represents a virtual orbital"""
function is_virtual(orbind::OrbitalIndex)
  return get_class(orbind) == "virtual"
end


#### Creation operator ####

"""Struct representing a fermionic creation operator"""
struct CreationOperator
  index::OrbitalIndex
end

"""Constructor of a creation operator"""
function CreationOperator(p::Symbol)
  return CreationOperator(OrbitalIndex(p,"general"))
end

"""Base constructor of a creation operator"""
function CreationOperator(p::Symbol, class::String)
  if class in ["general","inactive","active","virtual"]
    return CreationOperator(OrbitalIndex(p,class))
  else
    throw(ErrorException("operator class not known"))
  end
end

"""Shorthand constructor of a creation operator"""
function Fd(p::Symbol)
  return CreationOperator(OrbitalIndex(p,"general"))
end

"""Shorthand base constructor of a creation operator"""
function Fd(p::Symbol, class::String)
  return CreationOperator(OrbitalIndex(p,class))
end

"""Returns the latex code of a creation operator"""
function latex(op::CreationOperator)
  index = string(get_index(op))
  return "a^{\\dagger}_{$index}"
end


#### Annihilation operator ####

"""Struct representing a fermionic annihilation operator"""
struct AnnihilationOperator
  index::OrbitalIndex
end

"""Constructor of an annihilation operator"""
function AnnihilationOperator(p::Symbol)
  return AnnihilationOperator(OrbitalIndex(p,"general"))
end

"""Base constructor of an annihilation operator"""
function AnnihilationOperator(p::Symbol, class::String)
  if class in ["general","inactive","active","virtual"]
    return AnnihilationOperator(OrbitalIndex(p,class))
  else
    throw(ErrorException("operator class not known"))
  end
end

"""Shorthand constructor of an annihilation operator"""
function F(p::Symbol)
  return AnnihilationOperator(OrbitalIndex(p,"general"))
end

"""Shortrhand base constructor of an annihilation operator"""
function F(p::Symbol, class::String)
  return AnnihilationOperator(OrbitalIndex(p,class))
end

function latex(op::AnnihilationOperator)
  index = string(get_index(op))
  return "a_{$index}"
end

"""Union type describing a Fermionic operator"""
FermionicOperator = Union{CreationOperator,AnnihilationOperator}

get_index(op::FermionicOperator) = get_index(op.index)
get_class(op::FermionicOperator) = get_class(op.index)


#### functions on general elementary operators ####

function is_general(op::FermionicOperator)
  return get_class(op) == "general"
end

function is_inactive(op::FermionicOperator)
  return get_class(op) == "inactive"
end

function is_active(op::FermionicOperator)
  return get_class(op) == "active"
end

function is_virtual(op::FermionicOperator)
  return get_class(op) == "virtual"
end

function is_hole(op::FermionicOperator)
  return is_inactive(op) || is_active(op)
end

function is_particle(op::FermionicOperator)
  return is_active(op) || is_virtual(op)
end


"""Struct representing a product of Fermionic operators"""
mutable struct Product
  coefficient::Float64
  operators::Vector{FermionicOperator}
end

function Product(ops::Vector{FermionicOperator})
  return Product(1.0, ops)
end

get_coefficient(product::Product) = product.coefficient
get_operators(product::Product) = product.operators
length(product::Product) = length(product.operators)

function set_coefficient(product::Product, coefficient::Number)
  product.coefficient = coefficient
end

function latex(product::Product)
  coef = get_coefficient(product)
  out = string(coef)
  for op in get_operators(product)
    out *= latex(op)
  end
  return out
end

#### Algebraic rules for Fermionic operators ####

function *(op1::FermionicOperator, op2::FermionicOperator)
  return Product([op1,op2])
end

function *(product::Product, op::FermionicOperator)
  coef = get_coefficient(product)
  ops  = get_operators(product)
  return Product(coef,vcat(ops,op))
end

function *(op::FermionicOperator, product::Product)
  coef = get_coefficient(product)
  ops  = get_operators(product)
  return Product(coef,vcat(op,ops))
end

function *(product1::Product, product2::Product)
  coef = get_coefficient(product1) * get_coefficient(product2)
  ops  = vcat(get_operators(product1),get_operators(product2))
  return Product(coef,ops)
end

function *(x::Number, op::FermionicOperator)
  return Product(x, [op])
end

function *(op::FermionicOperator, x::Number)
  return Product(x, [op])
end

function *(product::Product, x::Number)
  coef = x * get_coefficient(product)
  ops  = get_operators(product)
  return Product(coef,ops)
end

function *(x::Number, product::Product)
  coef = x * get_coefficient(product)
  ops  = get_operators(product)
  return Product(coef,ops)
end


struct NO
  product::Product
end

# function latex(no::NO)
#   out = ""
#   for arg in NO.args
#     out *= latex(arg)
#   end
#   return "\\lbrace"*out*"\\rbrace"
# end


#### Kronecker delta function ####

"""Struct representing a Kronecker delta function"""
struct KroneckerDelta
  upper::OrbitalIndex
  lower::OrbitalIndex
end

function latex(δ::KroneckerDelta)
  up = string(δ.upper.index)
  dw = string(δ.lower.index)
  return "\\delta_{$up$dw}"
end

"""Constructor of Kronecker delta function"""
function KroneckerDelta(op1::FermionicOperator, op2::FermionicOperator)
  return KroneckerDelta(op1.index,op2.index)
end

"""Explicitly evaluate a Kronecker delta function"""
function eval(δ::KroneckerDelta)
  if δ.upper.class == δ.lower.class
    return 1.0
  else
    return 0.0
  end
end


"""Struct representing an antisymmetric tensor"""
struct AntiSymmetricTensor
  symbol::Symbol
  upper::Vector{OrbitalIndex}
  lower::Vector{OrbitalIndex}
end

function latex(ast::AntiSymmetricTensor)
  sym = string(ast.symbol)
  ups = ""
  dws = ""
  for up in ast.upper
    ups *= string(up.index.index)
  end
  for dw in ast.lower
    dws *= string(dw.index.index)
  end
  return "$sym^{$ups}_{$dws}"
end


# """Struct representing a commutator"""
# struct Commutator <: SQObject
#   A::SQObject
#   B::SQObject
# end

# """Explicitly evaluate a commutator"""
# function eval(comm::Commutator)
#   return comm.args[1]*comm.args[2] + (-1.0)*comm.args[2]*comm.args[1]
# end

# """Struct representing a anticommutator"""
# struct AntiCommutator <: SQObject
#   A::SQObject
#   B::SQObject
# end

# function AntiCommutator(A::CreationOperator, B::CreationOperator)
#   return Zero()
# end

# function AntiCommutator(A::AnnihilationOperator, B::AnnihilationOperator)
#   return Zero()
# end

# function AntiCommutator(A::CreationOperator, B::AnnihilationOperator)
#   return KroneckerDelta(A,B)
# end


"""Struct representing a n-particle reduced density matrix"""
struct Gamma
  upper::Vector{OrbitalIndex}
  lower::Vector{OrbitalIndex}
end

function latex(γ::Gamma)
  ups = ""
  dws = ""
  for up in γ.upper
    ups *= string(up.index.index)
  end
  for dw in γ.lower
    dws *= string(dw.index.index)
  end
  return "\\gamma^{$ups}_{$dws}"
end

"""Struct representing a n-hole reduced density matrix"""
struct Eta
  upper::Vector{OrbitalIndex}
  lower::Vector{OrbitalIndex}
end

function latex(η::Eta)
  ups = ""
  dws = ""
  for up in η.upper
    ups *= string(up.index.index)
  end
  for dw in η.lower
    dws *= string(dw.index.index)
  end
  return "\\eta^{$ups}_{$dws}"
end



"""Struct representing a contraction of Fermionic operators"""
mutable struct Contraction
  operators::Vector{FermionicOperator}
end


"""Applies to quasi-operators"""
function contraction(op1::CreationOperator, op2::CreationOperator)
  return 0
end

"""Applies to quasi-operators"""
function contraction(op1::AnnihilationOperator, op2::AnnihilationOperator)
  return 0
end

"""Applies to quasi-operators"""
function contraction(op1::CreationOperator, op2::AnnihilationOperator)
  if is_inactive(op2) || is_general(op2)
    if is_inactive(op1) || is_general(op1)
      return KroneckerDelta(op1,op2)
    end
  elseif is_active(op2) && is_active(op1)
    # 1-particle RDM
    return Gamma([op1.index],[op2.index])
  end
  # all other cases are zero
  return 0
end


"""Applies to quasi-operators"""
function contraction(op1::AnnihilationOperator, op2::CreationOperator)
  if is_virtual(op2) || is_general(op2)
    if is_virtual(op1) || is_general(op1)
      return KroneckerDelta(op1,op2)
    end
  elseif is_active(op2) && is_active(op1)
    # 1-hole RDM
    return Eta([op1.index],[op2.index])
  end
  # all other cases are zero
  return 0
end


"""Apply Wicks' theorem"""
function wicks(op::FermionicOperator)
  return op
end


"""Apply Wicks' theorem"""
function wicks(normal_ordered_string)
  return normal_ordered_string
end


"""Apply Wicks' theorem"""
function wicks(product::Product)
  # vector containing all terms
  result = Product[]
  push!(result,NO(product))

  # loop over all operators
  for i=1:length(product.operators)-1
    for j=i+1:length(product.operators)

      # compute contraction between the two Fermionic operators
      contracted_pair = contraction(product.operators[i],product.operators[j])

      sign  = (-1.0)^((j - i + 1) % 2)
      coeff = sign*contracted_pair

      # remaining operators to get contractions from
      # we exclude operators product[1:i-1] to avoid
      # double counting
      oplist = vcat(product.operators[i+1:j-1], product.operators[j+1:end])

      if length(oplist) > 0
        push!(result,coeff*NO(product.operators[1:i-1])*wicks(oplist))
      else
        if length(product.operators[1:i-1]) > 0
          push!(result,coeff*NO(product.operators[1:i-1]))
        else
          push!(result,coeff)
        end
      end

    end
  end

  return result
end
