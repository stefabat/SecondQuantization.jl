
import Base.*
import Base.+
import Base.print
import Base.eval


"""Abstract type representing a second-quantized object"""
abstract type SQObject end


abstract type ElementaryOperator <: SQObject end


"""Struct representing a multiplication object"""
struct Mul <: SQObject
  args::Vector{SQObject}
  coef::Float64
end

function latex(mul::Mul)
  coef = mul.coef
  out = "$coef"
  for arg in mul.args
    out *= latex(arg)
  end
  return out
end

"""Default constructor with coefficient 1.0 for the Mul struct"""
function Mul(args::Vector{SQObject})
  return Mul(args, 1.0)
end

"""Default constructor with coefficient 1.0 for the Mul struct"""
function Mul(args::Vector{ElementaryOperator})
  return Mul(args, 1.0)
end


"""Struct representing an addition object"""
struct Add <: SQObject
  args::Vector{SQObject}
end


function latex(add::Add)
  out = ""
  for arg in add.args
    out *= latex(arg) * " + "
  end
  return out
end


"""Struct representing the zero value"""
struct Zero <: SQObject end

function latex(zero::Zero)
  return "0"
end

function *(obj1::SQObject, obj2::SQObject)
  return Mul([obj1,obj2])
end

function *(obj::SQObject, zero::Zero)
  return Zero()
end

function *(zero::Zero, obj::SQObject)
  return Zero()
end

function *(obj::SQObject, mul::Mul)
  return Mul(vcat(obj,mul.args), mul.coef)
end

function *(mul::Mul, obj::SQObject)
  return Mul(vcat(mul.args,obj), mul.coef)
end

function *(num::Number, mul::Mul)
  return Mul(mul.args, num*mul.coef)
end

function *(mul::Mul, num::Number)
  return Mul(mul.args, num*mul.coef)
end

function *(num::Number, obj::SQObject)
  return Mul([obj], num)
end

function *(obj::SQObject, num::Number)
  return Mul([obj], num)
end

function *(mul1::Mul, mul2::Mul)
  return Mul(vcat(mul1.args,mul2.args), mul1.coef*mul2.coef)
end

function +(obj1::SQObject, obj2::SQObject)
  return Add([obj1,obj2])
end

function +(obj::SQObject, zero::Zero)
  return obj
end

function +(zero::Zero, obj::SQObject)
  return obj
end

function +(obj::SQObject, add::Add)
  return Add(vcat(obj,add.args))
end

function +(add::Add, obj::SQObject)
  return Add(vcat(add.args,obj))
end

function +(add1::Add, add2::Add)
  return Add(vcat(add1.args,add2.args))
end


"""Struct representing an orbital index"""
struct OrbitalIndex
  index::Symbol
  class::String
end

"""Returns `true` if the index represents an inactive orbital"""
function is_inactive(index::OrbitalIndex)
  return index.class == "inactive"
end

"""Returns `true` if the index represents an active orbital"""
function is_active(index::OrbitalIndex)
  return index.class == "active"
end

"""Returns `true` if the index represents a virtual orbital"""
function is_virtual(index::OrbitalIndex)
  return index.class == "virtual"
end


"""Struct representing a fermionic creation operator"""
struct CreationOperator <: ElementaryOperator
  index::OrbitalIndex
end


function latex(op::CreationOperator)
  index = string(op.index.index)
  return "a^{\\dagger}_{$index}"
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


"""Struct representing a fermionic annihilation operator"""
struct AnnihilationOperator <: ElementaryOperator
  index::OrbitalIndex
end


function latex(op::AnnihilationOperator)
  index = string(op.index.index)
  return "a_{$index}"
end

function is_general(op::ElementaryOperator)
  return op.index.class == "general"
end

function is_inactive(op::ElementaryOperator)
  return op.index.class == "inactive"
end

function is_active(op::ElementaryOperator)
  return op.index.class == "active"
end

function is_virtual(op::ElementaryOperator)
  return op.index.class == "virtual"
end

function is_hole(op::ElementaryOperator)
  return is_inactive(op) || is_active(op)
end

function is_particle(op::ElementaryOperator)
  return is_active(op) || is_virtual(op)
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



"""Struct representing a normal-ordered string of operators"""
struct NO <: SQObject
  args::Vector{SQObject}
end

#function NO(args::Vector{SQObject})
#  return NO(args)
#end

function latex(no::NO)
  out = ""
  for arg in NO.args
    out *= latex(arg)
  end
  return "\\lbrace"*out*"\\rbrace"
end


"""Struct representing a Kronecker delta function"""
struct KroneckerDelta <: SQObject
  upper::OrbitalIndex
  lower::OrbitalIndex
end


function latex(δ::KroneckerDelta)
  up = string(δ.upper.index)
  dw = string(δ.lower.index)
  return "\\delta_{$up$dw}"
end



"""Constructor of Kronecker delta function"""
function KroneckerDelta(op1::ElementaryOperator, op2::ElementaryOperator)
  return KroneckerDelta(op1.index,op2.index)
end

"""Explicitly evaluate a Kronecker delta function"""
function eval(δ::KroneckerDelta)
  if δ.upper.class == δ.lower.class
    return δ
  else
    return Zero()
  end
end


"""Struct representing an antisymmetric tensor"""
struct AntiSymmetricTensor <: SQObject
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



"""Struct representing a commutator"""
struct Commutator <: SQObject
  A::SQObject
  B::SQObject
end

"""Explicitly evaluate a commutator"""
function eval(comm::Commutator)
  return comm.args[1]*comm.args[2] + (-1.0)*comm.args[2]*comm.args[1]
end

"""Struct representing a anticommutator"""
struct AntiCommutator <: SQObject
  A::SQObject
  B::SQObject
end

function AntiCommutator(A::CreationOperator, B::CreationOperator)
  return Zero()
end

function AntiCommutator(A::AnnihilationOperator, B::AnnihilationOperator)
  return Zero()
end

function AntiCommutator(A::CreationOperator, B::AnnihilationOperator)
  return KroneckerDelta(A,B)
end


"""Struct representing a n-particle reduced density matrix"""
struct Gamma <: SQObject
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
struct Eta <: SQObject
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


"""Applies to quasi-operators"""
function contraction(op1::CreationOperator, op2::CreationOperator)
  return Zero()
end

"""Applies to quasi-operators"""
function contraction(op1::AnnihilationOperator, op2::AnnihilationOperator)
  return Zero()
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
  return Zero()
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
  return Zero()
end


"""Apply Wicks' theorem"""
function wicks(op::ElementaryOperator)
  return op
end


"""Apply Wicks' theorem"""
function wicks(normal_ordered_string::NO)
  return normal_ordered_string
end


"""Apply Wicks' theorem"""
function wicks(add::Add)
  return Add([wicks(arg) for arg in add.args])
end


"""
In theory, the mul object enetering here should only be composed by elementary
operators.
"""
function wicks(mul::Mul)
  # vector containing all terms
  result = SQObject[]
  push!(result,NO(mul))

  # loop over all operators
  for i=1:length(mul)-1
    for j=i+1:length(mul)

      # compute contraction between the two Fermionic operators
      contracted_pair = contraction(opstring[i],opstring[j])

      # check that the contraction is non-zero
      #if typeof(contracted_pair) != Zero
        sign  = (-1.0)^((j - i + 1) % 2)
        coeff = sign*contracted_pair
      #end

      # remaining operators to get contractions from
      # we exclude operators opstring[1:i-1] to avoid
      # double counting
      oplist = vcat(opstring[i+1:j-1], opstring[j+1:end])

      if length(oplist) > 0
        push!(result,coeff*NO(opstring[1:i-1])*get_contractions(oplist))
      else
        if length(opstring[1:i-1]) > 0
          push!(result,coeff*NO(opstring[1:i-1]))
        else
          push!(result,coeff)
        end
      end

    end
  end

  return Add(result)
end


"""Apply Wicks' theorem"""
function wicks(mul::Mul)
end





