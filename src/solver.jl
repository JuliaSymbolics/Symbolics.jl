using Symbolics

inverseOps = Dict(
                  sin => asin,
                  cos => acos,
                  tan => atan,
                  asin => sin,
                  acos => cos,
                  atan => tan,
                  exp => log,
                  log => exp)

function solve_single_eq(eq::Equation,var)
    
    eq = (eq.lhs-eq.rhs ~ 0)#move everything to the left side
    while(true)
        oldState = eq
        #println("eq $eq var $var")
        if(istree(eq.lhs))
            op = operation(eq.lhs)
            
            if(isequal(degree(eq.lhs,var),2) && op == +)
                coeffs = polynomial_coeffs(eq.lhs,[var])
                a = coeffs[1][var^2]
                b = haskey(coeffs[1],var) ? [var] : 0
                c = coeffs[2]-eq.rhs

                return [var ~ SymbolicUtils.Div(-b+term(sqrt,b^2-4*a*c),2*a) , var ~ SymbolicUtils.Div(-b-term(sqrt,b^2-4*a*c),2*a) ]

            elseif(op in [+,*])#N argumented types
                elements = arguments(eq.lhs)
                stays = []
                move = []
                for i in 1:length(elements)
                    hasVar = SymbolicUtils._occursin(var,elements[i])
                    if(hasVar)
                        push!(stays, elements[i])
                    else
                        push!(move, elements[i])
                    end
                end
            
                if(op == +)#reverse addition
                    eq = (length(stays) == 0 ? 0 : +(stays...)) ~ -( length(move) == 0 ? 0 : +(move...) )+eq.rhs
                elseif(op == *)#reverse multiplication
                    eq = (length(stays) == 0 ? 1 : *(stays...)) ~ SymbolicUtils.Div(eq.rhs , (length(move) == 0 ? 1 : *(move...) ) )
                end
            elseif(op == /)#reverse division
                eq = eq.lhs.num-eq.lhs.den*eq.rhs ~ 0
            elseif(op == ^)#reverse powers
                pow = eq.lhs

                baseHasVar = SymbolicUtils._occursin(var,pow.base)
                expoHasVar = SymbolicUtils._occursin(var,pow.exp)

                if(baseHasVar && !expoHasVar)
                    twoSolutions = isequal(pow.exp%2,zero(pow.exp))
                    
                    if(twoSolutions) 

                        eq1 = solve_single_eq( pow.base ~ eq.rhs^(SymbolicUtils.Div(1,pow.exp)) , var)
                        eq2 = solve_single_eq( pow.base ~ -eq.rhs^(SymbolicUtils.Div(1,pow.exp)) , var)
                    
                        return [eq1,eq2]
                    else
                        eq = pow.base ~ eq.rhs^(SymbolicUtils.Div(1,pow.exp))
                    end
                elseif(!baseHasVar && expoHasVar)
                    eq = pow.exp ~ SymbolicUtils.Div( term(log,eq.rhs) , term(log,pow.base) )
                end
            elseif haskey(inverseOps,op)
                inverseOp = inverseOps[op]
                inner = arguments(eq.lhs)[1]
                eq = inner~term(inverseOp,eq.rhs)
            else#pattern matched equation rules go here
                
            end

        end
        if(isequal(eq.lhs, var))
            return eq
        end


        if(isequal(eq,oldState))
            return nothing#unsolvable with these methods
        end
    end

end
