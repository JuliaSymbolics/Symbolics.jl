using Symbolics

#lambert w function is the inverse of x*exp(x)
function lambertw(x)#only the upper branch
	if x isa Number
		g = log(x+1.0)*0.75#initial guess, this function is close to lambert w function
		
		#newtons method
		for i in 1:7
			g = g-(g*exp(g)-x)*exp(-g)/(1+g)
		end
		
		return g
	else
		return term(lambertw,x;type = Number)
	end
end

#reduce the system of equations by one equation
function remove_eq(equs::Vector{Equation},var,removed::Vector{Equation})
    for i in 1:length(equs)
        solution = solve_single_eq(equs[i],var,true)
        if(solution isa Vector)
            solution = solution[1]
        end
        solution == nothing && continue

        push!(removed,solution)
        deleteat!(equs,i)
        
        for j in 1:length(equs)
            equs[j] = substitute(equs[j],Dict(solution.lhs => solution.rhs))
        end
        return
    end
end

#returns a dictionary of the solutions
function solve_system_eq(equs::Vector{Equation},vars)
    removed = Vector{Equation}()

    reduced::Vector{Equation} = copy(equs)
    
    for i in 1:length(vars)
        remove_eq(reduced,vars[i],removed)
    end
    
    #solve last equation in the reduced set
    
    solutions = Dict()

    for i in length(removed):-1:1
        current_eq = substitute(removed[i],solutions)
        solutions[current_eq.lhs]=current_eq.rhs
    end

    return solutions
end

inverseOps = Dict(
                  sin => asin,
                  cos => acos,
                  tan => atan,
                  asin => sin,
                  acos => cos,
                  atan => tan,
                  exp => log,
                  log => exp)

function solve_single_eq(eq::Equation,var,single_solution = false)
    eq = (SymbolicUtils.add_with_div(eq.lhs+-1*eq.rhs) ~ 0)#move everything to the left side
    while(true)
        oldState = eq
        
        if(istree(eq.lhs))
            op = operation(eq.lhs)
            
            
            if(isequal(degree(eq.lhs,var),2) && op == +) #quadratic
                coeffs = polynomial_coeffs(eq.lhs,[var])
                a = coeffs[1][var^2]
                b = haskey(coeffs[1],var) ? coeffs[1][var] : 0
                c = coeffs[2]-eq.rhs

				if(single_solution)
					return var ~  simplify_fractions(SymbolicUtils.Div(-b+term(sqrt,b^2-4*a*c),2*a))
				else
                	return [var ~  simplify_fractions(SymbolicUtils.Div(-b+term(sqrt,b^2-4*a*c),2*a)) , var ~  simplify_fractions(SymbolicUtils.Div(-b-term(sqrt,b^2-4*a*c),2*a)) ]
                end
            
            elseif(op in [+,*])#N argumented types
                elements = arguments(eq.lhs)
                if (op == +) && length(elements) == 2 && sum(istree.(elements))==length(elements) && isequal(operation.(elements),[sqrt for el=1:length(elements)]) #check for sqrt(a)+sqrt(b)=c form , to solve this sqrt(a)+sqrt(b)=c -> 4*a*b-full_expand((c^2-b-a)^2)=0 then solve using quadratics
                	
                	#grab values
                	a = (elements[1]).arguments[1]
                	b = (elements[2]).arguments[1]
                	c = eq.rhs
                	
                	
                	eq = expand(2*b*a)-expand(a^2)-expand(b^2)-expand(c^4)+expand(2*a*c^2)+expand(2*b*c^2)~0
                	
                elseif (op == +) && length(elements) == 2 && istree(elements[2]) && operation(elements[2]) == sqrt #a+sqrt(b)=c -> b-a^2+2*a*c-c^2=0
                	a = elements[1]
                	b = (elements[2]).arguments[1]
                	c = eq.rhs
                	
                	eq = b-a^2+2*a*c-c^2~0
                elseif (op == +) && isequal(eq.rhs,0) && length(elements) == 2 && sum(istree.(elements))==length(elements) && length(arguments(elements[1])) == 2 && isequal(arguments(elements[1])[1],-1) && istree(arguments(elements[1])[2]) && operation(elements[2]) == operation(arguments(elements[1])[2])#-f(y)+f(x)=0 -> x-y=0
                	
                	x = arguments(elements[2])[1]
                	y = arguments(arguments(elements[1])[2])[1]
                	
                	eq = x-y~0
                elseif (op == *) && istree(elements[2]) && operation(elements[2]) == log && isequal(elements[1],arguments(elements[2])[1])#x*ln(x)=y -> x=e^lambertW(y)
                	x = elements[1]
                	y = eq.rhs
                	eq = x~term(exp,term(lambertw,y))
                else#standard procedure
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
                end
                
            elseif(op == /)#reverse division
                eq = eq.lhs.num-eq.lhs.den*eq.rhs ~ 0
            elseif(op == ^)#reverse powers
                pow = eq.lhs

                baseHasVar = SymbolicUtils._occursin(var,pow.base)
                expoHasVar = SymbolicUtils._occursin(var,pow.exp)

                if(baseHasVar && !expoHasVar)
                    twoSolutions = !single_solution && isequal(pow.exp%2,zero(pow.exp))
                    
                    if(twoSolutions) 

                        eq1 = solve_single_eq( pow.base ~ eq.rhs^(SymbolicUtils.Div(1,pow.exp)) , var)
                        eq2 = solve_single_eq( pow.base ~ -eq.rhs^(SymbolicUtils.Div(1,pow.exp)) , var)
                    
                        return [eq1,eq2]
                    else
                        eq = pow.base ~ eq.rhs^(SymbolicUtils.Div(1,pow.exp))
                    end
                elseif(!baseHasVar && expoHasVar)
                    eq = pow.exp ~ SymbolicUtils.Div( term(log,eq.rhs) , term(log,pow.base) )
                elseif(baseHasVar && expoHasVar)
                	if isequal(pow.exp,pow.base)
                		eq = pow.exp ~ SymbolicUtils.Div( term(log,eq.rhs) , term(lambertw,term(log,eq.rhs)) )
                	else
                		eq = pow.exp*term(log,pow.base) ~ term(log,eq.rhs)
                	end
                end
         	elseif(op == sqrt)
         		inner = arguments(eq.lhs)[1]
         		eq = inner~(eq.rhs)^2
         	elseif(op == lambertw)
         		inner = arguments(eq.lhs)[1]
         		eq = inner~eq.rhs*term(exp,eq.rhs)
            elseif haskey(inverseOps,op)
                inverseOp = inverseOps[op]
                inner = arguments(eq.lhs)[1]
                eq = inner~term(inverseOp,eq.rhs)
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
