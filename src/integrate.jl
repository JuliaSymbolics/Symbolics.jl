using Symbolics

#converts polynomial into array of coefficients
#ex  poly_extract(2+x^2+5*x^3) -> [2,0,1,5]
#its more useful in certain cituations than polynomial_coeffs

function poly_extract(expr,v)
	!SymbolicUtils._occursin(v,expr) && return [expr]#constant
	
	if isequal(expr,v)
		return [0,1]
	elseif(expr isa SymbolicUtils.Pow && expr.exp isa Int && expr.exp > 0 && isequal(expr.base,v))#power
		d = expr.exp
		out = fill!(Vector{Any}(undef,d+1),0)
		out[d+1] = 1
		return out
	elseif(expr isa SymbolicUtils.Mul)
		d = degree(expr,v)
		parts = copy(arguments(expr))
		if(d>0)
			for i in 1:length(parts)
				if SymbolicUtils._occursin(v,parts[i])
					deleteat!(parts,i)
					break
				end
			end
			
			out = fill!(Vector{Any}(undef,d+1),0)
			out[d+1] = *(parts...)
			
			return out
		end
	elseif(expr isa SymbolicUtils.Add)
		sum_parts = arguments(expr)
		d = degree(expr,v)
		
		out = fill!(Vector{Any}(undef,d+1),0)
		
		for sum_part in sum_parts
			if(!SymbolicUtils._occursin(v,sum_part))
				out[1]+=sum_part
			elseif(isequal(sum_part,v))
				out[2]+=1
			elseif(sum_part isa SymbolicUtils.Pow && sum_part.exp isa Int && sum_part.exp > 0 && isequal(sum_part.base,v))
				d2 = d = sum_part.exp
				out[d2+1]+=1
			elseif(sum_part isa SymbolicUtils.Mul)
				d2 = degree(sum_part,v)
				parts = copy(arguments(sum_part))
				if(d2>0)
					for i in 1:length(parts)
						if SymbolicUtils._occursin(v,parts[i])
							deleteat!(parts,i)
							break
						end
					end
					out[d2+1] += *(parts...)
				end
			end
			
			
		end
		
		return out
	end
end

#does the reverse of poly_extract ex to_poly([1,0,2],x) -> 1+2*x^2
function to_poly(coeffs,v)
	length(coeffs) == 0 && return 0
	out = coeffs[1]
	
	for i in 2:length(coeffs)
		d = i-1
		out+=coeffs[i]*v^d
	end
	
	return out
end


function poly_div(expr,v)
	if expr isa SymbolicUtils.Div
		numPoly = poly_extract(expr.num,v)
		denPoly = poly_extract(expr.den,v)
		
		if numPoly != nothing && denPoly != nothing && length(numPoly) >= length(denPoly)
			outPart , remainPart = poly_div_c(numPoly,denPoly)
			
			outPart = to_poly(outPart,v)
			remainPart = to_poly(remainPart,v)
			
			expr = outPart+SymbolicUtils.Div(remainPart,expr.den)
		end
	end
	return expr
end


function poly_div_c(num,den) #returns output + remainder
	remain = copy(num)
	out = []
	
	while length(remain) >= length(den)
		if isequal(den[length(den)],0) #avoid divide by zero situation
			pushfirst!(out,0)
			pop!(remain)#pop last element
			continue
		end
			
			
		coef = SymbolicUtils.Div(remain[length(remain)],den[length(den)])
		pushfirst!(out,coef)
		
		coef = -coef
			
			
		for i in (length(remain)-length(den)+1):length(remain)-1#we can skip last one since we know it will be deleted
			remain[i] += den[i-(length(remain)-length(den))] *coef
		end
		
		pop!(remain)#pop last element
			
			
	end
	while length(remain)>0 && isequal(remain[length(remain)],0) #clean zeros off end
		pop!(remain)#pop last element
	end
	
	return (out,remain);
end


function partial_fraction(expr, var)
    expr = unwrap(expr)
    var = unwrap(var)

    expr isa SymbolicUtils.Div && expr.den isa SymbolicUtils.Mul || return expr
    den, num = expr.den, expr.num

    den_total_deg = 0
    den_args = unsorted_arguments(den)
    sol = []

    if operation(den) == (*)
        for factor in den_args
            if factor isa SymbolicUtils.Pow
                b, a, islinear = linear_expansion(factor.base, var)
                islinear || return expr
                isinteger(factor.exp) || return expr
                den_total_deg += factor.exp
                push!(sol, (b, -a // b))
            elseif factor isa SymbolicUtils.Add
                b, a, islinear = linear_expansion(factor, var)
                if islinear
                    push!(sol, -a // b)
                    den_total_deg += 1
                else
                    return expr
                end
            end
        end
        den_total_deg <= 1 && return expr
        degree(num) < den_total_deg || return expr

        out = 0
        for (s, factor) in zip(sol, den_args)
            # f = (ax + b) * g
            g = simplify_fractions(factor * expr)
            if factor isa Add
                out += substitute(g, var=>s) / factor
            else # factor isa Pow
                b, s = s
                count = 0
                factorial = 1
                current_exp = factor.exp
                while (current_exp > 0)
                    numer = substitute(g, var=>s)
                    out += numer / (factor.base^current_exp * factorial * b^count)
                    current_exp -= 1
                    count += 1
                    factorial *= count
                    g = simplify_fractions(derivative(g, var))
                end
            end
        end
        return out
    else
        return expr
    end
end

#integration by parts scoring, LIATE method
BAD = 0
GOOD = 1
GREAT = 2

function rate_term(expr,x)
	isone(degree(expr,x)) && return GOOD
	
	expr isa SymbolicUtils.Pow && expr.exp isa Int && expr.exp > 0 && isone(degree(expr.base,x)) && return GOOD
	
	istree(expr) && operation(expr) == log && return GREAT
	return BAD
end

function rewrite_inv_to_frac(expr)#converts negative exponents to fractions
	if(expr isa SymbolicUtils.Pow && expr.exp isa Int && expr.exp < 0)# convert x^-1 to 1/x
		expr = 1/(expr.base)^(-expr.exp)
	elseif(expr isa SymbolicUtils.Mul)
		parts = arguments(expr)
		numer = []
		denom = []
		for i in 1:length(parts)
			current = parts[i]
			if(current isa SymbolicUtils.Pow && current.exp isa Int && current.exp < 0)
				push!(denom,(current.base)^(-current.exp))
			else
				push!(numer,current)
			end
		end
		
		if(length(denom)>0)
			expr = (*(numer...))/(*(denom...))
		end
	elseif istree(expr)
		parts = convert(Array{Any},arguments(expr))
		for i in 1:length(parts)
			parts[i] = rewrite_inv_to_frac(parts[i])
		end
		
		return (operation(expr))(parts...)
	end
	return expr
end

#measure max depth of tree
function complexity(expr)
	if(istree(expr))
		m = 0#max depth
		for el in arguments(expr)
			m = max(m,complexity(el))
		end
		return m+1	
	else
		return 1
	end
end

#get the next inner function of an expression
function inner_func(expr,x)
	if(istree(expr))
		parts = arguments(expr)
		best = nothing
		m = 0#max complexity
		for part in parts
			if (SymbolicUtils._occursin(x,part))
				c = complexity(part)
				if(c>m)
					m = c
					best = part
				end
			end
		end
		
		return best
	else
		return x	
	end
	
end

function integration_by_parts(expr,x)
	if expr isa SymbolicUtils.Mul#integration by parts
    	parts = arguments(expr)
    	derivativeChoice = nothing
    	highest = BAD
    	index_choice = 0
    	
    	for i in 1:length(parts)
    		part = parts[i]
    		score = rate_term(part,x)
    		if(score>highest)
    			derivativeChoice = part
    			index_choice = i
    			highest = score
    		end
    	end
    	
    	if derivativeChoice != nothing
    	
    		deleteat!(parts,index_choice)
    		
			integrateChoice = integrate(*(parts...),x)
			
			if integrateChoice != nothing && derivativeChoice != nothing
				recursiveIntegral = integrate(Symbolics.unwrap(Symbolics.derivative(derivativeChoice,x)) *integrateChoice,x)
				
				if 	recursiveIntegral != nothing
					return derivativeChoice*integrateChoice-recursiveIntegral
				end
			
			end
    	
    	end
    	
    end
    return nothing
end

function fix_form(expr)#convert *(element) -> element  makes sure its not a sum or product with only 1 element
	if(istree(expr))
		parts = arguments(expr)
		if(operation(expr) in [+,*] && length(parts) == 1)
			return parts[1]
		end
		for i in 1:length(parts)
			parts[i] = fix_form(parts[i])
		end
	end
	return expr
end

function u_sub(expr,x)
	@syms u_sub_var
	
	isequal(u_sub_var,x) && return nothing
	
	u = fix_form(expr)
	
	while !isequal(x,u)
		sub = u => u_sub_var
		#println(u)
		deriv = Symbolics.derivative(u,x)
		temp_expr = substitute(SymbolicUtils.Div(substitute(expr,sub), Symbolics.unwrap(deriv)),sub)#for some reason derivative returns a wrapped version of the expression
		
		if(!SymbolicUtils._occursin(x,temp_expr))
			attempt = integrate(temp_expr,u_sub_var)
			attempt != nothing && return substitute(attempt,u_sub_var=>u)
		else#try solving for the variable
			sol = solve_single_eq(u~u_sub_var,x,true)
			if(sol != nothing)
				sub2 = x => sol.rhs
				attempt = substitute(temp_expr,sub2)
				println("temp_expr = $temp_expr    attemp = $attempt")
			end
		end
		
		u = inner_func(u,x)
	end
	
	return nothing
end

function algebraicIntegration(expr,x)

	#rubi 1.1.1.1
	
	powerRule = @rule (~a)^(~n) => !SymbolicUtils._occursin(x,~n) && isone(degree(~a,x)) ? (~a)^(~n+1)/((~n+1)*polynomial_coeffs(~a,[x])[1][x] ) : nothing
	
	invRule = @rule (~k)/(~a) => !SymbolicUtils._occursin(x,~k) && isone(degree(~a,x)) ? ((~k)*log(~a))/(polynomial_coeffs(~a,[x])[1][x] ) : nothing
    	
    invPowerRule = @rule (~k)/(~a)^(~n) => !SymbolicUtils._occursin(x,~k) && !SymbolicUtils._occursin(x,~n) ? (-(~k))/( ((~n)-1)*(~a)^((~n)-1)*(polynomial_coeffs(~a,[x])[1][x]) ) : nothing
    	
    invPowerRule2 = @rule (~k)/((~q)*(~a)^(~n)) => !SymbolicUtils._occursin(x,~k) && !SymbolicUtils._occursin(x,~n) ? (-(~k))/( (~q)*((~n)-1)*(~a)^((~n)-1)*(polynomial_coeffs(~a,[x])[1][x]) ) : nothing 
	
	
	#rubi 1.1.1.2
	
	
	
	cases = [powerRule,invRule,invPowerRule,invPowerRule2]
	
	for case in cases
        ans = case(expr)
        ans != nothing && return expand(ans)
    end
	
	return nothing
end

function exponentialIntegration(expr,x)
	exponentialRule = @rule (~n)^(~a) => !SymbolicUtils._occursin(x,~n) && isone(degree(~a,x))  ? SymbolicUtils.Div((~n)^(~a),term(log,~n)*polynomial_coeffs(~a,[x])[1][x]  ) : nothing
	
end

function trigIntegration(expr,x)

end

function integrate(expr, x)
	expr = Symbolics.unwrap(expr)#for some reason it sometimes decides to wrap itself in a num type. Not sure why
	expr = rewrite_inv_to_frac(expr)
	
	#println(expr)
	
	possibleAns = algebraicIntegration(expr,x)
    if(possibleAns != nothing)
    	return possibleAns
    end
	
	if !SymbolicUtils._occursin(x,expr)
		return expr*x
    elseif isequal(expr,x)
        return x^2/2
    elseif expr isa SymbolicUtils.Mul#seperate coeff
        parts = arguments(expr)
        constants = []
        hasVar = []
        
        for part in parts
            if(SymbolicUtils._occursin(x,part))
                push!(hasVar,part)
            else
                push!(constants,part)
            end
        end
    
    	if length(constants) > 0
        	return (*(constants...))*integrate(*(hasVar...),x)
        end
    elseif expr isa SymbolicUtils.Div &&  !SymbolicUtils._occursin(x,expr.den)
    	return integrate(expr.num,x)/expr.den
    elseif expr isa SymbolicUtils.Add
        parts = arguments(expr)

        return +((integrate.(parts,x))...)
    else#integral of certain functions
        expRule = @rule exp(~a) => isone(degree(~a,x)) ? exp(~a)/( polynomial_coeffs(~a,[x])[1][x]  ) : nothing
        logRule = @rule log(~a) => isone(degree(~a,x)) ? (x*log(~a)-x) : nothing
        #trig
        sinRule = @rule sin(~a) => isone(degree(~a,x)) ? -cos(~a)/( polynomial_coeffs(~a,[x])[1][x]  ) : nothing
        cosRule = @rule cos(~a) => isone(degree(~a,x)) ? sin(~a)/( polynomial_coeffs(~a,[x])[1][x]  ) : nothing
        tanRule = @rule tan(~a) => isone(degree(~a,x)) ? -log(cos(~a))/( polynomial_coeffs(~a,[x])[1][x]  ) : nothing
        
    
        cases = [expRule,
                logRule,
                sinRule,
                cosRule,
                tanRule
                ]

        for case in cases
            ans = case(expr)
            ans != nothing && return ans
        end

    end
    
    if expr isa SymbolicUtils.Div #partial fractions integration
    	if degree(expr.den,x) > degree(expr.num,x)
    		newExpr = partial_fraction(expr,x)
    		
    		if !isequal(newExpr,expr)
    			expr = newExpr
    			return integrate(expr,x)
    		end
    	else
    		newExpr = poly_div(expr,x)
    		
    		if !isequal(newExpr,expr)
    			expr = newExpr
    			return integrate(expr,x)
    		end
    	end
    end
    
   	#possible_ans = integration_by_parts(expr,x)
   	#possible_ans != nothing && return possible_ans
   	#possible_ans = u_sub(expr,x)
   	#possible_ans != nothing && return possible_ans
    
    @warn "unable to integrate expression"
    return nothing
    #return Integrate(expr,x)
end

