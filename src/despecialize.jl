struct __Despecialize1 end
struct __Despecialize2 <: Real end
struct __Despecialize3 <: AbstractVector{Int} end

hide_lhs(::__Despecialize1) = false
hide_lhs(::__Despecialize2) = false
hide_lhs(::__Despecialize3) = false
