struct Polynomial{C_T, M <: AbstractMonomial}
    monomials::Vector{M}
    coeffs::Vector{C_T}
    monoid::AbstractMonoid
    coeff_type::Type{C_T}
        # Constructor to filter zero coefficients
        function Polynomial(monomials::Vector{M}, coeffs::Vector{C_T}, monoid::AbstractMonoid) where {C_T, M <: AbstractMonomial}
            # Filter out zero coefficients and corresponding monomials
            filtered = filter(x -> x[2] != 0, collect(zip(monomials, coeffs)))
            new_monomials = [x[1] for x in filtered]
            new_coeffs = [x[2] for x in filtered]
            new{C_T, M}(new_monomials, new_coeffs, monoid,C_T)
        end
end

function Polynomial(m::AbstractMonoid)
    monomial_type=typeof(one(m))
    return Polynomial(Vector{monomial_type}(),Vector{Number}(),m)
end
function Polynomial{C_T}(m::AbstractMonoid) where C_T
    monomial_type=typeof(one(m))
    return Polynomial(Vector{monomial_type}(),Vector{C_T}(),m)
end

    
function Polynomial(m::AbstractMonomial,n)
    p=Polynomial([m],[n],m.monoid)
    return p
end
function Polynomial(m::AbstractMonomial)
    p=Polynomial([m],[1],m.monoid)
    return p
end



Polynomial(p::Polynomial)=p
Polynomial(v::Variable)=Polynomial(monomial(v))
Polynomial(v::Variable,n)=Polynomial(monomial(v),n)
Base.hash(p::Polynomial, h::UInt)=hash(p.monoid, hash(p.monomials, hash(p.coeffs,hash(0x7d6979235cb005d0, h))))
Base.zero(p::Polynomial)=Polynomial(p.monoid)


Base.iterate(p::Polynomial)=iterate(zip(p.monomials,p.coeffs))
Base.iterate(p::Polynomial, state)=iterate(zip(p.monomials,p.coeffs),state)
Base.length(p::Polynomial)=length(p.monomials)
function Base.show(io::IO, mime::MIME"text/plain", p::Polynomial)
    return _show(io, mime, p)
end
function Base.show(io::IO, mime::MIME"text/print", p::Polynomial)
    return _show(io, mime, p)
end
Base.print(p::Polynomial)=_show(stdout,MIME"text/plain",p)
function _show(io::IO, mime, p::Polynomial)
    if p==0
        print(io, "0")
        return
    end
    for i in 1:length(p.monomials)
        coeff = p.coeffs[i]
        # Convert rational coefficients to float and format to 4 decimal points
        coeff = coeff isa Rational ? Float64(coeff) : coeff
        coeff = coeff isa Float64 ? round(coeff, digits=5) : coeff

        if i == 1
            # For the first term, handle the sign explicitly
            if coeff isa Real 
                if coeff < 0
                    print(io, "-", abs(coeff), p.monomials[i])
                else
                    print(io, coeff, p.monomials[i])
                end
            else
                print(io, "(",coeff, ")" ,p.monomials[i])
            end
        else
            # For subsequent terms, include "+" or "-" as appropriate
            if coeff isa Real
                if coeff >= 0
                    print(io, " + ", coeff, p.monomials[i])
                elseif coeff < 0
                    print(io, " - ", abs(coeff), p.monomials[i])
                end
            else
                print(io, " + (", coeff, ")", p.monomials[i])
            end
        end
    end
end

monomials(p::Polynomial)=p.monomials

function is_number(p::Polynomial)
    return iszero(p) || (length(p.monomials) == 1 && is_identity(p.monomials[1]))
end

is_number(x::Number) = true
is_number(x) = false
is_number(m::AbstractMonomial) = is_identity(m)

function coefficient(p::Polynomial, t::AbstractMonomial)
    coeff_index = findfirst(==(t), p.monomials)
    if isnothing(coeff_index)
        return 0
    else
        return p.coeffs[coeff_index]
    end
end

function coefficient(p::Polynomial, T)
    return [coefficient(p, t) for t in T]
end

# function terms(p::Polynomial)
#     return zip(p.coeffs,p.monomials)
# end

function Base.:(==)(p::Polynomial,q::Polynomial)
    if (p.monoid!=q.monoid)
        p,q=get_root_polynomials(p,q)
        return ==(p,q)
    end
    return (p.monomials==q.monomials) && (p.coeffs==q.coeffs)
end
Base.:(==)(p::Polynomial,m::AbstractMonomial)= ==(p,Polynomial(m))
Base.:(==)(m::AbstractMonomial,p::Polynomial)= ==(Polynomial(m),p)
Base.:(==)(p::Polynomial,q::Number)=((q==0) && iszero(p)) || (!iszero(p) && is_number(p) && p.coeffs[1]==q)
Base.:(==)(q::Number,p::Polynomial)= p==q


function Base.conj(p::Polynomial{C_T}) where C_T
    T=MA.promote_operation(conj,C_T)
    q=Polynomial{T}(p.monoid)
    for i in 1:length(p.monomials)
        push!(q.monomials,conj(p.monomials[i]))
        push!(q.coeffs,conj(p.coeffs[i]))
    end
    return q
end
Base.adjoint(p::Polynomial{C_T}) where C_T = Base.conj(p)
Base.:*(x,m::AbstractMonomial)=x==0 ? zero(Polynomial(m,x)) : Polynomial(m,x)
Base.:*(m::AbstractMonomial,x)=x*m

function Base.:+(x::C_T_X, p::Polynomial{C_T_P}) where {C_T_X,C_T_P}
    T=MA.promote_operation(*,C_T_X,C_T_P)
    scaled_poly = Polynomial{T}(p.monoid)  # Create a new Polynomial with the same monoid
    for index in 1:length(p.coeffs)
        res_coeff = T(1) * p.coeffs[index]  # Scale the coefficient by x
        res_coeff == 0 && continue  # Skip if the coefficient is zero
        push!(scaled_poly.monomials, p.monomials[index])  # Copy the monomial
        push!(scaled_poly.coeffs, res_coeff)  # Scale the coefficient by x
    end
    return scaled_poly + x*one(p.monoid)
end

function Base.:*(x::C_T_X, p::Polynomial{C_T_P}) where {C_T_X,C_T_P}
    T=MA.promote_operation(*,C_T_X,C_T_P)
    scaled_poly = Polynomial{T}(p.monoid)  # Create a new Polynomial with the same monoid
    for index in 1:length(p.coeffs)
        res_coeff = x * p.coeffs[index]  # Scale the coefficient by x
        res_coeff == 0 && continue  # Skip if the coefficient is zero
        push!(scaled_poly.monomials, p.monomials[index])  # Copy the monomial
        push!(scaled_poly.coeffs, res_coeff)  # Scale the coefficient by x
    end
    return scaled_poly
end

Base.:+(p::Polynomial,x)=x+p

Base.:*(p::Polynomial,x)=x*p

Base.:/(p::Polynomial,f::FE) where {FE<:Number}= (1/f)*p
    

Base.:+(p::Polynomial,q::Polynomial)=add_poly(p,q)      

general_add(m::AbstractMonomial,n::AbstractMonomial)=Polynomial(m)+Polynomial(n)

    

function add_poly(p::Polynomial{C_T_1}, q::Polynomial{C_T_2}) where {C_T_1,C_T_2}
    if p.monoid == q.monoid
        T=MA.promote_operation(+,C_T_1,C_T_2)
        r = Polynomial{T}(p.monoid)  # Create a new polynomial with the same monoid

        # cmn_mons=intersect(p.monomials,q.monomials)
        # for mon in cmn_mons
        #     push!(r.monomials,mon)
        #     index_p=findfirst(m->m==mon,p.monomials)
        #     index_q=findfirst(m->m==mon,q.monomials)
        #     push!(r.coeffs,p.coeffs[index_p]+q.coeffs[index_q])
        # end
        # for i in 1:length(p)
        #     if !(p.monomials[i] in cmn_mons)
        #         push!(r.monomials,p.monomials[i])
        #         push!(r.coeffs,p.coeffs[i])
        #     end
        # end
        # for j in 1:length(q)
        #     if !(q.monomials[j] in cmn_mons)
        #         push!(r.monomials,q.monomials[j])
        #         push!(r.coeffs,q.coeffs[j])
        #     end
        # end
        i, j = 1, 1  # Pointers for p and q

        # Merge the monomials and coefficients
        while i <= length(p.monomials) && j <= length(q.monomials)
            if p.monomials[i] == q.monomials[j]
                # Monomials are the same, add coefficients
                coeff = p.coeffs[i] + q.coeffs[j]
                if coeff != 0
                    push!(r.monomials, p.monomials[i])
                    push!(r.coeffs, coeff)
                end
                i += 1
                j += 1
            elseif p.monomials[i] > q.monomials[j]
                # Monomial in p is larger, add it to r
                push!(r.monomials, p.monomials[i])
                push!(r.coeffs, p.coeffs[i])
                i += 1
            else
                # Monomial in q is larger, add it to r
                push!(r.monomials, q.monomials[j])
                push!(r.coeffs, q.coeffs[j])
                j += 1
            end
        end

        # Add remaining monomials and coefficients from p
        while i <= length(p.monomials)
            push!(r.monomials, p.monomials[i])
            push!(r.coeffs, p.coeffs[i])
            i += 1
        end

        # Add remaining monomials and coefficients from q
        while j <= length(q.monomials)
            push!(r.monomials, q.monomials[j])
            push!(r.coeffs, q.coeffs[j])
            j += 1
        end

        return r
    else
        # Handle the case where monoids differ
        return sum(get_root_polynomials(p, q))
    end
end

function mult_poly(p::Polynomial{C_T_1}, q::Polynomial{C_T_2}) where {C_T_1,C_T_2}
    if q==0 || p==0
        return zero(p)
    end
    T= MA.promote_operation(*,C_T_1,C_T_2)
    res_vec = Vector{Polynomial{T}}()  # Vector to store monic polynomials
    # Multiply all monomials in p with all monomials in q
    for i in 1:length(p.monomials)
        for j in 1:length(q.monomials)
            # Multiply monomials and coefficients
            new_monomial = p.monomials[i] * q.monomials[j]

            new_coeff = p.coeffs[i] * q.coeffs[j]
            (new_coeff==0 || new_monomial==0) && continue  # Skip if coefficient is zero
            # Create a monic polynomial for the result
            if isa(new_monomial, Polynomial)
                push!(res_vec,new_coeff*new_monomial)
            else
                monic_poly = Polynomial{T}(new_monomial.monoid)
                push!(monic_poly.monomials, new_monomial)
                push!(monic_poly.coeffs, new_coeff)

                # Store the monic polynomial in the result vector
                push!(res_vec, monic_poly)
            end
        end
    end

    # Sum all the monic polynomials in the result vector
    return sum(res_vec;init=Polynomial{T}(p.monoid))
end
Base.one(p::Polynomial)=one(p.monoid)

Base.:+(M::AbstractMonomial,n::Number)=+(M,n*one(M))
Base.:+(n::Number,M::AbstractMonomial)=+(n*one(M),M)
Base.:+(m1::AbstractMonomial,m2::AbstractMonomial)=1*m1+1*m2

Base.:-(m1::AbstractMonomial,m2::AbstractMonomial)=1*m1-1*m2
Base.:-(m::AbstractMonomial,n::Number)=m-n*one(m)
Base.:-(n::Number,m::AbstractMonomial)=n*one(m)-m

Base.:*(p::Polynomial,q::Polynomial)=mult_poly(p,q)
Base.:*(p::Polynomial,x::AbstractMonomial)=mult_poly(p,Polynomial(x))
Base.:*(x::AbstractMonomial,p::Polynomial)=mult_poly(Polynomial(x),p)
Base.:+(p::Polynomial,x::AbstractMonomial)=add_poly(p,Polynomial(x))
Base.:+(x::AbstractMonomial,p::Polynomial)=add_poly(Polynomial(x),p)
Base.:-(p::Polynomial,q::Polynomial)=add_poly(p,-q)
Base.:-(p::Polynomial,x::AbstractMonomial)=add_poly(p,-Polynomial(x))
Base.:-(x::AbstractMonomial,p::Polynomial)=add_poly(Polynomial(x),-p)
Base.:+(p::Polynomial,x::Number)=add_poly(p,Polynomial(one(p.monoid),x))
Base.:+(x::Number,p::Polynomial)=add_poly(Polynomial(one(p.monoid),x),p)
Base.:-(p::Polynomial,x::Number)=add_poly(p,Polynomial(one(p.monoid),-x))
Base.:-(x::Number,p::Polynomial)=add_poly(Polynomial(one(p.monoid),x),-p)
Base.:-(p::Polynomial)=-1*p


MA.promote_operation(::F,::Type{Number},::Type{N}) where {F<:Function,N<:Number} =Number
MA.promote_operation(::F,::Type{N},::Type{Number}) where {F<:Function,N<:Number} =Number
MA.promote_operation(::F, ::Type{Number}, ::Type{Number}) where F<:Function = Number

variables(p::Polynomial)=Vector{Variable}(reduce(union, [variables(i) for i in p.monomials];init=Variable[]))

function real_rep(p::Polynomial{C_T}) where C_T
    result = Polynomial{C_T}(p.monoid)

    for (m,c) in p
        (rm, s, im) = conj_minmax(m)
        result+=real(c)*rm
        result+=s*imag(c)*im
    end

    return result
end

real_rep(m::AbstractMonomial) = conj_min(m)

function degree(p::Polynomial)
    return maximum([degree(m) for m in p.monomials])
end

Base.:^(e::A,p::Int) where {A<: Union{AbstractMonomial,Polynomial}}= e==0 ? one(e) : prod([e for _ in 1:p])