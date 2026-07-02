abstract type AbstractMonoid end
abstract type AbstractMonomial end

Base.isless(a::AbstractMonoid, b::AbstractMonoid) = isless(a.name, b.name)

Base.show(io::IO, mime::MIME"text/plain", monoid::AbstractMonoid) = print(io, monoid.name)
Base.show(io::IO, mime::MIME"text/print", monoid::AbstractMonoid) = print(io, monoid.name)
Base.show(
    io::IO,
    mime::MIME"text/plain",
    monoid::Tuple{Vararg{Monoid}},
) where {Monoid <: AbstractMonoid} = print(io, Tuple((m.name for m in monoid)))
Base.show(
    io::IO,
    mime::MIME"text/print",
    monoid::Tuple{Vararg{Monoid}},
) where {Monoid <: AbstractMonoid} = print(io, Tuple((m.name for m in monoid)))
# Base.print(io::IO, mime::MIME"text/plain", monoid::Tuple{Vararg{Monoid}}) where Monoid<:AbstractMonoid = print(io,Tuple((m.name for m in monoid)))
# Base.print(io::IO, mime::MIME"text/print", monoid::Tuple{Vararg{Monoid}}) where Monoid<:AbstractMonoid = print(io,Tuple((m.name for m in monoid)))
# Base.show(io::IO, mime::MIME"text/plain", monoid::Vector{Monoid}) where Monoid<:AbstractMonoid = print(io,[m.name for m in monoid])
# Base.show(io::IO, mime::MIME"text/print", monoid::Vector{Monoid}) where Monoid<:AbstractMonoid = print(io,[m.name for m in monoid])
