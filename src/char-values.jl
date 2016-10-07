
export Characteristicλ,
       CharacteristicA,
       CharacteristicB
""" 
Characteristicλ(q,ν,k)

Characteristic value λ_(ν+k) for Mathieu's equation 

y'' + (λ_(ν+k) - 2 q cos( 2z )) y = 0

where

q ∈ ℝ       - parameter
ν ∈ [-1,1]  - fractional part of the non-integer order 
k ∈ ℤ⁺      - range of integer parts of the order

"""
function Characteristicλ(q::Real, nu::Real)
    Characteristicλ(q, nu:nu)
end

function Characteristicλ(q::Real, k_::FloatRange)
    @assert isinteger(step(k_)) && round(Int,step(k_))==1 "Range of eigenvalue indices must have unit step"

    k = floor(Int,first(k_)):floor(Int,last(k_))
    nu_ = first(k_)-floor(Int,last(k_))
    #nu = reduced ? rem(nu_+1,2)-1 : nu_;
    nu = rem(nu_+1,2)-1;
    
    # Set matrix size using formula from Shirts paper (1993), Eqs. (2.1)-(2.2).
    nu0 = nu + maximum(k);
    C = (8.46 + 0.444*nu0)/(1 + 0.085*nu0);
    D =  (0.24 + 0.0214*nu0)/(1 + 0.059*nu0);
    N = round(Int,ceil((nu0 + 2 + C*abs(q)^D)/2)); # matrix size is 2N+1
    
    d0 = (2*[-N:N;]-nu).^2;
    d1 = q*ones(2N);
    A = SymTridiagonal(d0,d1)
    a = eigvals(A,k)
    return a
end

""" 
CharacteristicA(q, k)

Characteristic value A_k for Mathieu's equation 

y'' + (A_k - 2 q cos( 2z )) y = 0

where

q ∈ ℝ  - parameter
k ∈ ℤ⁺ - eigenvalue index, or range of indices

"""
function CharacteristicA(q::Real, k::Int)
    CharacteristicA(q,k:k)
end

function CharacteristicA(q::Real, k::UnitRange)
    @assert all(k.>=0) "Indices must be non-negative integers."

    # Boolean indices of even and odd n values
    ie = map(iseven, k);
    io = !ie;

    a = Array(Float64,length(k));
    a[ie] = Characteristicλ(abs(q),k+1)[ie];
    if q>=0
        a[io] = Characteristicλ(q,k+1)[io]; 
    else
        if 0 in k # maybe not the cleanest way to do it
            a[io] = Characteristicλ(abs(q),k[2]:last(k))[io[2:end]]
        else
            a[io] = Characteristicλ(abs(q),k)[io]; 
        end
    end
    return a
end

""" 
CharacteristicB(q,k)

Characteristic value B_k for Mathieu's equation 

y'' + (B_k - 2 q cos( 2z )) y = 0

where

q ∈ ℝ  - parameter
k ∈ ℤ  - eigenvalue index or range of eigenvalues indices

"""
function CharacteristicB(q::Real, k::Integer)
    CharacteristicB(q,k:k)
end

function CharacteristicB(q::Real, k::UnitRange)
    @assert all(k.>0) "Indices must be positive integers."
    # Boolean indices of even and odd n values
    ie = map(iseven, k);
    io = !ie;

    b = Array(Float64,length(k));
    b[ie] = Characteristicλ(q,k)[ie]; 
    if q>=0
        b[io] = Characteristicλ(q,k)[io];
    else
        b[io] = Characteristicλ(abs(q),k+1)[io]; 
    end

    return b
end


