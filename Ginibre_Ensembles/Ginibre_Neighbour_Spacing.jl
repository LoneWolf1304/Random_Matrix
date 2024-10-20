using Plots, LinearAlgebra, Random, Statistics, CSV, DataFrames, LaTeXStrings, StatsBase,StatsPlots, Distributions


function comp_H(n)
    H = complex(zeros(n,n))
    @inbounds for j in 1:n
        @inbounds for i in 1:n
            H[i,j] = (1/sqrt(2))*(randn()+ randn()im)
        end
    end
    return H
end

function dist(R1,R2)
    return sqrt((R1[1]-R2[1])^2+(R1[2]-R2[2])^2)
end

function NN(A)
    s=[]
    @inbounds for i= 1:length(A)
        dista=[]
        @inbounds for j =1:length(A)
            if i!=j
                push!(dista, dist(A[i],A[j]))
            end
        end
        push!(s, dista[argmin(dista)])
    end
    return s
end


Spacing=[]
n=2
@inbounds for i =1:10000
    A=eigvals(comp_H(n));
    Evals=collect(zip(real(A), imag(A)))
    s= NN(Evals);
    append!(Spacing, s)
end




        
function XDYD(X0)
    min=minimum(X0)
    max=maximum(X0)
    n=100
    weights=zeros(n-1)
    X=LinRange(min,max,n)
    @inbounds for i in 1:length(X)-1
        @inbounds for j in X0
            if j>=X[i] && j<=X[i+1]
                weights[i]+=1
            end
        end
    end
    return X[1:end-1], weights
end





A,B= XDYD(Spacing)
df=DataFrame(spacing=A, density=B)
CSV.write("Neighbour_spacing_N=$n, non-sparse.csv", df)

df = CSV.read("Neighbour_spacing_N=$n, non-sparse.csv", DataFrame)
x = df[!,1]
y=df[!,2]
global a=0.0
for i =1:length(y)-1
    global a+=0.5*(y[i]+y[i+1])*(x[i+1]-x[i])
end
print(a)

plot(x, y./a)
savefig("Neighbour_spacing_N=$n, non-sparse, normalised.svg")
df=DataFrame(spacing=x, density=y./a)
CSV.write("Neighbour_spacing_N=$n, non-sparse, normalised.csv", df)