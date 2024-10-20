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

function comp_sparse(H, p)
    n=size(H)[1]
    @inbounds for i in 1:n
       while(true)
            p=rand(1:n)
            q=rand(1:n)
            if p!=q && H[p,q]!=0
                break
            end
        end
            H[p,q]=0
            H[q,p]=0
    end
    return H
end

f(x) = sqrt.(1 .- x.^2)

n=1000
A=eigvals(comp_H(n))./sqrt(n);
Evals=collect(zip(real(A), imag(A)))
scatter(real(A), imag(A), xlabel=L"\Re(\lambda)", ylabel=L"\Im(\lambda)", label="", color="blue")
plot!(-1:0.01:1, f(-1:0.01:1), label="", color="red", lw=2)
plot!(-1:0.01:1, -f(-1:0.01:1), label="", color="red", lw=2)
savefig("circular.svg")

# function dist(R1,R2)
#     return sqrt((R1[1]-R2[1])^2+(R1[2]-R2[2])^2)
# end

# function NN(A)
#     s=[]
#     for i= 1:length(A)
#         dista=[]
#         for j =1:length(A)
#             if i!=j
#                 push!(dista, dist(A[i],A[j]))
#             end
#         end
#         push!(s, dista[argmin(dista)])
#     end
#     return s
# end
        

# s= NN(Evals);
# savefig("neighbour_ginibre.svg")
