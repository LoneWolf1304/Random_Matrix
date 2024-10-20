using Plots, LinearAlgebra, Random, Statistics, CSV, DataFrames, LaTeXStrings, StatsBase,StatsPlots, Distributions

function comp_H(n)
    H = complex(zeros(n,n))
    @inbounds for i in 1:n
        for j in 1:n
            H[i,j] = (1/sqrt(2))*(randn()+ randn()im)
        end
    end
    return H
end


function comp_sparse(H, p)
    n=size(H)[1]
    @inbounds for i in 1:p
        m=0
        q=0
       while(true)
            m=rand(1:n)
            q=rand(1:n)
            if m!=q && H[m,q]!=0
                break
            end
        end
            H[m,q]=0
            H[q,m]=0
    end
    return H
end


f(x) = @. sqrt(1 - x^2)

n=500
p=100000
prob=p/(n*(n-1)/2)
A=eigvals(comp_sparse(comp_H(n),p))./sqrt(n);
scatter(real(A), imag(A), xlabel=L"\Re(\lambda)", ylabel=L"\Im(\lambda)", label="", color="blue")
plot!(-1:0.01:1, f(-1:0.01:1), label="", color="red", lw=2)
plot!(-1:0.01:1, -f(-1:0.01:1), label="", color="red", lw=2)
savefig("circular_p=$p.svg")



