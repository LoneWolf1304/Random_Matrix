
using LinearAlgebra
using Plots
using Distributions, Random, Statistics, StatsPlots
using CSV, DataFrames

#defining GOE hamiltonian of dimension nxn
function H_GOE(n)
    A=rand(Normal(0,1), n,n)
    Q=(A+transpose(A))/2
    Q/sqrt(n)
end
     

#defining GUE hamiltonian of dimension nxn
function H_GUE(n)
    H=complex(zeros(n,n));
    for i = 1 : n
        for j = 1 : n
            H[i,j]=rand(Normal(0,1))+ rand(Normal(0,1))im
        end
    end
    #Making a hermitian matrix from H
    Hs=(H+H')/2
    Hs= Hs.*sqrt(n/tr(Hs*Hs))
    Hs
end

function spacings_ratio(n, b)
    l=[]
    A=zeros(n,n)
    if b==1
        A = H_GOE(n)
    else
        A=H_GUE(n)
    end
        E=real(eigvals(A))
        p=sort(E)
        for i =2:length(p)-1
            s=(p[i+1]-p[i])/((p[i]-p[i-1]))
            append!(l,s)
        end
    l
end
#function which returns array of spacing ratio r~ for GOE
function spacings_ratio_tilde(n,b)
    l=[]
    A=zeros(n,n)
    if b==1
        A = H_GOE(n)
    else
        A=H_GUE(n)
    end
        E=real(eigvals(A))

        p=sort(E)
        for i =2:length(p)-1
            r_ti=minimum([p[i+1]-p[i],p[i]-p[i-1]])/maximum([p[i+1]-p[i],p[i]-p[i-1]])
            append!(l,r_ti)
        end
    l
end




for cntt =1:4000
n=10
t=100
b=2
arr=["GOE", "GUE"]
q=arr[b]
l=spacings_ratio(n,b)
df = DataFrame(l1=l)
CSV.write("N=$n $q ratio $cntt.csv", df)

for i=1:t
	df1=CSV.read("N=$n $q ratio $cntt.csv", DataFrame)
	l=spacings_ratio(n,b)
	col=hcat(df1,l, makeunique=true)
	CSV.write("N=$n $q ratio $cntt.csv",col)
end

end


