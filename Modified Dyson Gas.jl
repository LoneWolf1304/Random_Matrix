
using Random, Distributions, LinearAlgebra, Plots, CSV, DataFrames

function mod_rep(X,i, p)
    s=0.0
    for j in 1:length(X)
        if rand()<p
            J = 0
        else
            J= 1
        end
        if j!=i
            s+=J/(X[i]-X[j])
        end
    end
    s
end




function mod_dys_gas(n,ti, step, p)
    dt=step
    t=dt:dt:ti
    l=length(t)+1
    P=zeros(l,n)
    P[1,:] = rand(sqrt(10^(-3))*Normal(0,1),n)
    #P[1,:] = eigvals(H_GOE(n))
    for (i,j) in enumerate(t)
        X=P[i,:]
        for m =1:n
            P[i+1,m]=P[i,m] +sqrt(2*dt/n)*randn()+sqrt(2)*dt*mod_rep(X,m, p)/n
        end
        P[i+1,: ] = sort(P[i+1,:])
    end
    P
end


function ratio_p(n,trial, dt, ti, p)
    t=ti
    r=[]
    for m=1:trial
        gas=mod_dys_gas(n, ti,dt, p)
        pf = gas[size(gas)[1],:]
        for i=8:n-8
            s=minimum([pf[i+1]-pf[i],pf[i]-pf[i-1]])/maximum([pf[i+1]-pf[i],pf[i]-pf[i-1]])
            append!(r,s)
        end
    end
    mean(r)
end
    
    
n=50
t=2
M = [ratio_p(n,1000, 0.001, t, p) for p in LinRange(0,1,100)];


# df0= CSV.read("r_vs_p=$p, n=$n.csv", DataFrame)
df1=DataFrame(n500=M)
# col = hcat(df0, df1, makeunique=true)
# CSV.write("r_vs_p=$p, n=$n.csv", df1)

CSV.write("r_vs_p, n=$n.csv", df1)


# n=50
# ti=1
# step=0.001
# p=0.8
# A = mod_dys_gas(n,ti, step, p)
# for i=8:n-8
# plot!(0:step:ti, A[1:size(A)[1], i], label="")
# end
# display(plot!(xlabel="Time", ylabel="Position"))
# savefig("p=$p.svg")