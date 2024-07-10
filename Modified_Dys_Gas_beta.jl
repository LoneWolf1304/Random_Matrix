
using Random, Distributions, LinearAlgebra, Plots, CSV, DataFrames



function mod_rep(X,i, p, J)
    s=0.0
    for j in 1:length(X)
        if j!=i
            s+=J[i,j]/(X[i]-X[j])
        end
    end
    s
end

function J(n, p)
    J=ones(n,n)
    for i =1:n 
        for j=1:i-1 
            if rand()<p
                J[i,j] = 0
                J[j,i] = 0
            end
        end
    end
    return J
end

function mod_dys_gas(n,ti, step, p, beta)
    dt=step
    t=dt:dt:ti
    l=length(t)+1
    P=zeros(l,n)
    P[1,:] = rand(sqrt(10^(-3))*Normal(0,1),n)
    #P[1,:] = eigvals(H_GOE(n))
    for (i,j) in enumerate(t)
        X=P[i,:]
        Jt= J(n,p)
        for m =1:n
            P[i+1,m]=P[i,m] +sqrt(2*dt/(beta*n))*randn()+sqrt(2)*dt*mod_rep(X,m, p, Jt)/n
        end
        P[i+1,: ] = sort(P[i+1,:])
    end
    P
end


function ratio_p(n,trial, dt, ti, p, beta)
    t=ti
    r=[]
    for m=1:trial
        gas=mod_dys_gas(n, ti,dt, p, beta)
        pf = gas[size(gas)[1],:]
        for i=8:n-8
            s=minimum([pf[i+1]-pf[i],pf[i]-pf[i-1]])/maximum([pf[i+1]-pf[i],pf[i]-pf[i-1]])
            append!(r,s)
        end
    end
    mean(r)
end
    
N=40
beta=2
M = [ratio_p(N,500, 0.001, 1, p,beta) for p in LinRange(0,1,70)];


# df0= CSV.read("r_vs_p_$N.csv", DataFrame)
# df1=DataFrame(n=M)
# col = hcat(df0, df1, makeunique=true)
# CSV.write("r_vs_p_$N.csv", col)

df1=DataFrame(n=M)
CSV.write("r_vs_p_symm coupling_N=$N beta=$beta.csv", df1)

