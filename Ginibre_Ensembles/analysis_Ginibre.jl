# using CSV, DataFrames, Plots, QuadGK

# df = CSV.read("Neighbour_spacing_N=50.csv", DataFrame)
# n=50
# x = df[!,1]
# y=df[!,2]
# global a=0.0
# for i =1:length(y)-1
#     global a+=0.5*(y[i]+y[i+1])*(x[i+1]-x[i])
# end
# print(a)

# plot(x, y./a)
# savefig("Normalised_neighbour_ginibre_n=$n.svg")
# df=DataFrame(spacing=x, density=y./a)
# CSV.write("Neighbour_spacing_N=$n, normalised.csv", df)


using CSV, DataFrames, Plots, LaTeXStrings
n=2

df1 = CSV.read("Sparse_Neighbour_spacing_n=2, p = 0, normalised, mean_space 1.csv", DataFrame)
# df2 = CSV.read("Sparse_Neighbour_spacing_n=500, p = 1000, normalised, mean_space 1.csv", DataFrame)
# df3= CSV.read("Sparse_Neighbour_spacing_n=500, p = 50000, normalised, mean_space 1.csv", DataFrame)
# df4= CSV.read("Sparse_Neighbour_spacing_n=500, p = 100000, normalised, mean_space 1.csv", DataFrame)
q=1
function data_point(df)
    return df[!,1], df[!,2]
end
# p1=10/q
# p2=1000/q
# p3=50000/q
# p4=100000/q
# p=[p1, p2, p3,p4]
p=[0]
plot()
@inbounds for i =1:1
    q=round(p[i], digits=2)
    x,y=data_point(eval(Meta.parse("df$i")))
    plot!(x,y, label="p=$q", lw=2)
end

x=LinRange(0:0.01:3)
# f(x) = (81*pi^2/(2^7)).*(x.^3).*(exp.(-(9*pi/16).*x.^2))

f(x) =@. (81*pi^2/(2^7))*(x^3)*(exp(-(9*pi/16)*x^2))


plot!(xlabel="Nearest neighbour distance s", ylabel="Density P(s)")
plot!(x, f(x), lw=2, c="red")
savefig("density_sparse,n=$n.svg")

