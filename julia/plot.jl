using PyPlot
using CSV, DataFrames

cd(@__DIR__)
logname_ac = "118_ac_feasibility.csv"
data_ac = CSV.read(string(logname_ac), DataFrame, header = 0)
logname_dc = "118_dc_feasibility.csv"
data_dc = CSV.read(string(logname_dc), DataFrame, header = 0)

diff_ac = Vector{Float64}()
diff_dc = Vector{Float64}()
@assert size(data_ac, 2) == size(data_dc, 2) "Warning: Different solution files might have been used."
for i in 1:size(data_ac, 2)
    push!(diff_dc, (data_dc[1, 1] - data_dc[1, i]) / data_dc[1, 1])
    if data_ac[2, i] != -1
        push!(diff_ac, (data_ac[2, 1] - data_ac[2, i]) / data_ac[2, 1])
    else
        push!(diff_ac, 0)
    end
end

diff_ac = diff_ac .* 100
diff_dc = diff_dc .* 100

fig = figure(figsize=(5, 1.8))
rc("text", usetex = true)
rc("text.latex", preamble = "\\usepackage{amsmath}\\usepackage{calrsfs}")
rc("font",family="serif",style="italic", size = 11)
rc("mathtext",fontset="dejavuserif")
rc("lines",linewidth=1)

#ax = fig.add_axes([0.092,0.24,0.9,0.7])
ax = fig.add_axes([0.095,0.24,0.9,0.7])

ax.tick_params(direction="in",top=true,right=true,width=1.4)
grid(color="lightgray",linewidth=0.3, ls = "dashed")
xlabel("\$\\Delta z_{5}\$ [\\%]")
ylabel("\$\\Delta z_{2}\$ [\\%]")
ax.set_axisbelow(true)

for axis in ["top", "bottom", "left", "right"]
        ax.spines[axis].set_linewidth(0.0)
end

#ylim(bottom = 0.001, top=5.1)
#ylabel("gap \$[\\%]\$")

msize = 15

nf_init = false
f_init = false
for i in 1:length(diff_ac)
    if data_ac[1, i] != 1
        if !nf_init
            scatter(diff_ac[i], diff_dc[i], s = msize, edgecolor = "red", facecolor = "none", label="\$\\lnot\$ AC-feasible")
            #plot(diff_ac[i], diff_dc[i], lw = 0, mew = 0, marker =  "o", mfc = "red", label="\$\\lnot\$ AC-feasible")
            nf_init = true
        else
            scatter(diff_ac[i], diff_dc[i], s = msize, edgecolor = "red", facecolor = "none")
        end
    else
        if !f_init
            scatter(diff_ac[i], diff_dc[i], s = msize, edgecolor = "blue", facecolor = "none", label="AC-feasible")
            f_init = true
        else
            scatter(diff_ac[i], diff_dc[i], s = msize, edgecolor = "blue", facecolor = "none")
        end
    end
end

legend(loc = "lower left",fancybox=false,edgecolor="black", ncol = 1)
savefig(string("plot_118_ac_feasibility.pdf"), format = :pdf)
