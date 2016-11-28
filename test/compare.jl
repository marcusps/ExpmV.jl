
master = load("master-bench.jld")["bench"]
optim = load("optimize-bench.jld")["bench"]

expmv_t = [ratio(median(optim[i]["expmv"]), median(master[i]["expmv"])).time for i in 1:6]

expmv_m = [ratio(median(optim[i]["expmv"]), median(master[i]["expmv"])).memory for i in 1:6]

time = reduce(hcat,[[optim[i]["d"],median(optim[i]["expm"]).time,median(optim[i]["expokit"]).time,median(optim[i]["expmv"]).time] for i in 1:6])'

memory = reduce(hcat,[[optim[i]["d"],median(optim[i]["expm"]).memory,median(optim[i]["expokit"]).memory,median(optim[i]["expmv"]).memory] for i in 1:6])'
