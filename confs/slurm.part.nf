process {
executor="slurm"
clusterOptions = "--partition=Bird "
cache="lenient"
maxForks=100
errorStrategy= "finish"
}

