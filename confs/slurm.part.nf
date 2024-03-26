// --exclude=bigmem001
process {
executor="slurm"
clusterOptions = "--partition=Bird --exclude=bird002"
cache="lenient"
maxForks=100
errorStrategy= "finish"
}

