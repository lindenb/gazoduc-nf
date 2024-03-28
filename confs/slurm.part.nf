// --exclude=bigmem001
// --exclude=bird002
process {
executor="slurm"
clusterOptions = "--partition=Bird"
cache="lenient"
maxForks=100
errorStrategy= "finish"
}

