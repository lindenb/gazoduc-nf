// -l h='!(gknzwd2*)
// -l h='!(3t4jb5j*)'
process {
	executor="sge"
	clusterOptions = "-S /bin/bash -q max-24h.q  -l h='!(gkpwwd2*)'"
	cache = "lenient"
	penv = "smp"
}
