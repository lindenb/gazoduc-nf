// -l h='!(gknzwd2*)
// -l h='!(3t4jb5j*)'
// -l h='!(gkr1xd2*)'

process {
	executor="sge"
	clusterOptions = "-S /bin/bash -q max-24h.q -l h='!(2t4jb5j*)' "
	cache = "lenient"
	penv = "make"
}
