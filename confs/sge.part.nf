process {
	executor="sge"
	clusterOptions = "-S /bin/bash -q max-24h.q"
	cache = "lenient"
	penv = "smp"
}
