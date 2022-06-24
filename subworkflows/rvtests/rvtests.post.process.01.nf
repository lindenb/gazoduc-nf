/** extract test name from filename */
String extractTestName(f) {
    String s = f.getName();
    if(s.endsWith(".assoc")) s=s.substring(0,s.length()-6);
    int dot  = s.lastIndexOf('.');
    s = s.substring(dot+1);
    return s;
    }


workflow postProcessRvtest {
	take:
		meta
		reference
		subtitle
		list
	main:
		R_ch = installPackages(meta)
		assoc2file_ch = list.splitCsv(header: false,sep:'\t',strip:true).
			map{T->T[0]}.
			map{F->[extractTestName(file(F)),F]}.
			groupTuple()
		
		version_ch = Channel.empty()

		group_ch = groupByTest(meta, assoc2file_ch)
		version_ch = version_ch.mix(group_ch.version)

		plot_ch = plotIt(meta, R_ch.lib, subtitle, group_ch.assoc)
		version_ch = version_ch.mix(plot_ch.version)

		zip_ch = zipIt(meta, plot_ch.plots.concat(group_ch.assoc).collect())

		version_ch = merge_version(meta,"rvtests","rvtests",version_ch.collect())
	emit:
		zip = zip_ch.zip
		version = version_ch.version
	}
