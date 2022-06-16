process sqrt_file {
executor "local"
tag "${file(files).name}"
input:
      	val(meta)
        val(files)
output:
        path("clusters.list"),emit:clusters
	path("version.xml"),emit:version
script:
       	def min_file_split = getKeyValue(meta,"min_file_split","-1")
        """
        SQRT=`awk 'END{X=NR;if(${min_file_split} > 0 && X <= ${min_file_split}){print(X);} else {z=sqrt(X); print (z==int(z)?z:int(z)+1);}}' "${files}"`

        split -a 9 --additional-suffix=.list --lines=\${SQRT} "${files}" chunck.

      	find \${PWD} -type f -name "chunck.*.list" > clusters.list

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="Name">${task.process}</entry>
		<entry key="Description">Split file into parts</entry>
		<entry key="Input">${files}</entry>
		<entry key="N">$${SQRT}</entry>
		<entry key="Output">clusters.list</entry>
	</properties>
	EOF
        """
	}

