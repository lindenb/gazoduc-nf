process GHOSTSCRIPT_MERGE {
        tag "${title}"
        label "process_short"
        conda "${moduleDir}/../../../conda/ghostscript.yml"
	afterScript "rm -rf TMP"
        input:
                tuple val(title),path("PDF/*")
        output:
                tuple val(title),path("${title}.pdf"),emit:output
        script:
		def cmd = task.ext.args?:""
        """
        hostname 1>&2
	mkdir -p TMP

        find PDF/ -type l -name "*.pdf" ${cmd} > TMP/jeter.txt

        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=TMP/jeter.pdf @TMP/jeter.txt
	
	mv -v TMP/jeter.pdf "${title}.pdf"
        """
        }
