process DOWNLOAD_ORAD {
label 'process_single'
afterScript "rm -rf jeter.tar.gz"
input:
        val(meta)
output:
        tuple val(meta),path("orad.*.linux"),emit:oradir
        path("versions.yml"),emit:versions
script:
        def version = task.ext.version?:"2.7.0"
        def url = task.ext.url?:"https://s3.amazonaws.com/webdata.illumina.com/downloads/software/dragen-decompression/orad.${version}.linux.tar.gz";
"""
curl -L -o jeter.tar.gz "${url}"
tar xvfz jeter.tar.gz
rm jeter.tar.gz


cat <<-END_VERSIONS > versions.yml
"${task.process}":
    orad: ${version}
    url: ${url}
END_VERSIONS

"""
stub:
""" 
touch orad.stub.linux
touch versions.yml
"""
}
