/*

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
include {CIRCOS as APPLY_CIRCOS                 } from '../../modules/circos/apply'
include {CYTOBAND_TO_KARYOTYPE                  } from '../../modules/circos/cytoband2karyotype'
include {DOWNLOAD_CYTOBAND                      } from '../../modules/ucsc/download.cytobands'
include {isBlank                                } from '../../modules/utils/functions'
include {IF_EMPTY as IF_EMPTY1                  } from '../../subworkflows/nf/if_empty'


workflow CIRCOS {
    take:
        workflow_metadata
        fasta
        fai
        dict
        cytobands //or empty channel
        data_channel // heterogenous source of data, first must be a meta with 'circos' property. One one '(bamcov)'
    main: 
        versions = Channel.empty()


        //download cytobands if needed
        DOWNLOAD_CYTOBAND(
            cytobands
                .count()
                .combine(fasta)
                .filter{count,_meta,_fasta->count==0}
                .map{_count,meta,fasta->[meta,fasta]},
            fai,
            dict
            )
        versions = versions.mix(DOWNLOAD_CYTOBAND.out.versions)
        cytobands = IF_EMPTY1(cytobands, DOWNLOAD_CYTOBAND.out.bed)

        circos_data = Channel.empty()
        circos_conf = Channel.empty()

        dispatch= data_channel.branch{v->
            blank: isBlank(v[0].tracktype)
            bamcov : v[0].tracktype.equals("bamcov")
            other:true
            }

    	dispatch.blank
		    .mix(dispatch.other)
            .map{throw new IllegalArgumentException("unknow or undefined tracktype in ${it}")}


        BAMCOV_TO_CIRCOS(dispatch.bamcov)
        circos_data = circos_data.mix(BAMCOV_TO_CIRCOS.out.datafile.map{it[1]})
        circos_conf = circos_conf.mix(BAMCOV_TO_CIRCOS.out.datafile.map{it[0].plus(filename:it[1].name)})

        CYTOBAND_TO_KARYOTYPE(
            fai,
            cytobands.ifEmpty([[id:"noctyoband"],[]])
            )
        versions = versions.mix(CYTOBAND_TO_KARYOTYPE.out.versions)

        
        
         plot_string = circos_conf
            .collect(flat:false)
            .flatMap{v->
                def L=[];
                def n= v.size();
                def sum_weight=0.0;
                for(x=0;x< n;x++) {
                    def vx = v[x];
                    def weight = ((vx.weight?:1.0) as double)
                    sum_weight += weight;
                    }
                def r_max= 0.9;
                def r_min= 0.4;
                def r=r_max;
                for(x=0;x< n;x++) {
                    def vx = v[x];
                    def weight = ((vx.weight?:1.0) as double)
                    def dr = (weight/sum_weight)*(r_max-r_min);
                    def order = vx.order?:java.lang.String.format("%03d",x+1);
                    L.add(vx.plus([
                        circos_index:x,
                        r0:(r-dr*0.98),
                        r1:r,
                        order: order
                        ]));
                    r-=dr;
                    }
                return L.sort{a,b->a.order.compareTo(b.order)};
                }
            .map{
                """
                <plot>
                type = ${it.tracktype.matches("bamcov")?"histogram":"todo"}
                file = ${it.filename}
                color = ${it.color?:"grey"}
                orientation = ${it.orientation?:"in"}
                r0 = ${it.r0}r
                r1 = ${it.r1}r
                ${it.min?"min = ${it.min}":""}
                ${it.max?"max = ${it.max}":""}
                fill_color = ${it.fill_color?:"green"}
                thickness = ${it.thickness?:"1"}
                </plot>
                """
                }
            .collect()
            .map{ "<plots>\n"+it.sort().join("\n")+"</plots>\n"}
            .map{[[id:"plots"],it]}
            .view()
       // circos_conf = circos_conf.mix(IDEOGRAM_CONF.out.config)
        //circos_data = circos_data.mix(IDEOGRAM_CONF.out.data)

        MAKE_CONF(
            fai,
            CYTOBAND_TO_KARYOTYPE.out.karyotype,
            plot_string
            )


        APPLY_CIRCOS(
            MAKE_CONF.out.config,
            CYTOBAND_TO_KARYOTYPE.out.karyotype,
            circos_data.collect().map{[[id:"circos"],it]}
            )
        versions = versions.mix(APPLY_CIRCOS.out.versions)
        
    emit:
        versions
        svg = APPLY_CIRCOS.out.svg
        png = APPLY_CIRCOS.out.png
}

process MAKE_CONF {
tag "${meta.id}"
label "process_single"
input:
    tuple val(meta),path(fai)
    tuple val(meta2),path(karyotype)
    tuple val(meta2),val(plot_str)
output:
    tuple val(meta),path("*.conf"),emit:config
script:
    def regex=task.ext.regex?:"^(chr)?[0-9XY]+\$"
"""

cat << EOF > jeter.conf
karyotype = ${karyotype.name}
chromosomes_units = 1000000
chromosomes_display_default = false

<ideogram>
<spacing>
default = 10u
</spacing>

	## basic
	radius = 0.8r
	thickness = 80p
	fill = yes
	stroke_thickness = 2p
	stroke_color = black

	## label
	show_label = yes
	label_radius = 1r + 20p
	label_size = 10
	label_parallel = yes
	label_color = 
	label_format = eval(uc var(label))

	## bands
	show_bands = yes
	fill_bands = yes
	band_stroke_thickness = 0.5
	band_stroke_color     = vlgray
	band_transparency = 3
</ideogram>
<image>
<<include etc/image.conf>>
</image>
EOF

cut -f 1 ${fai} |\\
    grep -E '${regex}' |\\
    paste -s -d ';' |\\
    sed 's/^/chromosomes = /' >> jeter.conf

cat << EOF >> jeter.conf
<<include etc/colors_fonts_patterns.conf>>
EOF

cat  << 'EOF' >> jeter.conf
${plot_str}
EOF

cat << EOF >> jeter.conf
<<include etc/housekeeping.conf>>
data_out_of_range* = trim
EOF

mv jeter.conf karyotype.${meta.id}.conf

"""
}

process BAMCOV_TO_CIRCOS {
tag "${meta.id}"
label "process_single"
input:
    tuple val(meta),path(cov)
output:
    tuple val(meta),path("*.txt"),emit:datafile
script:
"""
awk -F '\t' '{printf("%s\t%s\t%s\t%f\\n",\$1,\$2,\$3,int(\$4)/(int(\$3)-int(\$2)));}' ${cov} > ${cov.name}.circos.txt
"""
}