/*

Copyright (c) 2026 Pierre Lindenbaum

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
process PDF_NAVIGATION {
label "process_single"
input:
	tuple val(meta),path("PDF/*")
output:
	tuple val(meta),path("*.zip"),optional:true,emit:zip
    tuple val(meta),path("index.html"),emit:html
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:(meta.id?meta.id:"archive")
    def dir = prefix
    def with_zip = (task.ext.with_zip?:true).toBoolean()
    def title = task.ext.title?:"PDFs"
"""
hostname 1>&2
mkdir -p "${dir}"



cat << EOF > index.html
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<meta http-equiv="author" content="Pierre Lindenbaum Phd ">
<title>${title}</title>
<script>
var files=[
EOF

find PDF -name "*.pdf" -printf "\\"%f\\"\\n" | sort -V -T . | paste -sd ','  >>  index.html 

cat << EOF >> index.html
];
var page =0;

function goTo(dx) {
    if(files.length==0) return;
    page = page+dx;
    if(page<0) page = files.length-1;
    if(page>=files.length) page=0;
    document.getElementById("id1").src = files[page];
    document.getElementById("h1").textContent = files[page]+ " ("+(page+1)+"/"+files.length+")";
    }

 

window.addEventListener('load', (event) => {
  var frame = document.getElementById("id1");
  frame.style.height=(frame.contentWindow.document.body.scrollHeight+20)+'px';

  var sel = document.getElementById("select")
  for(var i=0;i< files.length;++i) {
    var opt = document.createElement("option");
    sel.appendChild(opt);
    opt.appendChild(document.createTextNode(files[i]));
    }
  sel.addEventListener("change", (event) => { goTo(sel.selectedIndex);})

  goTo(0);
});

</script>
</head>
<body>
<div>
    <button onclick="goTo(-1);">PREV</button>
    <span id="h1"></span>
    <button onclick="goTo(1);">NEXT</button>
    <select id="select"></select>
</div>

<iframe id="id1" style="height:800px;width:100%;" src="blank_">
</iframe>

<div>
    <button onclick="goTo(-1);">PREV</button>
    <span>navigation</span>
    <button onclick="goTo(1);">NEXT</button>
</div>


</body>
</html>
EOF

if ${with_zip}
then
    cp -v PDF/*.pdf ${dir}/
    cp index.html ${dir}/
    zip -9r  "${prefix}.zip" "${dir}"
    rm  -rf "${dir}"
fi

touch versions.yml
"""

stub:
	def prefix = task.ext.prefix?:(meta.id?meta.id:"archive")
"""
touch versions.yml ${prefix}.zip index.html
"""
}

