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
	tuple val(meta),path("DATA/*")
output:
	tuple val(meta),path("*.zip"),optional:true,emit:zip
    tuple val(meta),path("index.html"),emit:html
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:(meta.id?meta.id:"archive")
    def extension = task.ext.extension?:"pdf"
    def dir = prefix
    def with_zip = (task.ext.with_zip?:true).toBoolean()
    def title = task.ext.title?:"${meta.id?:""}"
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

find DATA -name "*.${extension}" -printf "\\"%f\\"\\n" |\\
    sort -V -T . |\\
    paste -sd ','  >>  index.html

cat << EOF >> index.html
];
var the_page_index =0;


function updateSamplesheet() {
    var s="";
    for(var i in files) {
        item = files[i];
        if(item.status==null) continue;
        s+= item.file +"\t"+item.status+"\t"+item.comment+"\\n";
        }
    document.getElementById("samplesheet").textContent=s;
    }

function moveR(shift,status) {
    goTo(the_page_index + shift,status)
    }

function goTo(index,status) {
    if(files.length==0 ) return;

    var txfield = document.getElementById("comment");
    var item = files[the_page_index];
    if(status!=null) item.status = status;
    item.comment = txfield.value;

    index = index%files.length;
        
    the_page_index = index;

    
    item = files[index];
    
    txfield.value = item.comment;
    
    var E = document.getElementById("theimg");
    E.setAttribute("src",item.file);
    
    document.getElementById("navigation").textContent = item.file + " "+(index+1)+"/"+files.length+" "+(item.status==null?"":" Status:"+ item.status);

    updateSamplesheet();
 }
 

window.addEventListener('load', (event) => {
  var frame = document.getElementById("theimg");
  frame.style.height=(frame.contentWindow.document.body.scrollHeight+20)+'px';

  for(var i in files) {
    files[i]={"file":files[i],"status":null,"comment":""};
    }
  
  var sel = document.getElementById("select")
  for(var i in files) {
    var opt = document.createElement("option");
    sel.appendChild(opt);
    opt.appendChild(document.createTextNode(files[i].file));
    }
  sel.addEventListener("change", (event) => { goTo(sel.selectedIndex,null);})

  goTo(0,null);
});

</script>
</head>
<body>
<div>
    <button onclick="moveR(-1,null);">&#x1F519;</button>
    &nbsp;
    <span id="navigation">1/x</span>
    &nbsp;
    <button onclick="moveR(1,null);">&#128284;</button>
    &nbsp;
    <select id="select"></select>
    &nbsp;
    <label for="comment">Comment:</label><input id="comment"></input>
    &nbsp;
    <button style="background-color: #008CBA;" onclick="moveR(1,'AMBIGOUS')">AMBIGOUS</button>
    <button style="background-color: #f44336;" onclick="moveR(1,'BAD')">BAD</button>
    <button style="background-color: #04AA6D;" onclick="moveR(1,'OK')">OK</button>
</div>

<iframe id="theimg" style="height:800px;width:100%;" src="blank_">
</iframe>

<div>
    <button onclick="moveR(-1,null)">&#x1F519;</button>
    &nbsp;
    <span></span>
    &nbsp;
    <button onclick="moveR( 1,null)">&#128284;</button>
</div>
<br/>
<pre style="background-color: gray;" id="samplesheet"></pre>
</body>
</html>
EOF

if ${with_zip}
then
    cp -v DATA/*.${extension} ${dir}/
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

