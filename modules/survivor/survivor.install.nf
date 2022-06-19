process SURVIVOR_INSTALL {
input:
	val(meta)
output:
	path("SURVIVOR/Debug/survivor"),emit:executable
	path("version.xml"),emit:version
script:
"""
git clone "https://github.com/fritzsedlazeck/SURVIVOR.git"
cd SURVIVOR/Debug && make
"""
}
