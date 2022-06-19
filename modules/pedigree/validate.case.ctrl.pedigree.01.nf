
process VALIDATE_CASE_CONTROL_PED_01 {
executor "local"
input:
	val(meta)
	val(pedigree)
output:
	path("rvtest.ped"),emit:rvtest
	path("validated.ped"),emit:pedigree
	path("cases.txt"),emit:cases_list
	path("controls.txt"),emit:controls_list
script:
"""
# check no space
grep -F " " "${pedigree}"| head -n 1 | cat > jeter.txt
test ! -s jeter.txt

# check num columns
awk -F '\t' '(NF<6)' "${pedigree}" > jeter.txt
cat -n jeter.txt
test ! -s jeter.txt

# check no extra white spaces
tr "\t" "${pedigree}" | awk '{S=\$1;gsub(/[ ]/,"",S); if(S!=\$1) print;}' > jeter.txt
cat -n jeter.txt
test ! -s jeter.txt

# 6th column must be case or control
awk -F '\t' '!(\$6=="case" || \$6=="control")' "${pedigree}" > jeter.txt
cat -n jeter.txt
ctest ! -s jeter.txt

# make list of cases
awk -F '\t' '(\$6=="case") {print \$2;}' "${pedigree}" | sort | uniq > cases.list

# make list of controls
awk -F '\t' '(\$6=="control") {print \$2;}' "${pedigree}" | sort | uniq > controls.list

cp -v "${pedigree}"  "validated.ped"
rm jeter.txt
"""
}

