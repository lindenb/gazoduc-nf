import sys

def read_tsv(filename):
    """Read a TSV file without a header and return a list of rows."""
    rows = []
    with open(filename, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) != 6:
                raise ValueError(f"Expected 6 columns, got {len(cols)}: {line}")
            rows.append(cols)
    return rows

def update_status(status):
    status = status.lower()
    if status == "affected":
        return "case"
    elif status == "unaffected":
        return "control"
    else:
        return status   

def update_sex(sex):
    sex2 = sex.upper()
    if sex2 == "XX" or sex2 == "F"  or sex2 == "2":
        return "female"
    elif sex2 == "XY" or sex2 == "M" or sex2 == "1":
        return "male"
    else:
        return sex

def process_table(rows):
    # Check for duplicate IDs
    id_to_row = {}
    for row in rows:
        id, father, mother, sex, status , pop = row
        if id in id_to_row:
            raise ValueError(f"Duplicate id found: {id}")
        id_to_row[id] = row

    # First, update sex values
    for row in rows:
        row[3] = update_sex(row[3])
        row[4] = update_status(row[4])
        if is_empty(row[1]):
           row[1] = "0"
        if is_empty(row[2]):
            row[2] = "0"

    # Second, update father/mother if absent
    ids = set(id_to_row.keys())
    for row in rows:
        parent = row[1]
        # father
        if is_empty(parent) or parent not in ids:
            row[1] = "0"
        # mother
        parent = row[2]
        if is_empty(parent) or parent not in ids:
            row[2] = "0"

    # Third, check father/mother sex
    for row in rows:
        id, father, mother, sex, status, pop = row
        if is_empty(father)==False:
            father_row = id_to_row[father]
            father_sex = father_row[3]
            if father_sex == "female":
                raise ValueError(f"Father {father} for individual {id} has sex 'female'")
        if is_empty(mother)==False:
            mother_row = id_to_row[mother]
            mother_sex = mother_row[3]
            if mother_sex == "male":
                raise ValueError(f"Mother {mother} for individual {id} has sex 'male'")

    return rows


def is_empty(s):
    return s=="" or s=="." or s=="0"

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.tsv>", file=sys.stderr)
        sys.exit(1)
    filename = sys.argv[1]
    rows = read_tsv(filename)
    updated_rows = process_table(rows)

    samples = [row[0] for row in rows if row[3] == "male"]
    if samples:
        with open("males.txt", "w") as f:
            for sample_id in samples:
                f.write(sample_id + "\n")
    samples = [row[0] for row in rows if row[3] == "female"]
    if samples:
        with open("females.txt", "w") as f:
            for sample_id in samples:
                f.write(sample_id + "\n")
    samples = [row[0] for row in rows if row[4] == "case"]
    if samples:
        with open("cases.txt", "w") as f:
            for sample_id in samples:
                f.write(sample_id + "\n")
    samples = [row[0] for row in rows if row[4] == "control"]
    if samples:
        with open("control.txt", "w") as f:
            for sample_id in samples:
                f.write(sample_id + "\n")
    with open("sample2collection.tsv", "w") as f:
        for row in rows:
            if is_empty(row[3]) == False:
                f.write(row[0] + "\t" + row[3] + "\n")
            if is_empty(row[4]) == False:
                 f.write(row[0] + "\t" + row[4] + "\n")
            if is_empty(row[5]) == False:
                 f.write(row[0] + "\t" + row[5] + "\n")
    # write pedigree for GATK
    with open("pedigree4gatk.ped","w") as f:
        for row in rows:
            f.write(row[0])
            f.write("\t")
            f.write(row[0])
            f.write("\t")
            f.write(row[1])
            f.write("\t")
            f.write(row[2])
            f.write("\t")
            if row[3]=="male":
                f.write("1")
            elif row[3]=="female":
                f.write("2")
            else:
                f.write("0")
            f.write("\t")
            if row[4]=="case":
                f.write("2")
            elif row[4]=="control":
                f.write("1")
            else:
                f.write("0")
            f.write("\n")
if __name__ == "__main__":
    main()
