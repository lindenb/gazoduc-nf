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
            if len(cols) != 5:
                raise ValueError(f"Expected 5 columns, got {len(cols)}: {line}")
            rows.append(cols)
    return rows

def update_sex(sex):
    if sex == "XX":
        return "female"
    elif sex == "XY":
        return "male"
    else:
        return sex

def process_table(rows):
    # Check for duplicate IDs
    id_to_row = {}
    for row in rows:
        id, father, mother, sex, status = row
        if id in id_to_row:
            raise ValueError(f"Duplicate id found: {id}")
        id_to_row[id] = row

    # First, update sex values
    for row in rows:
        row[3] = update_sex(row[3])

    # Second, update father/mother if absent
    ids = set(id_to_row.keys())
    for row in rows:
        id, father, mother, sex, status = row
        # father
        if father != "0" and father not in ids:
            row[1] = "0"
        # mother
        if mother != "0" and mother not in ids:
            row[2] = "0"

    # Third, check father/mother sex
    for row in rows:
        id, father, mother, sex, status = row
        if father != "0":
            father_row = id_to_row[father]
            father_sex = father_row[3]
            if father_sex == "female":
                raise ValueError(f"Father {father} for individual {id} has sex 'female'")
        if mother != "0":
            mother_row = id_to_row[mother]
            mother_sex = mother_row[3]
            if mother_sex == "male":
                raise ValueError(f"Mother {mother} for individual {id} has sex 'male'")

    return rows

def write_tsv(rows,filename):
    with open(filename,"w") as f:
        for row in rows:
            print("\t".join(row))

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.tsv>", file=sys.stderr)
        sys.exit(1)
    filename = sys.argv[1]
    rows = read_tsv(filename)
    updated_rows = process_table(rows)
    write_tsv(updated_rows,"gatk.ped")

if __name__ == "__main__":
    main()
