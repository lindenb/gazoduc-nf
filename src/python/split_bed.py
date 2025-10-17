import sys
import hashlib

class Interval:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
    def __repr__(self):
        return f"Interval({self.chrom}, {self.start}, {self.end})"
    def length(self):
        return self.end - self.start
    def split(self, split_length):
        """Split this interval into two at split_length (relative to start).
        Returns (part1, part2) where part1 is of size split_length, part2 is the remainder (or None if not needed)."""
        if split_length >= self.length():
            return (self, None)
        part1 = Interval(self.chrom, self.start, self.start + split_length)
        part2 = Interval(self.chrom, self.start + split_length, self.end)
        return (part1, part2)

class Cluster:
    def __init__(self, intervals=None):
        self.intervals = intervals if intervals is not None else []
    def add(self, interval):
        self.intervals.append(interval)
    def total_length(self):
        return sum(iv.length() for iv in self.intervals)
    def _md5_of_intervals(self):
        bed_string = ''.join(f"{iv.chrom}\t{iv.start}\t{iv.end}\n" for iv in self.intervals)
        return hashlib.md5(bed_string.encode('utf-8')).hexdigest()
    def save_as_bed(self,prefix):
        filename = prefix + self._md5_of_intervals() + ".bed"
        with open(filename, "w") as out:
            for iv in self.intervals:
                out.write(f"{iv.chrom}\t{iv.start}\t{iv.end}\n")
        return filename

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <factor>  <prefix> <bed_file>", file=sys.stderr)
        sys.exit(1)

    try:
        factor = int(sys.argv[1])
        if factor <= 1:
            raise ValueError
    except ValueError:
        print("Error: 'factor' must be an integer greater than one.", file=sys.stderr)
        sys.exit(1)
    prefix = sys.argv[2]
    bed_path = sys.argv[3]
    intervals = []
    try:
        with open(bed_path, "r") as bed_file:
            for line in bed_file:
                if line.strip() == "" or line.startswith("#"):
                    continue
                fields = line.rstrip().split()
                if len(fields) < 3:
                    continue # skip malformed lines
                chrom, start, end = fields[0], fields[1], fields[2]
                intervals.append(Interval(chrom, start, end))
    except Exception as e:
        print(f"Error reading BED file: {e}", file=sys.stderr)
        sys.exit(1)

    # Compute genome length
    genome_length = sum(iv.length() for iv in intervals)
    print(f"Total genome length: {genome_length}", file=sys.stderr)
    if genome_length <= 1:
        print("Genome length <= 1, exiting successfully.", file=sys.stderr)
        sys.exit(0)

    expected_length = genome_length // factor
    print(f"Expected cluster length: {expected_length}", file=sys.stderr)

    # Sort intervals by length (descending)
    intervals.sort(key=lambda iv: iv.length())

    clusters = []
    while intervals:
        cluster = Cluster()
        cluster_length = 0
        # Fill cluster until total_length >= expected_length or no intervals left
        while intervals and cluster_length < expected_length:
            iv = intervals.pop(0)
            iv_len = iv.length()
            remaining = expected_length - cluster_length
            if iv_len > remaining:
                # Split interval
                part1, part2 = iv.split(remaining)
                print(f"Splitting interval {iv} into {part1} (to fill cluster) and {part2} (to be processed)", file=sys.stderr)
                cluster.add(part1)
                cluster_length += part1.length()
                # Put the second part back at the front
                if part2:
                    intervals.insert(0, part2)
            else:
                cluster.add(iv)
                cluster_length += iv_len
        clusters.append(cluster)
        print(f"Cluster with {len(cluster.intervals)} intervals, total length {cluster.total_length()}", file=sys.stderr)

    # Save clusters
    for cluster in clusters:
        filename = cluster.save_as_bed(prefix)
        print(f"Cluster saved as: {filename}", file=sys.stderr)

if __name__ == "__main__":
    main()
