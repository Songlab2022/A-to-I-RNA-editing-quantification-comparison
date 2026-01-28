import sys
if not '.gz' in sys.argv[1]:
    fastq_file = open(sys.argv[1], 'r')
else:
    import gzip
    fastq_file = gzip.open(sys.argv[1], 'r')
fastq_out = open(sys.argv[2], 'w')
length = int(sys.argv[3])

a = 1
for line in fastq_file:
    if '.gz' in sys.argv[1]: line = line.decode()
    if a % 4 in [1,3]:fastq_out.write(line)
    else:
        fastq_out.write(line.strip()[:length]+'\n')
    a += 1
fastq_out.close()
