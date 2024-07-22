from Bio import SeqIO
from Bio.Seq import Seq

# 参考序列文件路径
fasta_file = "./simulation.fasta"  
fastq_file = './simulation_0001.fastq'
maf_file = './simulation_0001.maf'

# 记录序列和位置信息
alignments = {}
reverse_complement_flags = {}

# 读取MAF文件，记录序列名和在参考上的位置，并记录是否需要反向互补
with open(maf_file, 'r') as file:
    for line in file:
        if line.startswith('a'):
            ref_start = None
            ref_end = None
            segment_name = None
        elif line.startswith('s'):
            fields = line.split()
            if fields[1].startswith('ref'):
                ref_start = int(fields[2])
                ref_end = int(fields[2]) + int(fields[3])
            else:
                segment_name = fields[1]
                reverse_complement = fields[4] == '-'
            
            if ref_start is not None and segment_name is not None:
                alignments[segment_name] = (ref_start, ref_end)
                reverse_complement_flags[segment_name] = reverse_complement
                ref_start = None
                segment_name = None

# 读取FASTQ文件
fastq_sequences = []
with open(fastq_file, "r") as file:
    for record in SeqIO.parse(file, "fastq"):
        fastq_sequences.append(str(record.seq))

# 将结果写入11111.seq文件
with open("1.seq", "w") as output_file:
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            reference_sequence = str(record.seq)
            for segment_name, (ref_start, ref_end) in alignments.items():
                if len(fastq_sequences) > 0:
                    ref_segment = reference_sequence[ref_start:ref_end]
                    
                    fastq_seq = fastq_sequences.pop(0)
                    
                    # 如果需要反向互补
                    if reverse_complement_flags[segment_name]:
                        fastq_seq = str(Seq(fastq_seq).reverse_complement())
                    
                    output_file.write(f"{ref_segment}\n")
                    output_file.write(f"{fastq_seq}\n")
