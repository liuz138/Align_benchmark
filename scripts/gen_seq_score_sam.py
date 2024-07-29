from Bio import SeqIO
from Bio.Seq import Seq
import re

# 输入输出文件路径
fasta_file = "500L/20%/500L-0.2_0001.ref"
fastq_file = '500L/20%/500L-0.2_0001.fastq'
maf_file = '500L/20%/500L-0.2_0001.maf'
seq_output_file = "500L/20%/500L-0.2-sam.seq"
score_output_file = "500L/20%/score-sam.ans"

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

# 生成双序列文件1.seq
with open(seq_output_file, "w") as output_file:
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

# 生成CIGAR
def parse_maf_and_compute_cigar(maf_file_path):
    sequences = []
    with open(maf_file_path, 'r') as file:
        for line in file:
            if line.startswith('s'):
                sequences.append(line.split()[6])  # 根据maf格式，序列在第7个元素

    cigars = []
    sam_cigars = []  # 新增：用于存储SAM格式的CIGAR字符串
    for i in range(0, len(sequences), 2):
        if i + 1 < len(sequences):
            ref_seq = sequences[i]
            query_seq = sequences[i + 1]
            cigar = compute_cigar(ref_seq, query_seq)
            sam_cigar = compute_sam_cigar(ref_seq, query_seq)  # 新增：生成SAM格式的CIGAR字符串
            cigars.append((cigar, sam_cigar))  # 修改：同时保存标准CIGAR和SAM格式CIGAR
    return cigars

def compute_cigar(ref_seq, query_seq):
    cigar = []
    M, I, D, X = 0, 0, 0, 0
    for i in range(len(ref_seq)):
        ref_char = ref_seq[i]
        query_char = query_seq[i]

        if ref_char != '-' and query_char != '-':
            if ref_char == query_char or ref_char == 'N' or query_char == 'N':
                if X > 0 or I > 0 or D > 0:
                    flush_counts(cigar, M, X, I, D)
                    M, X, I, D = 0, 0, 0, 0
                M += 1
            else:
                if M > 0 or I > 0 or D > 0:
                    flush_counts(cigar, M, X, I, D)
                    M, X, I, D = 0, 0, 0, 0
                X += 1
        elif ref_char == '-' and query_char != '-':
            if M > 0 or X > 0 or D > 0:
                flush_counts(cigar, M, X, I, D)
                M, X, I, D = 0, 0, 0, 0
            I += 1
        elif query_char == '-' and ref_char != '-':
            if M > 0 or X > 0 or I > 0:
                flush_counts(cigar, M, X, I, D)
                M, X, I, D = 0, 0, 0, 0
            D += 1

    flush_counts(cigar, M, X, I, D)
    return ''.join(cigar)

# 新增：计算SAM格式CIGAR字符串的函数
def compute_sam_cigar(ref_seq, query_seq):
    sam_cigar = []
    M, I, D = 0, 0, 0
    for i in range(len(ref_seq)):
        ref_char = ref_seq[i]
        query_char = query_seq[i]

        if ref_char != '-' and query_char != '-':
            if I > 0 or D > 0:
                flush_sam_counts(sam_cigar, M, I, D)
                M, I, D = 0, 0, 0
            M += 1
        elif ref_char == '-' and query_char != '-':
            if M > 0 or D > 0:
                flush_sam_counts(sam_cigar, M, I, D)
                M, I, D = 0, 0, 0
            I += 1
        elif query_char == '-' and ref_char != '-':
            if M > 0 or I > 0 :
                flush_sam_counts(sam_cigar, M, I, D)
                M, I, D = 0, 0, 0
            D += 1

    flush_sam_counts(sam_cigar, M, I, D)
    return ''.join(sam_cigar)

# 修改：新增用于刷新SAM格式CIGAR计数的函数
def flush_counts(cigar, M, X, I, D):
    if M > 0:
        cigar.append(f'{M}M')
    if X > 0:
        cigar.append(f'{X}X')
    if I > 0:
        cigar.append(f'{I}I')
    if D > 0:
        cigar.append(f'{D}D')

# 新增：刷新SAM格式CIGAR计数的函数
def flush_sam_counts(sam_cigar, M, I, D):
    if M > 0:
        sam_cigar.append(f'{M}M')
    if I > 0:
        sam_cigar.append(f'{I}I')
    if D > 0:
        sam_cigar.append(f'{D}D')

# 计算得分
def calculate_scores(cigar_string):
    edit_score = 0
    gap_linear_score = 0
    gap_affine_score = 0
    gap_affine2p_score = 0

    matches = re.findall(r'(\d+)([MIDX])', cigar_string)

    for match in matches:
        count, operation = match
        count = int(count)

        if operation != 'M':
            edit_score += count

        if operation == 'M':
            pass
        elif operation == 'X':
            gap_linear_score += 4 * count
            gap_affine_score += 4 * count
            gap_affine2p_score += 4 * count
        elif operation == 'I' or operation == 'D':
            gap_linear_score += 2 * count

            if count == 1:
                gap_affine_score += 8
                gap_affine2p_score += 8
            else:
                gap_affine_score += 8 + 2 * (count - 1)
                if count <= 10:
                    gap_affine2p_score += 8 + 2 * (count - 1)
                else:
                    gap_affine2p_score += 46 + 2 * (count - 10)

    return edit_score, gap_linear_score, gap_affine_score, gap_affine2p_score

# 修改：计算得分并写入文件，使用SAM格式的CIGAR字符串
def calculate_score_and_write(cigar_tuples, score_file_path):
    with open(score_file_path, 'w') as score_file:
        for cigar_string, sam_cigar in cigar_tuples:
            edit_score, gap_linear_score, gap_affine_score, gap_affine2p_score = calculate_scores(cigar_string)
            score_file.write(f"{edit_score} -{gap_linear_score} -{gap_affine_score} -{gap_affine2p_score} {sam_cigar}\n")

# 解析MAF文件并生成CIGAR字符串
cigar_tuples = parse_maf_and_compute_cigar(maf_file)

# 计算得分并写入文件
calculate_score_and_write(cigar_tuples, score_output_file)
