import re

def parse_maf_and_compute_cigar(maf_file_path, output_cigar_path):
    sequences = []
    with open(maf_file_path, 'r') as file:
        for line in file:
            if line.startswith('s'):
                sequences.append(line.split()[6])  # 根据maf格式，序列在第7个元素

    with open(output_cigar_path, 'w') as f_cigar:
        for i in range(0, len(sequences), 2):
            if i + 1 < len(sequences):
                ref_seq = sequences[i]
                query_seq = sequences[i + 1]
                cigar = compute_cigar(ref_seq, query_seq)
                f_cigar.write(cigar + '\n')

def compute_cigar(ref_seq, query_seq):
    cigar = []
    M, I, D, X = 0, 0, 0, 0
    for i in range(len(ref_seq)):
        ref_char = ref_seq[i]
        query_char = query_seq[i]

        if ref_char != '-' and query_char != '-':
            if ref_char == query_char:
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

def flush_counts(cigar, M, X, I, D):
    if M > 0:
        cigar.append(f'{M}M')
    if X > 0:
        cigar.append(f'{X}X')
    if I > 0:
        cigar.append(f'{I}I')
    if D > 0:
        cigar.append(f'{D}D')

# 脚本使用
if __name__ == "__main__":
    maf_file = "simulation_0001.maf"  # 输入文件路径:maf文件
    output_cigar_file = "cigar.ans"  # 输出文件路径：答案cigar文件
    parse_maf_and_compute_cigar(maf_file, output_cigar_file)
