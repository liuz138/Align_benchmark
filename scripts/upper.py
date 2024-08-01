#!/usr/bin/env python3

def convert_to_uppercase(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                outfile.write(line)  # 保留FASTA标题行
            else:
                outfile.write(line.upper())  # 将序列行转换为大写

if __name__ == "__main__":
    input_file = 'hg38.chr12.fa'
    output_file = 'hg38.chr121.fa'
    convert_to_uppercase(input_file, output_file)
    print(f"Conversion complete. The output file is {output_file}")
