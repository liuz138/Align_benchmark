import sys

# 检查是否提供了两个输入文件
if len(sys.argv) != 3:
    print("Usage: python compare_scores.py 1.ans 2.txt")
    sys.exit(1)

# 输入文件
file1 = sys.argv[1]
file2 = sys.argv[2]
output_file = "100L/20%/com-score_result.txt"

# 初始化计数器
total = 0
count_A = 0
count_B = 0
count_C = 0

# 读取文件内容
with open(file1, 'r') as f1, open(file2, 'r') as f2:
    lines1 = f1.readlines()
    lines2 = f2.readlines()

    # 确保两个文件的行数相同
    if len(lines1) != len(lines2):
        print("Error: The files do not have the same number of lines.")
        sys.exit(1)

    # 比较每一行的分数
    for line1, line2 in zip(lines1, lines2):
        score1 = abs(int(line1.split()[2]))
        score2 = abs(int(line2.split()[0]))
        total += 1
        if score1 == score2:
            count_A += 1
        #edit:>;gap_affine:< （abs:都是 >）
        elif score1 > score2:
            count_B += 1
        else:
            count_C += 1

# 计算占比
percent_A = (count_A / total) * 100
percent_B = (count_B / total) * 100
percent_C = (count_C / total) * 100

# 输出结果到文件
with open(output_file, 'a') as out_file:
    out_file.write(f"{file1} {file2}\n")
    out_file.write(f"总数：{total}条\n")
    out_file.write(f"A {count_A}条 {percent_A:.2f}%\n")
    out_file.write(f"B {count_B}条 {percent_B:.2f}%\n")
    out_file.write(f"C {count_C}条 {percent_C:.2f}%\n")
    out_file.write(f"{percent_A:.2f}%/{percent_B:.2f}%/{percent_C:.2f}%\n")
    out_file.write("----------------\n")

print(f"Results written to {output_file}")
