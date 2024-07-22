import re

def calculate_scores(cigar_string):
    """
    计算给定CIGAR字符串的四种得分：edit, gap-linear, gap-affine, gap-affine2p。

    参数:
    cigar_string (str): CIGAR字符串

    返回值:
    tuple: 四种得分 (edit_score, gap_linear_score, gap_affine_score, gap_affine2p_score)
    """
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

def calculate_score_and_write(cigar_file_path, score_file_path):
    """
    读取CIGAR字符串文件，计算每个字符串的得分，并将得分和原始CIGAR字符串写入文件。

    参数:
    cigar_file_path (str): 包含CIGAR字符串的输入文件路径
    score_file_path (str): 写入得分和CIGAR字符串的输出文件路径
    """
    with open(cigar_file_path, 'r') as cigar_file, open(score_file_path, 'w') as score_file:
        for line in cigar_file:
            cigar_string = line.strip()
            edit_score, gap_linear_score, gap_affine_score, gap_affine2p_score = calculate_scores(cigar_string)
            score_file.write(f"{edit_score} -{gap_linear_score} -{gap_affine_score} -{gap_affine2p_score} {cigar_string}\n")

# 脚本使用
if __name__ == "__main__":
    cigar_file = "cigar.ans"  # 输入文件路径
    score_file = "score.ans"  # 输出文件路径
    calculate_score_and_write(cigar_file, score_file)
