import argparse,subprocess
import svgwrite

def parse_args():
    parser = argparse.ArgumentParser(description='Argument Parser for file paths and save paths')
    parser.add_argument('--ref', help='File path for reference file')
    parser.add_argument('--vec', help='File path for vector file')
    parser.add_argument('--out', help='File path for out file')
    parser.add_argument('--con',help='File path for consensus fasta file')
    return parser.parse_args()

def draw(ref, out,vec,con):
    input_file_path = ref
    vec_file_path=vec
    output_file_path = out
    consensus_path = con
    read_len=140
    #1 将不符合要求的行去除
    awk_command1 = f"awk 'length($4) >= {read_len}' {input_file_path} > {output_file_path}"
    subprocess.run(awk_command1, shell=True)

    #2 将参考序列写入最上面一行
    sequence = ""
    with open(vec_file_path, "r") as input_file:
        lines = input_file.readlines()
        if len(lines) >= 2:
            sequence = lines[1].strip()[0:140]
    with open(output_file_path, "r") as output_file:
        original_content = output_file.read()
    with open(output_file_path, "w") as output_file:
        output_file.write(sequence + "\n" + original_content)

    # 3 按照cigar值将每一行进行排序
    with open(output_file_path, "r") as f:
        lines = f.readlines()
    header = lines[0]
    sorted_lines = sorted(lines[1:], key=lambda x: int(x.split()[2].split('S')[0]), reverse=True)
    with open(output_file_path, "w") as f:
        f.write(header)
        f.writelines(sorted_lines)

        # 绘图1 将cigar中比对上和没比对上的结果区分开
    with open(output_file_path, 'r') as file:
        lines = file.readlines()
        # 准备一个列表来存储处理后的行
    processed_lines = [lines[0]]  # 保留第一行
    # 处理每一行
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != 4:
            # 如果行的格式不正确，就原样保留
            processed_lines.append(line)
            continue
        id_, start, cigar, sequence = parts
        match_length = int(cigar.split('S')[0])  # 提取并转换M前的数字
        modified_sequence = sequence[:match_length] + ' ' + sequence[match_length:]
        processed_lines.append(f"{id_} {start} {cigar} {modified_sequence}\n")
    with open(output_file_path, 'w') as file:
        file.writelines(processed_lines)

    # 读取文件内容
    with open(output_file_path, 'r') as file:
        lines = file.readlines()
    # 提取第一行并确定碱基序列的长度
    first_line = lines[0].strip()
    base_length = len(first_line.split('\t')[-1])
    # 处理剩余行
    processed_lines = [first_line]  # 包含第一行
    for line in lines[1:]:
        parts = line.strip().split('\t')
        sequence = parts[-1]
        space_index = sequence.find(' ')
        # 调整序列以与第一行对齐
        adjusted_sequence = ' ' * (base_length - space_index) + sequence
        parts[-1] = adjusted_sequence
        # 将调整后的行添加到结果中
        processed_lines.append('\t'.join(parts))
    # 将处理后的内容写回文件
    with open(output_file_path, 'w') as file:
        for line in processed_lines:
            file.write(line + '\n')
    with open(output_file_path,"r")as file:
        line=file.readlines()
        header=line[0]
        content=[line.split(None,3)[3]for line in lines[1:]]
    with open(output_file_path,"w")as file:
        file.write(header)
        for i in content:
            file.write(i)
    # 读取文件
    with open(output_file_path, 'r') as file:
        lines = file.readlines()
    # 提取第一行碱基
    first_line = lines[0].strip()
    # 计算第一行碱基的长度
    first_line_length = len(first_line)
    # 处理后的碱基列表
    processed_bases = []
    # 处理每一行的碱基片段
    for line in lines[1:]:
        # 分割每行的碱基和空格
        parts = line.strip().split(' ')
        # 获取空格的位置
        space_position = len(parts[0])
        # 计算需要移动的空格数
        shift_amount = first_line_length - space_position
        # 移动碱基片段并添加到处理后的列表
        processed_bases.append(' ' * shift_amount + parts[0] + ' ' + parts[1])
    # 将处理后的碱基写回文件
    with open(output_file_path, 'w') as file:
        base_count = len(first_line)
        spaces = ' ' * len(first_line)
    #添加consensus序列
    with open(consensus_path,'r') as file:
        line=file.readlines()
        consensus_line=''.join(line[1:])
        start_position = base_count - len(consensus_line)
        space = spaces[:start_position] + consensus_line + spaces[start_position + len(consensus_line):]
        space = space.rstrip(' \n')
    with open(output_file_path,'a') as file:
        file.write(f" {space} {first_line}\n")
        file.writelines('\n'.join(processed_bases))
        file.write("\n")
        file.write(spaces+' '+first_line)
    with open(output_file_path, 'r') as file:
        seq2 = file.read()

        dwg = svgwrite.Drawing("pic.svg", profile="full", width=500, height=500)
        x=0  # 将x的初始化移到循环外部
        y = 10
        counter = 0
        bin=2
        for i in seq2:
            if i == "\n":
                y += 4  # 换行时，y坐标增加以表示移动到新行
                counter = 0  # 换行后重置counter
                continue
            if i == " ":
                counter += 1  # 空格时，只增加counter
                continue
            x=counter*bin
            if i == "A":
                color = "red"
            elif i == "G":
                color = "blue"
            elif i == "C":
                color = "black"
            elif i == "T":
                color = "green"
            else:
                color = "yellow"

            dwg.add(dwg.text(i, insert=(x, y), fill=color, font_family="Courier New", font_size=3))
            counter += 1  # 每处理一个字符，counter增加
        dwg.save()

def run_draw(ref, out, vec,con):
    if ref and out:
        draw(ref, out, vec,con)

if __name__ == "__main__":
    args = parse_args()
    run_draw(args.ref, args.out, args.vec,args.con)