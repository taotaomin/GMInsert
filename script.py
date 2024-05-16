import subprocess
import os

def run_script(reads1, reads2, tdna, ref):
    # 建立基因组文件的BWA索引
    subprocess.run(["bwa", "index", ref])
    # 建立TDNA文件的BWA索引
    subprocess.run(["bwa", "index", tdna])

    # 进行bwa比对
    subprocess.run(["bwa", "mem", "-t", "80", ref, reads1, reads2], stdout=open("mem_ref.sam", "w"),
                   stderr=subprocess.PIPE)
    # 提取
    subprocess.run(["samtools", "view", "-b", "-F", "8", "-f", "4", "mem_ref.sam"], stdout=open("map_ref.bam", "w"),
                   stderr=subprocess.PIPE)
    subprocess.run(["samtools", "fastq", "-n", "map_ref.bam"], stdout=open("map_ref.fastq", "w"),
                   stderr=subprocess.PIPE)

    # 与载体质粒比对
    subprocess.run(["bwa", "mem", tdna, "map_ref.fastq"], stdout=open("map_ref_mem_T.sam", "w"), stderr=subprocess.PIPE)
    subprocess.run(["samtools", "view", "-bF", "4", "map_ref_mem_T.sam"], stdout=open("map_ref_map_T.bam", "w"),
                   stderr=subprocess.PIPE)
    subprocess.run("samtools view map_ref_map_T.bam | cut -f1 > id.txt", shell=True)

    with open('id.txt', 'r') as f:
        id_list = [line.strip() for line in f]

    # 逐行读取SAM文件，并提取相关行
    with open('mem_ref.sam', 'r') as input_file, open('output.sam', 'w') as output_file:
        for line in input_file:
            if line.startswith('@'):
                # 处理SAM文件头部
                output_file.write(line)
            else:
                # 提取当前行的ID信息
                fields = line.split('\t')
                read_id = fields[0]

                # 检查ID是否存在于txt文件中
                if read_id in id_list:
                    # 保存相关行到输出文件中
                    output_file.write(line)

    subprocess.run(["samtools", "view", "-bS", "output.sam", "-o", "output.bam"], check=True)
    subprocess.run(["samtools", "sort", "output.bam", "-o", "output_sort.bam"], check=True)
    subprocess.run(["samtools", "index", "output_sort.bam"], check=True)

    ####################################start——end#####################################
    import pysam
    import re

    # 加载bam文件
    bam_file = pysam.AlignmentFile('output_sort.bam', 'rb')

    # 新建txt文件并打开
    with open('output.txt', 'w') as f:
        # 遍历bam文件的每一行
        for read in bam_file:
            # 如果质量值为0，则跳过
            if int(read.mapping_quality) == 0:
                continue

            # 提取所需信息
            read_id = read.query_name
            chromosome = read.reference_name
            start_pos = read.reference_start
            cigar_string = read.cigarstring

            # 提取CIGAR字符串中匹配上的数值
            match_length = sum([int(i[:-1]) for i in re.findall("\d+M", cigar_string)])

            # 计算新的结束位置
            end_pos = start_pos + match_length

            # 将信息写入txt文件中
            f.write(f"{read_id}\t{chromosome}\t{start_pos}\t{end_pos}\t{match_length}\n")

    # 关闭bam文件
    bam_file.close()

    ###########################取并集###############################################3
    def merge_intervals(intervals):
        intervals.sort(key=lambda x: x[0])
        merged = []
        for interval in intervals:
            if not merged or merged[-1][1] < interval[0]:
                merged.append(interval)
            else:
                merged[-1][1] = max(merged[-1][1], interval[1])
        return merged

    # 读取文件数据
    data = {}
    with open('output.txt', 'r') as f:
        for line in f:
            id, chromosome, start, end, gap = line.strip().split()
            interval = [int(start), int(end)]
            if chromosome in data:
                data[chromosome].append((interval, set([id]), int(gap)))
            else:
                data[chromosome] = [(interval, set([id]), int(gap))]

    # 对每个染色体编号的序列位置进行合并操作
    result = []
    for chromosome in data.keys():
        intervals = [x[0] for x in data[chromosome]]
        merged_intervals = merge_intervals(intervals)
        for interval in merged_intervals:
            seq_ids = set()
            seq_starts = []
            seq_ends = []
            for entry in data[chromosome]:
                if entry[0][0] <= interval[1] and entry[0][1] >= interval[0]:
                    seq_ids |= entry[1]
                    seq_starts.extend([str(entry[0][0])] * len(entry[1]))
                    seq_ends.extend([str(entry[0][1])] * len(entry[1]))
            result.append([chromosome, interval[0], interval[1], '|'.join(seq_ids),
                           '|'.join([f"{start}-{end}" for start, end in zip(seq_starts, seq_ends)])])

    # 将结果写入文件
    with open('bingji_new.txt', 'w') as f:
        for row in result:
            f.write('\t'.join(map(str, row)) + '\n')

    # 打开文件并读取内容
    with open('bingji_new.txt', 'r') as f:
        lines = f.readlines()

    # 添加新的一列
    chromosome_dict = {}
    new_lines = []
    current_chromosome = None
    current_index = 0

    for line in lines:
        elements = line.strip().split('\t')
        chromosome = elements[0]
        if chromosome != current_chromosome:
            current_chromosome = chromosome
            current_index = 1
        else:
            current_index += 1
        new_line = '\t'.join(
            [elements[0], f'C{current_index}', elements[1], elements[2], elements[3], elements[4]]) + '\n'
        new_lines.append(new_line)

    # 将修改后的内容写回文件
    with open('bingji_new.txt', 'w') as f:
        f.writelines(new_lines)

    # 打开文件并读取内容
    with open('bingji_new.txt', 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        elements = line.strip().split('\t')
        seq_ids = elements[4].split('|')
        subsets = elements[5].split('|')
        for i in range(len(seq_ids)):
            new_line = '\t'.join([elements[0], elements[1], elements[2], elements[3], seq_ids[i], subsets[i]])
            new_lines.append(new_line + '\n')

    # 将修改后的内容写回文件
    with open('bingji_new.txt', 'w') as f:
        f.writelines(new_lines)

    # 打开文件并读取内容
    with open('bingji_new.txt', 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        elements = line.strip().split('\t')
        positions = elements[-1].split('-')
        new_line = '\t'.join(elements[:-1] + positions)
        new_lines.append(new_line + '\n')

    # 将修改后的内容写回文件
    with open('bingji_new.txt', 'w') as f:
        f.writelines(new_lines)

    ########################################寻找gap#################################################
    def merge_intervals(intervals):
        intervals.sort(key=lambda x: x[0])
        merged = []
        for interval in intervals:
            if not merged or merged[-1][1] < interval[0]:
                merged.append(interval)
            else:
                merged[-1][1] = max(merged[-1][1], interval[1])
        return merged

    # 读取文件数据
    data = {}
    with open('output.txt', 'r') as f:
        for line in f:
            id, chromosome, start, end, gap = line.strip().split()
            interval = [int(start), int(end)]
            if chromosome in data:
                data[chromosome].append((interval, set([id]), int(gap)))
            else:
                data[chromosome] = [(interval, set([id]), int(gap))]

    # 对每个染色体编号的序列位置进行合并操作
    result = []
    for chromosome in data.keys():
        intervals = [x[0] for x in data[chromosome]]
        merged_intervals = merge_intervals(intervals)
        for interval in merged_intervals:
            seq_ids = set()
            seq_starts = []
            seq_ends = []
            for entry in data[chromosome]:
                if entry[0][0] <= interval[1] and entry[0][1] >= interval[0]:
                    seq_ids |= entry[1]
                    seq_starts.extend([str(entry[0][0])] * len(entry[1]))
                    seq_ends.extend([str(entry[0][1])] * len(entry[1]))
            result.append([chromosome, interval[0], interval[1], '|'.join(seq_ids),
                           '|'.join([f"{start}-{end}" for start, end in zip(seq_starts, seq_ends)])])

    # 将结果写入文件
    with open('bingji.txt', 'w') as f:
        for row in result:
            f.write('\t'.join(map(str, row)) + '\n')

    # 读取文件内容
    with open("bingji.txt", "r") as f:
        lines = f.readlines()

    # 添加新列并写入文件
    with open("bingji.txt", "w") as f:
        for i, line in enumerate(lines):
            columns = line.strip().split("\t")
            new_column_value = "C" + str(i + 1)
            new_line = "\t".join([columns[0], new_column_value, columns[1], columns[2], columns[3], columns[4]]) + "\n"
            f.write(new_line)

    # 打开文件
    with open('bingji.txt', 'r') as f:
        # 读取每一行数据
        lines = f.readlines()
        # 遍历每一行数据
        for i in range(len(lines)):
            line = lines[i].strip().split('\t')
            # 计算起始位置和终止位置之间的范围
            start = int(line[2])
            end = int(line[3])
            length = end - start + 1
            # 在第四列后面添加新列（范围）
            line.insert(4, str(length))
            # 将处理后的数据写回文件
            lines[i] = '\t'.join(line) + '\n'

    # 写入结果到文件
    with open('bingji.txt', 'w') as f:
        f.writelines(lines)
    # 打开bingji_new.txt文件
    with open('bingji.txt', 'r') as f:
        lines = f.readlines()

    # 初始化变量
    last_end = 0
    chromosome = ''
    results = []

    # 遍历每行数据
    for line in lines:
        # 按列分割数据
        data = line.rstrip().split('\t')
        # 获取染色体编号
        current_chromosome = data[0]
        # 如果当前行和上一行不在同一个染色体，则重置起始位置
        if current_chromosome != chromosome:
            last_end = 0
            chromosome = current_chromosome
        # 计算差值
        diff = int(data[2]) - last_end
        # 将结果添加到列表中
        results.append([current_chromosome, str(diff), data[2], data[3], data[2], str(last_end)])
        # 更新上一行终止位置
        last_end = int(data[3])

    # 将结果写入gap.txt文件
    with open('gap.txt', 'w') as f:
        for result in results:
            f.write('\t'.join(result) + '\n')

    # 打开gap.txt文件，并读取所有行数据
    with open('gap.txt', 'r') as f:
        lines = f.readlines()

    # 遍历每行数据
    for i in range(len(lines)):
        # 按列分割数据
        data = lines[i].rstrip().split('\t')
        # 将最后两列位置互换
        data[-2], data[-1] = data[-1], data[-2]
        # 将修改后的行数据重新赋值给原始列表
        lines[i] = '\t'.join(data) + '\n'

    # 将修改后的数据写回到gap.txt文件中
    with open('gap.txt', 'w') as f:
        f.writelines(lines)

    with open('gap.txt', 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        for line in lines:
            columns = line.split()
            new_line = '\t'.join([columns[0], columns[1], columns[4], columns[5]])
            f.write(new_line + '\n')
        f.truncate()

    with open('gap.txt', 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        for line in lines:
            columns = line.split()
            if columns[2] != '0':
                f.write(line)
        f.truncate()

    with open('gap.txt', 'r+') as f:
        lines = f.readlines()
        f.seek(0)

        # 添加表头，并根据列宽对齐
        f.write('{:<10} {:<10} {:<10} {:<10}\n'.format("ref_ID", "gap_len", "gap_S", "gap_E"))

        for line in lines:
            f.write(line)
        f.truncate()
