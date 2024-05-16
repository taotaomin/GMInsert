import subprocess
import shutil
from datetime import datetime
import argparse
import os
import re

def parse_args():
    parser = argparse.ArgumentParser(description='Argument Parser for file paths and save paths')
    parser.add_argument('--ref', help='File path for reference file')
    parser.add_argument('--reads1', help='Reads1 file path')
    parser.add_argument('--reads2', help='Reads2 file path')
    return parser.parse_args()

def loop(ref,reads1,reads2):
    fasta_file_path = ref
    timestamp = datetime.now().strftime("%Y%m%d%H%M")
    with open(fasta_file_path, 'r') as fasta_file:
        # 跳过第一行（标题行）
        fasta_file.readline()

        # 读取第二行（DNA序列）
        dna_sequence = fasta_file.readline().strip()
        sequence_length = len(dna_sequence)

    subprocess.run(["bwa", "index", ref])
    subprocess.run(["bwa", "mem", "-t", "80", ref, reads1, reads2], stdout=open("mem_ref_GA.sam", "w"),
                   stderr=subprocess.PIPE)
    # 提取
    subprocess.run(["samtools", "view", "-bF", "4", "mem_ref_GA.sam"], stdout=open("map_ref_GA.bam", "w"),
                   stderr=subprocess.PIPE)
    new_filename4=f"loop_{timestamp}.bam"
    shutil.copyfile("map_ref_GA.bam",new_filename4)
    subprocess.run(["samtools", "view", "-h", "map_ref_GA.bam"], stdout=open("map_ref_GA.sam", "w"), stderr=subprocess.PIPE)
    awk_command = 'awk \'$1 !~ /^@/ && $6 ~ /^[0-9]+M[0-9]+S/\' map_ref_GA.sam'
    subprocess.run(awk_command, shell=True, stdout=open("filter.txt", "w"), stderr=subprocess.PIPE)
    # awk_command2 = r'''awk -F'\t' '{{ match($6, /[0-9]+M([0-9]+)S/, m); if ($4 + m[1] > {0}) print $0 }}' filter.txt > filter_1.txt'''.format(sequence_length)
    awk_command2 = r'''awk -F'\t' '{{ match($6, /([0-9]+)M[0-9]+S/, m); if ($4 + m[1] > {0}) print $0 }}' filter.txt > filter_1.txt'''.format(
        sequence_length)
    subprocess.run(awk_command2, shell=True, stdout=open("filter_1.txt", "w"), stderr=subprocess.PIPE)
    awk_command3 = r'''awk -F'\t' '{
        match($6, /^[0-9]+M/);
        split(substr($6, RSTART, RLENGTH - 1), arr, "M");
        $10 = substr($10, arr[1]+1);
        print $1 "\t" $4 "\t" $6 "\t" $10
    }' filter_1.txt > new.txt'''
    subprocess.run(awk_command3, shell=True)
    # sequence_length2 = 10  # 用实际值替换
    # awk_command4 = r'''awk -F'\t' '{{ print $1"\t"$2"\t"$3"\t"substr($4, 1, length($4)-{0}) }}' new.txt > new_1.txt'''.format(
    #     sequence_length2)
    # subprocess.run(awk_command4, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    awk_command5 = "awk '{print length($4)\"\\t\"$0}' new.txt"
    sort_command = "sort -nr"
    cut_command = "cut -f 2-"
    # 使用管道连接命令
    process1 = subprocess.Popen(awk_command5, shell=True, stdout=subprocess.PIPE)
    process2 = subprocess.Popen(sort_command, shell=True, stdin=process1.stdout, stdout=subprocess.PIPE)
    process1.stdout.close()  # 允许 process1 收到 SIGPIPE，如果 process2 退出
    process3 = subprocess.Popen(cut_command, shell=True, stdin=process2.stdout, stdout=open("new_1.txt", "w"))
    # 等待命令执行完成
    process3.communicate()

##############################################绘图文件获取#########################################################
    # 获取绘图文件 第一个命令
    command1 = "awk '{print $1, $4, $6, $10}' filter_1.txt | sort -k4,4n > filter_draw.txt"
    subprocess.run(command1, shell=True)

    new_filename_draw=f'draw_{timestamp}.txt'
    shutil.copyfile("filter_draw.txt",new_filename_draw)
#######################################################绘图文件获取完成#######################################################


    new_filename = f"new_{timestamp}.txt"
    shutil.copyfile("new_1.txt", new_filename)

    # 读取文件，获取碱基序列和对应的 IDs 和 cigars
    file_path = "new_1.txt"
    sequences = []
    ids_cigars = []

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                sequences.append(parts[3])
                ids_cigars.append((parts[1], parts[2]))
    with open(file_path, 'r+') as file:
        lines = file.readlines()  # 读取所有行

        # 将满足条件的行筛选出来
        lines_to_keep = [line for line in lines if line.split('\t')[3].strip()]  # 假设使用制表符分隔列

        # 定位文件指针到文件开头，清空文件内容
        file.seek(0)
        file.truncate()

        # 将满足条件的行写回原始文件
        file.writelines(lines_to_keep)

    # 寻找共识序列
    consensus = ''
    position = 0
    while True:
        base_count = {}
        base_max = 0
        base_consensus = None
        seqs_active = 0

        for seq in sequences:
            if position < len(seq):
                seqs_active += 1
                base = seq[position]
                base_count[base] = base_count.get(base, 0) + 1
                if base_count[base] > base_max:
                    base_max = base_count[base]
                    base_consensus = base

        if base_consensus is None or base_max <= seqs_active * 0.8:
            break

        consensus += base_consensus
        position += 1

    # 计算 mode 值
    counts = {}
    for (id, cigar) in ids_cigars:
        m = int(cigar.split('M')[0])
        id_m = int(id) + m
        counts[id_m] = counts.get(id_m, 0) + 1
    mode_val = max(counts, key=counts.get)

    # 将共识序列输出到 FASTA 文件
    output_file_path = "consensus.fasta"
    with open(output_file_path, 'w') as file:
        file.write(">gap_" + str(mode_val) + "-" + str(len(consensus)) + "-" + str(mode_val + len(consensus)) + "\n")
        file.write(consensus + "\n")
    # 读取 consensus.fasta 文件的第二行内容
    with open('consensus.fasta', 'r') as fasta_file:
        lines = fasta_file.readlines()
        second_line = lines[1].strip()

    new_filename2= f"consensus_{timestamp}.fasta"
    shutil.copyfile("consensus.fasta", new_filename2)

    txt_file_path = f"new_{timestamp}.txt"
    new_line = ["consensus_seq", "起始位置", "cigar", second_line]

    # 读取原始文件内容
    with open(txt_file_path, 'r') as original_file:
        original_content = original_file.readlines()

    # 在原始内容前添加新行，并写回文件
    with open(txt_file_path, 'w') as file_to_update:
        # 首先写入新行
        file_to_update.write('\t'.join(new_line) + '\n')

        # 接着写入原始内容
        file_to_update.writelines(original_content)

    # 读取原始参考文件
    with open(ref, 'r') as ref_handle:
        ref_lines = ref_handle.readlines()

    # 读取共识文件
    with open("consensus.fasta", 'r') as consensus_handle:
        consensus_lines = consensus_handle.readlines()
    # 从共识文件的第二行提取序列
    consensus_sequence = consensus_lines[1].strip()
    # 将共识序列追加到原始参考文件的第二行末尾
    ref_lines[1] = ref_lines[1].rstrip() + consensus_sequence + "\n"
    new_filename3 = f"combined_{timestamp}.fasta"
    with open(new_filename3, 'w') as new_file:
        new_file.writelines(ref_lines)
    # 将修改后的内容写回原始参考文件
    with open(ref, 'w') as ref_handle:
        ref_handle.writelines(ref_lines)

    # 源文件夹路径（假设在当前目录中）
    source_folder = os.path.abspath(os.path.dirname(ref))  # 获取源文件夹的绝对路径
    destination_folder = os.getcwd()  # 获取当前工作目录作为目标文件夹路径

    # 获取源文件夹下的所有文件和文件夹
    all_files = [f for f in os.listdir(source_folder) if os.path.isfile(os.path.join(source_folder, f))]

    # 使用正则表达式提取文件名中的时间戳
    timestamp_pattern = re.compile(r'\d{12}')  # 假设时间戳是12位数字

    # 存储已处理的时间戳及其对应的批次号
    processed_timestamps = {}

    # 对文件进行排序，确保按时间戳顺序处理
    all_files.sort()

    for file in all_files:
        timestamp_match = timestamp_pattern.search(file)
        if timestamp_match:
            timestamp = timestamp_match.group()

            # 检查这个时间戳是否已处理过
            if timestamp not in processed_timestamps:
                # 未处理过的时间戳，分配新的批次号
                loop_counter = len(processed_timestamps) + 1
                processed_timestamps[timestamp] = loop_counter
            else:
                # 已处理过的时间戳，使用现有批次号
                loop_counter = processed_timestamps[timestamp]

            # 根据时间戳和批次号创建目标文件夹
            destination_folder_name = f"{timestamp}_{loop_counter}"
            timestamp_folder = os.path.join(destination_folder, destination_folder_name)
            os.makedirs(timestamp_folder, exist_ok=True)

            # 移动文件到目标文件夹
            source_file = os.path.join(source_folder, file)
            destination_file = os.path.join(timestamp_folder, file)
            shutil.move(source_file, destination_file)

    # 获取所有目标文件夹的路径
    timestamp_folders = [f for f in os.listdir(destination_folder) if
                         os.path.isdir(os.path.join(destination_folder, f))]

    # 按照时间戳排序
    timestamp_folders.sort()

    # 重新命名文件夹
    for index, folder in enumerate(timestamp_folders, start=1):
        old_folder_path = os.path.join(destination_folder, folder)
        new_folder_name = f"{timestamp}_loop_{index}"
        new_folder_path = os.path.join(destination_folder, new_folder_name)
        os.rename(old_folder_path, new_folder_path)


def run_gap_attend(ref, reads1, reads2):
    if ref and reads1 and reads2:
            stop_condition = False
            while not stop_condition:
                # 调用 run_gap_attend 函数
                loop(ref, reads1, reads2)

                # 读取 consensus.fasta 文件，获取第二行的序列
                with open("consensus.fasta", 'r') as consensus_file:
                    consensus_lines = consensus_file.readlines()

                # 判断第二行序列是否为空
                if len(consensus_lines) >= 2 and len(consensus_lines[1].strip()) > 0:
                    # 如果不为空，继续循环
                    print("Continue loop.")
                else:
                    # 如果为空，停止循环
                    print("Stop loop.")
                    stop_condition = True

def main():
    args=parse_args()
    run_gap_attend(args.ref, args.reads1, args.reads2)

if __name__ == "__main__":
    main()
