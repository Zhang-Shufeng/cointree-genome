#!/bin/bash

# ======================= CONFIG ===========================
INPUT_LIST="sample_list.txt"  # 你的样本名列表，每行一个，例如：DZ1、DZ2 等
DATA_DIR="../01_data"
BLOCK_DIR="../03_blocks"
OUTPUT_DIR="../05_kaks_resoult"

THREADS=10
ALIGN_METHOD="mafft"
KAKS_METHOD="YN"

# ======================= START PIPELINE ====================
set -e  # 出错即终止

# 检查依赖软件是否存在（可选）
command -v ParaAT.pl >/dev/null 2>&1 || { echo "ParaAT.pl not found. Please check your environment."; exit 1; }
command -v KaKs_Calculator >/dev/null 2>&1 || { echo "KaKs_Calculator not found."; exit 1; }
command -v seqkit >/dev/null 2>&1 || { echo "seqkit not found."; exit 1; }

# 创建输出目录（如果不存在）
mkdir -p "${OUTPUT_DIR}"

# 读取样本列表并逐个处理
while read -r sample; do
    echo "==> Processing: ${sample}"

    # 拷贝 CDS 文件
    cp "${DATA_DIR}/${sample}_cds.fa" ./

    # 清理注释中的 locus 字段
    sed '/^>/s/\(>[^locus]*\)locus=.*\s*/\1/' "${sample}_cds.fa" > "${sample}_cds.fas"
    mv "${sample}_cds.fas" "${sample}_cds.fa"

    # 翻译为蛋白
    seqkit translate --trim "${sample}_cds.fa" > "${sample}_pro.fa"

    # 保留标题行（不含序列）
    awk '{print $1}' "${sample}_cds.fa" > "${sample}_new.cds.fa"
    awk '{print $1}' "${sample}_pro.fa" > "${sample}_new.pro.fa"

    # 替换标题中的点号
    sed -i "s/\./_/g" "${sample}_new.cds.fa"
    sed -i "s/\./_/g" "${sample}_new.pro.fa"

    # 拷贝 blocks 文件并清洗
    cp "${BLOCK_DIR}/${sample}.i1_resoult.blocks" ./
    sed -i "s/\./_/g" "${sample}.i1_resoult.blocks"

    # 写入线程数
    echo "${THREADS}" > proc.txt

    # Run ParaAT
    ParaAT.pl -h "${sample}.i1_resoult.blocks" \
              -a "${sample}_new.pro.fa" \
              -n "${sample}_new.cds.fa" \
              -p proc.txt -o "${sample}_align_out" \
              -m "${ALIGN_METHOD}" -f axt

    # 合并 axt 文件
    cat "${sample}_align_out"/*.axt > "${sample}_merge_align.axt"

    # 计算 KaKs
    KaKs_Calculator -m "${KAKS_METHOD}" \
                    -i "${sample}_merge_align.axt" \
                    -o "${sample}_resoult.txt"

    # 移动结果
    mv "${sample}_resoult.txt" "${OUTPUT_DIR}/"

    # 清理中间文件
    rm -rf "${sample}"*

    echo "==> ${sample} done."

done < "${INPUT_LIST}"

echo "==> All samples processed!"
