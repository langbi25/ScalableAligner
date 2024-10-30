#!/bin/bash

# 定义目标文件夹
directory="/mnt/data/22_liangjialang/dataset/chroms"

# 遍历目录下所有文件
for file in "$directory"/*; do
  if [ -f "$file" ]; then
    echo "Processing file: $file"
    # 执行命令，"$file" 作为参数传递给命令
    # your_command "$file"
    filename=$(basename "$file")
     ./bin/bwt_index $file /mnt/data/22_liangjialang/dataset/chroms_index1/$filename
  fi
done
