#bin/bash -l
module load bio
#metaphlan结果处理
# 把 S_前面的去掉,留下细菌物种丰度信息
grep -E "(s__)|(^ID)" merged_abundance_table.txt | grep -v "t__" | sed 's/^.*s__//g' > merged_abundance_table_species.txt

#将不同样品的丰度信息merge
merge_metaphlan_tables.py  *txt > merged_abundance_table.txt