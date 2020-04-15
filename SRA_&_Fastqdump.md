
#sra download from ncbi

```
while read id 
do
/home/xmixu/sratoolkit.2.9.6-1-centos_linux64/bin/prefetch.2.9.3 $id 
done < 2019_Moran_SRR_Acc_List.txt

```
#make the sra file into fastq file

```
while read id 
do
fastq-dump --outdir ./ --split-files /home/xmixu/ncbi/public/sra/$id.sra
done < /home/xmixu/infant_metadata/2019_Moran_cellhostmicrobes/2019_Moran_SRR_Acc_List.txt
```