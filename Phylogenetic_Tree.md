```
#generate an anvi’o contigs database for each one 
for i in *fa
do
	anvi-script-FASTA-to-contigs-db $i
done
```


```
anvi-get-sequences-for-hmm-hits
--external-genomes external-genomes.txt 
-o concatenated-proteins.fa 
--gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 
--return-best-hit 
--get-aa-sequences 
--concatenate 
--just-do-it
```
```
#Visualize @local laptop
ssh -L 8080:localhost:8080 dylu@192.168.156.54
```

