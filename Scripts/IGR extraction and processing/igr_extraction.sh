for seq in *.fas;
do;

cut -f 1,2 /Users/u5501917/bioinf_tools/ncbi-blast-2.14.0+/bin/MenY_db.fasta.fai > chrom.sizes;
blastn -db /Users/u5501917/bioinf_tools/ncbi-blast-2.14.0+/bin/MenYdb -query ${seq} -outfmt 6 -out coding_intermediates/${seq}.blast;

cut -f 2,9,10 coding_intermediates/${seq}.blast | awk '{if ($2>$3)print $1,$3,$2,".",".","-";else print $1,$2,$3,".",".","+";}' OFS='\t' |awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'|bedtools sort > coding_intermediates/${seq}.blast.bed;
bedtools getfasta -s -fi /Users/u5501917/bioinf_tools/ncbi-blast-2.14.0+/bin/MenY_db.fasta -bed coding_intermediates/${seq}.blast.bed -fo BLAST_out_coding/${seq}igr.fasta;
done
