
Per chromosome stats:
```bash
zcat aln.vcf.gz | vcf-nalleles | popstats 10 > pstats.bed
```

Per window stats:
```bash
bedtools makewindows -g genome.fa.fai -w 10000 > genome_10kb.bed
zcat aln.vcf.gz | vcf-nalleles | bedtools intersect -wo -a stdin -b genome_10kb.bed | awk -F'\t' ' BEGIN { OFS="\t" }; { print $1,$2,$3,$7":"$8+1"-"$9,$5,$6 } ' | popstats 10 4 > pstats_10kb.bed
```

Per gene:
```bash
zcat aln.vcf.gz | vcf-nalleles | bedtools intersect -wo -a stdin -b ann.genes.gff | awk -F'\t' ' BEGIN { OFS="\t" }; { print $1,$2,$3,$15,$5,$6 } ' | popstats 10 4 > pstats_genes.bed
```

```bash
bblast gene.fa outgroup_genome.fa > outgroup_seq.bed
zcat aln.vcf.gz | bedtools intersect -wo -a stdin -b outgroup_seq.bed 
```
