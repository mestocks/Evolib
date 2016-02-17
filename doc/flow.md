
Per chromosome stats:
```bash
zcat aln.vcf.gz | vcf-nalleles | popstats 10 > pstats.bed
```

Per window stats:
```bash
bedtools makewindows -g genome.fa.fai -w 10000 > genome_winsize.bed
zcat aln.vcf.gz | vcf-nalleles | bedtools intersect -wo -a stdin -b genome_winsize.bed | awk -F'\t' ' BEGIN { OFS="\t" }; { print $1,$2,$3,$7":"$8+1"-"$9,$5,$6 } ' | popstats 10 4 > pstats_10kb.bed
```

Per chromosome, genes only:
```bash
zcat <file.vcf.gz> | vcf-nalleles | bedtools intersect -a stdin -b <annotation.genes.gff> | popstats > <out.bed>
```
