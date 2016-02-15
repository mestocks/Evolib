
Per chromosome stats:
```bash
zcat <file.vcf.gz> | vcf-nalleles | popstats > <out.bed>
```

Per chromosome, genes only:
```bash
zcat <file.vcf.gz> | vcf-nalleles | bedtools intersect -a stdin -b <annotation.genes.gff> | popstats > <out.bed>
```
