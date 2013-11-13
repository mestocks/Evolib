def multi_row_iter(*args):
        """
        Allows simultaneous looping through multiple VCF files.
        
        multi_row_iter(vcf1, vcf2, ..., ordered_chrs)
        
        The VCF classes are given in the first arguments with the last
        argument reserved for a list of all chromosomes in the order in which they
        appear in the VCF files. All chromosomes do not have to appear in each of
        the VCF files themselves but the list must be in the order in which they occur.
        
        Usage:
        
            from SequenceFormats import VCF
            
            vcfI = VariantCallFormat('Aedit.vcf')
            vcfII = VariantCallFormat('BCedit.vcf')
            vcfIII = VariantCallFormat('Dedit.vcf')
            
            for I, II, III in gnrtr2(vI, vII, vIII, ref_chrs):
            
                 if I is not None:
                     print I['CHROM'], I['POS'],
                 else:
                     print 'None', 'None',
                     
                 if II is not None:
                     print II['CHROM'], II['POS'],
                 else:
                     print 'None', 'None',
                     
                 if III is not None:
                     print III['CHROM'], III['POS'],
                 else:
                     print 'None', 'None',
                     
             print '\n',
        """
        iterables = [it.row_iter() for it in args[: -1]]
        chrs_names = args[-1]
    
        first_time = True
        chr_index = 0
        while chr_index < 10:
            
            if first_time is True:
                rows = [next(row) for row in iterables]
                first_time = False
            
            current_chr = chrs_names[chr_index]
            
            chrs_list = [str(i['CHROM']) for i in rows]
            
            if current_chr not in chrs_list:
                chr_index += 1
            else:
                positions = [i['POS'] for i in rows]
                min_pos = min([rows[i]['POS'] for i in range(len(rows)) if str(rows[i]['CHROM']) == current_chr])
                
                rows_to_yield = []
                for index in range(len(rows)):
                    if chrs_list[index] == current_chr and positions[index] == min_pos:
                        rows_to_yield.append(rows[index])
                        rows[index] = next(iterables[index])
                    else:
                        rows_to_yield.append(None)
                
                yield rows_to_yield