bins = dna2codon popstats vcf-nalleles

HOME = $(shell echo $$HOME)
BASE = $(HOME)/.local/

##

bin = bin/
src = src/

BINS = $(addprefix $(bin),$(bins))

.PHONY:	all
all:	$(BINS)

$(bin)%:	$(src)%.c
	mkdir -p $(bin)
	gcc -I include/ -I $(BASE)include/librawk/ -L $(BASE)lib/ -o $@ $^ -lrawk -lm

.PHONY:	clean
clean:
	-rm $(BINS)
