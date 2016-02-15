HOME = $(shell echo $$HOME)
BASE = $(addprefix $(HOME)/,.local/)

bin = bin/
src = src/

#BIN = $(addprefix $(HOME),$(bin))
#SRC = $(addprefix $(HOME),$(src))

.PHONY:	all
all:	$(bin)popstats $(bin)vcf-nalleles

$(bin)%:	$(src)%.c
	gcc -I include/ -I $(BASE)include/librawk/ -L $(BASE)lib/ -o $@ $^ -lrawk -lm
