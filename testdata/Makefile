FETCH_SILVA = wget -O $@ https://zenodo.org/record/1172783/files/$@?download=1

all: gtdb silva unite

gtdb: GTDB_bac-arc_ssu_r86.fa.gz GTDB_dada2_assignment_species.fa.gz

silva: silva_nr_v132_train_set.fa.gz silva_species_assignment_v132.fa.gz

unite: sh_general_release_dynamic_all_02.02.2019.fasta.gz

sh_general_release_dynamic_all_02.02.2019.fasta.gz:
	wget https://files.plutof.ut.ee/public/orig/D6/96/D69658E99589D888A207805A744019DBA4EC0F603E67E53732767B3E03A5AA86.zip
	unzip D69658E99589D888A207805A744019DBA4EC0F603E67E53732767B3E03A5AA86.zip $(basename $@)
	gzip $(basename $@)
	rm D69658E99589D888A207805A744019DBA4EC0F603E67E53732767B3E03A5AA86.zip

silva_nr_v132_train_set.fa.gz:
	$(FETCH_SILVA)
	
silva_species_assignment_v132.fa.gz:
	$(FETCH_SILVA)

GTDB_bac-arc_ssu_r86.fa.gz:
	wget -O $@ https://zenodo.org/record/2541239/files/GTDB_bac-arc_ssu_r86.fa.gz?download=1

GTDB_dada2_assignment_species.fa.gz:
	wget -O $@ https://zenodo.org/record/2658728/files/GTDB_dada2_assignment_species.fa.gz?download=1
