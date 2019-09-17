FETCH_DECIPHER = wget -O $@ http://www2.decipher.codes/Classification/TrainingSets/$@
FETCH_SILVA    = wget -O $@ https://zenodo.org/record/1172783/files/$@?download=1

all: gtdb silva unite atacama_reads.qza

gtdb: GTDB_bac-arc_ssu_r86.fa.gz GTDB_dada2_assignment_species.fa.gz GTDB_r89-mod_June2019.RData

silva: silva_nr_v132_train_set.fa.gz silva_species_assignment_v132.fa.gz SILVA_SSU_r132_March2018.RData

unite: sh_general_release_dynamic_all_02.02.2019.fasta.gz UNITE_v2019_July2019.RData

sh_general_release_dynamic_all_02.02.2019.fasta.gz:
	wget https://files.plutof.ut.ee/public/orig/D6/96/D69658E99589D888A207805A744019DBA4EC0F603E67E53732767B3E03A5AA86.zip
	unzip D69658E99589D888A207805A744019DBA4EC0F603E67E53732767B3E03A5AA86.zip $(basename $@)
	gzip $(basename $@)
	rm D69658E99589D888A207805A744019DBA4EC0F603E67E53732767B3E03A5AA86.zip

UNITE_v2019_July2019.RData:
	$(FETCH_DECIPHER)

silva_nr_v132_train_set.fa.gz:
	$(FETCH_SILVA)
	
silva_species_assignment_v132.fa.gz:
	$(FETCH_SILVA)

SILVA_SSU_r132_March2018.RData:
	$(FETCH_DECIPHER)

GTDB_bac-arc_ssu_r86.fa.gz:
	wget -O $@ https://zenodo.org/record/2541239/files/GTDB_bac-arc_ssu_r86.fa.gz?download=1

GTDB_dada2_assignment_species.fa.gz:
	wget -O $@ https://zenodo.org/record/2658728/files/GTDB_dada2_assignment_species.fa.gz?download=1

GTDB_r89-mod_June2019.RData:
	$(FETCH_DECIPHER)

atacama_reads.qza:
	wget -O $@ https://docs.qiime2.org/2019.7/data/tutorials/atacama-soils/demux.qza
	@echo "Downloaded $@, you need to unpack it with 'qiime tools export --input-path $@ --output-path reads'."
	@echo "To make a selection of 17 large samples: 'find reads -name "*_R1_*" -size -100000c|xargs rm; for f in reads/*_R2_*; do if [ ! -e \$${f/_R2_/_R1_} ]; then rm \$$f; fi; done'"
