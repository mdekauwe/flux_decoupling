USER = $(shell whoami)
MAKE = /usr/bin/make
DATA_DIR = data/processed
SRC = src
PLT_DIR = plotting_scripts
TEX_DIR = /Users/$(USER)/Dropbox/Decoupling_paper/figures
FIG_DIR = /Users/$(USER)/Dropbox/Decoupling_paper/figures/figs


all: data figs paper
data: $(DATA_DIR)/omega_fluxnet_PM.csv \
	  $(DATA_DIR)/omega_fluxnet_screened_PM.csv
figs: $(FIG_DIR)/Fluxnet_decoupling_boxplot.pdf \
		$(FIG_DIR)/omega_vs_sprecip.pdf \
		$(FIG_DIR)/omega_vs_wind.pdf \
		$(FIG_DIR)/omega_vs_lai.pdf \
		$(FIG_DIR)/omega_ENF.pdf
paper: $(TEX_DIR)/figures.pdf

##
# Make data files
##
$(DATA_DIR)/omega_fluxnet_PM.csv:	$(SRC)/estimate_decoupling_from_fluxnet_MPI.py
	python $<

$(DATA_DIR)/omega_fluxnet_screened_PM.csv.csv:	$(SRC)/screen_fluxnet_omega.py \
												$(DATA_DIR)/g1_fluxnet.csv
	python $<

##
# Figures - set these up so that if any data changes, they are *ALL* remade
##
$(FIG_DIR)/Fluxnet_decoupling_boxplot.pdf:	$(PLT_DIR)/plot_decoupling_boxplot.py \
											$(DATA_DIR)/*.csv
	python $<

$(FIG_DIR)/omega_ENF.pdf:	$(PLT_DIR)/plot_omega_ENF.py \
							$(DATA_DIR)/*.csv
	python $<

$(FIG_DIR)/omega_vs_sprecip..pdf:	$(PLT_DIR)/plot_omega_vs_3mth_rainfall.py \
								$(DATA_DIR)/*.csv
	python $<

$(FIG_DIR)/omega_vs_wind.pdf:	$(PLT_DIR)/plot_omega_vs_wind.py \
								$(DATA_DIR)/*.csv
	python $<

$(FIG_DIR)/omega_vs_lai.pdf:	$(PLT_DIR)/plot_omega_vs_lai.py \
								$(DATA_DIR)/*.csv
	python $<



# figure.pdf
##
$(TEX_DIR)/figures.pdf:	$(TEX_DIR)/figures.tex \
						$(FIG_DIR)/Fluxnet_decoupling_boxplot.pdf \
						$(FIG_DIR)/omega_vs_sprecip.pdf \
						$(FIG_DIR)/omega_vs_wind.pdf \
						$(FIG_DIR)/omega_vs_lai.pdf \
						$(FIG_DIR)/omega_ENF.pdf
	cd $(TEX_DIR) && $(MAKE) clean && $(MAKE)

cleanfigs:
	rm -f $(FIG_DIR)/*.pdf

.PHONY : clean

clean:
	rm -f $(DATA_DIR)/*.csv $(FIG_DIR)/*.pdf
