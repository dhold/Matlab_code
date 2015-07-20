include ../make.config

latex2html=/Common/linux/encap/latex2html-2008/bin/latex2html

pdf: phi_ug.tex
	pdflatex phi_ug
	bibtex phi_ug
	pdflatex phi_ug
	pdflatex phi_ug

web:
	rm -rf ug/*
	$(latex2html) -local_icons -t "PHI $(PHIVERSION) User's Guide" -address "http://www.ks.uiuc.edu/Research/phi" -up_url "../" -up_title "Main PHI Page" -no_ps_images -image_type gif -discard -antialias -transparent -noaccent_images -dir ug -mkdir  phi_ug.tex
	tar czf ug.tar.gz ug

