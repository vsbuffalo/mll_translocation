# Makefile for running sets (or all) data, cleaning, etc.
# Vince Buffalo <vsbuffaloAAAAA@gmail.com> (with poly-A tail removed)


clean-bam:
	find data/ -type f -and -not -name  "*.sam" -and -not -wholename "*mll_template*" -exec rm {} \;
