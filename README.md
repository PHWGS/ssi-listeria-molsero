# ssi-listeria-molsero
molecular serotyping script for Listeria monocytogenes by Kristoffer Kiil

Requirements: python3 and bwa in path

Usage: Molsero.py -p primers.fsa contigfile.fsa

Output: The detected bands and their respective sizes are output, and the predicted serotype is printed. If the pattern doesn't match a serotype, Eureka is printed instead.
