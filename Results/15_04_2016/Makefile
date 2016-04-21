all: SDHC_FA.x SDHC_Divergence.x FA_1.dat div_1.dat # FA_2.dat div_1.dat div_2.dat

FA_1.dat:SDHC_FA.x
	./SDHC_FA.x ../Outputs/07_04_2016/snapshot_005.eigen_1 ../Outputs/07_04_2016/snapshot_005.eigen_2 ../Outputs/07_04_2016/snapshot_005.eigen_3 FA_1.dat

#FA_2.dat:SDHC_FA.x
#	./SDHC_FA.x ../Outputs/07_04_2016/snapshot_005.s2.00.eigenvec_1 ../Outputs/07_04_2016/snapshot_005.s2.00.eigenvec_2 ../Outputs/07_04_2016/snapshot_005.s2.00.eigenvec_3 FA_2.dat

div_1.dat:SDHC_Divergence.x
	./SDHC_Divergence.x ../Outputs/07_04_2016/snapshot_005.eigen_1 ../Outputs/07_04_2016/snapshot_005.eigen_2 ../Outputs/07_04_2016/snapshot_005.eigen_3 div_1.dat

#div_2.dat:SDHC_Divergence.x
#	./SDHC_Divergence.x ../Outputs/07_04_2016/snapshot_005.s2.00.eigenvec_1 ../Outputs/07_04_2016/snapshot_005.s2.00.eigenvec_2 ../Outputs/07_04_2016/snapshot_005.s2.00.eigenvec_3 div_2.dat

SDHC_FA.x:SDHC_read_eigenvalues_form_FA.c
	cc SDHC_read_eigenvalues_form_FA.c -lm -o SDHC_FA.x

SDHC_Divergence.x:SDHC_read_eigenvalues_form_Divergence.c
	cc SDHC_read_eigenvalues_form_Divergence.c -lm -o SDHC_Divergence.x
clean:
	rm *.dat *.pdf
