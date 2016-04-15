all: SDHC_FA.x test_1.dat test_2.dat

test_1.dat:SDHC_FA.x
	./SDHC_FA.x ../Outputs/07_04_2016/snapshot_005.s1.00.eigenvec_1 ../Outputs/07_04_2016/snapshot_005.s1.00.eigenvec_2 ../Outputs/07_04_2016/snapshot_005.s1.00.eigenvec_3 test_1.dat

test_2.dat:SDHC_FA.x
	./SDHC_FA.x ../Outputs/07_04_2016/snapshot_005.s2.00.eigenvec_1 ../Outputs/07_04_2016/snapshot_005.s2.00.eigenvec_2 ../Outputs/07_04_2016/snapshot_005.s2.00.eigenvec_3 test_2.dat

SDHC_FA.x:SDHC_read_eigenvalues_form_FA.c
	cc SDHC_read_eigenvalues_form_FA.c -lm -o SDHC_FA.x
clean:
	rm *.dat *.pdf
