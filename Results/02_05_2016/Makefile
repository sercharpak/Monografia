all: pipeline_script_1.py pipeline_script_2.py pipeline_histograms.py ../FoF/src/grid_FA_0.6_Trace_1.0.dat ../FoF/src/grid_FA_0.8_Trace_1.0.dat ../FoF/src/grid_FA_0.8_Trace_0.0.dat ../FoF/src/grid_FA_0.9_Trace_0.0.dat ../FoF/src/grid_FA_1.0_Trace_0.0.dat group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.dat group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.dat group_results_seeds_FA_0.6_Trace_1.0_search_FA_1.0_Trace_0.0.dat group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.8_Trace_0.0.dat group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.9_Trace_0.0.dat group_results_seeds_FA_0.8_Trace_1.0_search_FA_1.0_Trace_0.0.dat inertia_cases_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.png

../FoF/src/grid_FA_0.6_Trace_1.0.dat:pipeline_script_1.py
	python pipeline_script_1.py ../Outputs/07_04_2016/snapshot_005.eigen_1 ../Outputs/07_04_2016/snapshot_005.eigen_2 ../Outputs/07_04_2016/snapshot_005.eigen_3 0.6 1.0

../FoF/src/grid_FA_0.8_Trace_1.0.dat:pipeline_script_1.py
	python pipeline_script_1.py ../Outputs/07_04_2016/snapshot_005.eigen_1 ../Outputs/07_04_2016/snapshot_005.eigen_2 ../Outputs/07_04_2016/snapshot_005.eigen_3 0.8 1.0

../FoF/src/grid_FA_0.8_Trace_0.0.dat:pipeline_script_1.py
	python pipeline_script_1.py ../Outputs/07_04_2016/snapshot_005.eigen_1 ../Outputs/07_04_2016/snapshot_005.eigen_2 ../Outputs/07_04_2016/snapshot_005.eigen_3 0.8 0.0

../FoF/src/grid_FA_0.9_Trace_0.0.dat:pipeline_script_1.py
	python pipeline_script_1.py ../Outputs/07_04_2016/snapshot_005.eigen_1 ../Outputs/07_04_2016/snapshot_005.eigen_2 ../Outputs/07_04_2016/snapshot_005.eigen_3 0.9 0.0

../FoF/src/grid_FA_1.0_Trace_0.0.dat:pipeline_script_1.py
	python pipeline_script_1.py ../Outputs/07_04_2016/snapshot_005.eigen_1 ../Outputs/07_04_2016/snapshot_005.eigen_2 ../Outputs/07_04_2016/snapshot_005.eigen_3 1.0 0.0

#129 seeds
group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.dat:pipeline_script_2.py
	python pipeline_script_2.py ../FoF/src/grid_FA_0.8_Trace_0.0.dat ../FoF/src/grid_FA_0.6_Trace_1.0.dat ../FoF/src/values_FA_0.8_Trace_0.0.dat

group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.dat:pipeline_script_2.py
	python pipeline_script_2.py ../FoF/src/grid_FA_0.9_Trace_0.0.dat ../FoF/src/grid_FA_0.6_Trace_1.0.dat ../FoF/src/values_FA_0.9_Trace_0.0.dat

#group_results_seeds_FA_0.6_Trace_1.0_search_FA_1.0_Trace_0.0.dat:pipeline_script_2.py
#	python pipeline_script_2.py ../FoF/src/grid_FA_1.0_Trace_0.0.dat ../FoF/src/grid_FA_0.6_Trace_1.0.dat ../FoF/src/values_FA_1.0_Trace_0.0.dat

#10^4 seeds
group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.8_Trace_0.0.dat:pipeline_script_2.py
	python pipeline_script_2.py ../FoF/src/grid_FA_0.8_Trace_0.0.dat ../FoF/src/grid_FA_0.8_Trace_1.0.dat ../FoF/src/values_FA_0.8_Trace_0.0.dat

group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.9_Trace_0.0.dat:pipeline_script_2.py
	python pipeline_script_2.py ../FoF/src/grid_FA_0.9_Trace_0.0.dat ../FoF/src/grid_FA_0.8_Trace_1.0.dat ../FoF/src/values_FA_0.9_Trace_0.0.dat

#group_results_seeds_FA_0.8_Trace_1.0_search_FA_1.0_Trace_0.0.dat:pipeline_script_2.py
#	python pipeline_script_2.py ../FoF/src/grid_FA_1.0_Trace_0.0.dat ../FoF/src/grid_FA_0.8_Trace_1.0.dat ../FoF/src/values_FA_1.0_Trace_0.0.dat

#histograms
#Seeds FA<0.6
inertia_cases_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.dat
volumes_distr_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.dat
inertia_diff_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.8_Trace_0.0.dat

inertia_cases_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.dat
volumes_distr_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.dat
inertia_diff_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.dat

#inertia_cases_0.6_Trace_1.0_search_FA_1.0_Trace_0.0.png:pipeline_histograms.py
#	python pipeline_histograms.py group_results_seeds_FA_0.6_Trace_1.0_search_FA_1.0_Trace_0.0.dat
#volumes_distr_0.6_Trace_1.0_search_FA_1.0_Trace_0.0.png:pipeline_histograms.py
#	python pipeline_histograms.py group_results_seeds_FA_0.6_Trace_1.0_search_FA_1.0_Trace_0.0.dat

#Seeds FA<0.8
inertia_cases_0.8_Trace_1.0_search_FA_0.8_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.8_Trace_0.0.dat
volumes_distr_0.8_Trace_1.0_search_FA_0.8_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.8_Trace_0.0.dat
inertia_diff_0.8_Trace_1.0_search_FA_0.8_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.8_Trace_0.0.dat

inertia_cases_0.8_Trace_1.0_search_FA_0.9_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.9_Trace_0.0.dat
volumes_distr_0.8_Trace_1.0_search_FA_0.9_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.9_Trace_0.0.dat
inertia_diff_0.8_Trace_1.0_search_FA_0.9_Trace_0.0.png:pipeline_histograms.py
	python pipeline_histograms.py group_results_seeds_FA_0.8_Trace_1.0_search_FA_0.9_Trace_0.0.dat

#inertia_cases_0.8_Trace_1.0_search_FA_1.0_Trace_0.0.png:pipeline_histograms.py
#	python pipeline_histograms.py group_results_seeds_FA_0.8_Trace_1.0_search_FA_1.0_Trace_0.0.dat
#volumes_distr_0.8_Trace_1.0_search_FA_1.0_Trace_0.0.png:pipeline_histograms.py
#	python pipeline_histograms.py group_results_seeds_FA_0.8_Trace_1.0_search_FA_1.0_Trace_0.0.dat

clean:
	rm *.dat *.png
