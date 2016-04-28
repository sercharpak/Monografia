all: pipeline_script_1.py ../FoF/src/grid_FA_0.6_Trace_1.0.dat ../FoF/src/grid_FA_0.8_Trace_0.0.dat

../FoF/src/grid_FA_0.6_Trace_1.0.dat:pipeline_script_1.py
	python pipeline_script_1.py ../Outputs/07_04_2016/snapshot_005.eigen_1 ../Outputs/07_04_2016/snapshot_005.eigen_2 ../Outputs/07_04_2016/snapshot_005.eigen_3 0.6 1.0

../FoF/src/grid_FA_0.8_Trace_0.0.dat:pipeline_script_1.py
	python pipeline_script_1.py ../Outputs/07_04_2016/snapshot_005.eigen_1 ../Outputs/07_04_2016/snapshot_005.eigen_2 ../Outputs/07_04_2016/snapshot_005.eigen_3 0.8 0.0
#To be re-written
pipeline_script_1.py:pipeline_script_2.py
	python pipeline_script_2.py ../FoF/src/grid_FA_0.8_Trace_0.0.dat ../FoF/src/grid_FA_0.6_Trace_1.0.dat


clean:
	rm *.dat *.pdf
