clean:
	rm -rf run_experiment experiments*

experiments:
	g++ -O2 run_experiment.cc oracles/*.cc oracles/util/*.cc -lboost_program_options -o run_experiment

results: experiments
	./run_all_experiments.py --debug --output_dir=experiments

figures: clean results
	./generate_figures.py --results_dir=experiments/results
