import std.stdio;
import assimilation.EAKF;
import assimilation.likelihood.Likelihood;
import assimilation.likelihood.LikelihoodGetter;
import data.Vector;
import data.Ensemble;
import experiment.Analytics;
import experiment.Experiment;
import experiment.error.GaussianError;
import integrators.RK4;
import systems.System;
import systems.Lorenz63;

void main() {
	writeln("Experiment");
	RK4 rk4 = new RK4(new Lorenz63());
	Experiment process = new Experiment(rk4, new EAKF);
	process.getTruth(Vector(0, 0, 0), 0, 200, 0.01);
	process.setError(new GaussianError(Vector(0.1, 0.1, 0.1), process.truth, rk4));
	process.getObservations(0, 200, 0.1);
	process.setLikelihood(new LikelihoodGetter(process.observations, Vector(0.1, 0.1, 0.1)));
	process.getEnsembleTimeseries(0, 200, 0.01, 4, new Ensemble(Vector(0, 0, 0), 80, Vector(0.1, 0.1, 0.1)));
	writeln("Control RMSE is ", RMSE(process.ensembleSeries, process.truth));
}