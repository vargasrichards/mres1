
using MixSwitchStochEpi
using Random
using QuadGK
using CSV

# Build a default parametrisation and calibrate β to target R0 values.
params, init, tspan = MixSwitchStochEpi.default_parametrisation(n_activity=5,
                                                                 Ntot=20000,
                                                                 ξ = 0.,
                                                                 ε = 0.,
                                                                 gamma_mean=5.0,
                                                                 gamma_variance=5.0,
                                                                 p_inf=1,
                                                                 init_mode=:class,
                                                                 init_class=1)

pcal  = MixSwitchStochEpi.calibrate_parms(params, 1.8)
pcal2 = MixSwitchStochEpi.calibrate_parms(params, 1.2)
pcal3 = MixSwitchStochEpi.calibrate_parms(params, 0.8)

pswitch1 = MixSwitchStochEpi.SEIRSStochParams(pcal.β, pcal.σ, pcal.γ, pcal.ω, pcal.M,
               MixSwitchStochEpi.uniform_switching(0.5, pcal.n_activity),
               pcal.n_activity, pcal.class_sizes, pcal.act_levels, pcal.pop_size)
pswitchcal1 = MixSwitchStochEpi.calibrate_parms(pswitch1, 1.8)
pswitchcal2 = MixSwitchStochEpi.calibrate_parms(pswitch1, 1.2)
pswitchcal3 = MixSwitchStochEpi.calibrate_parms(pswitch1, 0.8)

a = plot_pext_initial(pcal;       n_sims=10000, init_class=1, threshold_fraction=0.01)
b = plot_pext_initial(pcal;       n_sims=10000, init_class=5, threshold_fraction=0.01)
d = plot_pext_initial(pcal2;      n_sims=10000, init_class=1, threshold_fraction=0.01)
e = plot_pext_initial(pcal2;      n_sims=10000, init_class=5, threshold_fraction=0.01)
f = plot_pext_initial(pcal3;      n_sims=10000, init_class=1, threshold_fraction=0.01)
g = plot_pext_initial(pcal3;      n_sims=10000, init_class=5, threshold_fraction=0.01)
h = plot_pext_initial(pswitchcal1; n_sims=10000, init_class=1, threshold_fraction=0.01)
i = plot_pext_initial(pswitchcal1; n_sims=10000, init_class=5, threshold_fraction=0.01)
j = plot_pext_initial(pswitchcal2; n_sims=10000, init_class=1, threshold_fraction=0.01)
k = plot_pext_initial(pswitchcal2; n_sims=10000, init_class=5, threshold_fraction=0.01)
l = plot_pext_initial(pswitchcal3; n_sims=10000, init_class=1, threshold_fraction=0.01)
m = plot_pext_initial(pswitchcal3; n_sims=10000, init_class=5, threshold_fraction=0.01)

combined = vcat(a, b, d, e, f, g, h, i, j, k, l, m)
CSV.write("extinction_threshold_all.csv", combined)
@info "Finished: results written to extinction_threshold_all.csv"
