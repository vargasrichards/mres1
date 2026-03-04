# tests dealing with the homogeneous case, where we can compare to a simple branching process approximation. This is a sanity check for the package's `pr_fadeout` routine, which is used in the main analyses but is not trivial to verify by hand. We also have some checks that the homogeneous case is actually homogeneous (i.e. that the contact matrix is uniform and all classes are equally active).
using Test
using MixSwitchStochEpi

@testset "Homogeneous activity scores: Pr(ext) equal regardless of activity class and regardless of switching rate" begin
    hom_params = SEIRSStochParams(1.8, 1/3, 1/4, 0.0,
                             ones(5, 5) ./ 5, # M
                             MixSwitchStochEpi.uniform_switching(0.0, 5), # W
                             5,
                             fill(Int(50000 / 5), 5),
                             ones(5), 
                             50000)
    # check that extinction probabilities are equal across activity classes
    #calibrate once to get the right R0, then loop over initial class and check that Pr(extinction) is the same regardless of which class we seed in.
    calibrated_hom  = MixSwitchStochEpi.calibrate_parms(hom_params, 1.8)
    for class in 1:5
        init_state = MixSwitchStochEpi.make_initial_state(hom_params.class_sizes, 1; init_mode = :class, init_class = class)
        p_ext = MixSwitchStochEpi.pr_fadeout(hom_params, init_state; n_sims=1000, tmax=50.0, threshold_value=0.05)
        
    end

    # check that extinction probabilities are equal across switching rates
    for switch_rate in 0.1 .* [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        hom_params = SEIRSStochParams(1.8, 1/3, 1/4, 0.0,
                                 ones(5, 5) ./ 5, # M
                                 MixSwitchStochEpi.uniform_switching(switch_rate, 5), # W
                                 5,
                                 fill(Int(50000 / 5), 5),
                                 ones(5),
                                 50000)
        calibrated_hom  = MixSwitchStochEpi.calibrate_parms(hom_params, 1.8)
        for class in 1:5
            if class == 1
            end
            init_state = MixSwitchStochEpi.make_initial_state(hom_params.class_sizes, 1; init_mode = :class, init_class = class)
            p_ext = MixSwitchStochEpi.pr_fadeout(hom_params, init_state; n_sims=1000, tmax=50.0, threshold_value=0.05)
        end
    end
end


