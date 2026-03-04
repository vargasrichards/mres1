# basic scalar parameters are stable

    Code
      basic_scalar_parms()
    Output
      $beta
      [1] 0.1
      
      $sigma
      [1] 0.3333333
      
      $gamma
      [1] 0.25
      

# computed NGM and R0 are consistent

    Code
      compute_r0(pars)
    Output
      [1] 6.968962

---

    Code
      compute_ngm(pars)
    Output
                [,1]     [,2]     [,3]
      [1,] 0.2064516 0.064516 0.129032
      [2,] 0.0645160 2.645160 1.290320
      [3,] 0.1290324 1.290324 6.580648

# neaten_state output structure is consistent

    Code
      names(nt)
    Output
       [1] "class1"  "class2"  "S1"      "S2"      "E1"      "E2"      "I1"     
       [8] "I2"      "R1"      "R2"      "s_tot"   "e_tot"   "i_tot"   "r_tot"  
      [15] "pop_tot" "time"   

---

    Code
      head(nt)
    Output
         class1 class2       S1    S2       E1    E2       I1    I2        R1    R2
          <num>  <num>    <num> <num>    <num> <num>    <num> <num>     <num> <num>
      1:    500    500 495.0000   500 5.000000     0 0.000000     0 0.0000000     0
      2:    500    500 494.9317   500 3.643611     0 1.252249     0 0.1724088     0
      3:    500    500 494.7718   500 2.748288     0 1.903486     0 0.5764456     0
      4:    500    500 494.5668   500 2.144210     0 2.194588     0 1.0943916     0
      5:    500    500 494.3446   500 1.725565     0 2.273641     0 1.6562242     0
      6:    500    500 494.1212   500 1.426282     0 2.231335     0 2.2211949     0
            s_tot    e_tot    i_tot     r_tot pop_tot  time
            <num>    <num>    <num>     <num>   <num> <num>
      1: 995.0000 5.000000 0.000000 0.0000000    1000     0
      2: 994.9317 3.643611 1.252249 0.1724088    1000     1
      3: 994.7718 2.748288 1.903486 0.5764456    1000     2
      4: 994.5668 2.144210 2.194588 1.0943916    1000     3
      5: 994.3446 1.725565 2.273641 1.6562242    1000     4
      6: 994.1212 1.426282 2.231335 2.2211949    1000     5

