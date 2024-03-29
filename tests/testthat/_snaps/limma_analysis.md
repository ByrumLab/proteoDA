# fit_limma_model gives consistent output

    Code
      fit_limma_model(a)
    Output
      $data
                 sample1  sample2  sample3  sample4  sample5  sample6
      protein1  0.000000 3.459432 4.392317 4.954196 5.357552 5.672425
      protein2  1.000000 3.584963 4.459432 5.000000 5.392317 5.700440
      protein3  1.584963 3.700440 4.523562 5.044394 5.426265 5.727920
      protein4  2.000000 3.807355 4.584963 5.087463 5.459432 5.754888
      protein5  2.321928 3.906891 4.643856 5.129283 5.491853 5.781360
      protein6  2.584963 4.000000 4.700440 5.169925 5.523562 5.807355
      protein7  2.807355 4.087463 4.754888 5.209453 5.554589 5.832890
      protein8  3.000000 4.169925 4.807355 5.247928 5.584963 5.857981
      protein9  3.169925 4.247928 4.857981 5.285402 5.614710 5.882643
      protein10 3.321928 4.321928 4.906891 5.321928 5.643856 5.906891
      
      $annotation
                uniprot_id
      protein1    protein1
      protein2    protein2
      protein3    protein3
      protein4    protein4
      protein5    protein5
      protein6    protein6
      protein7    protein7
      protein8    protein8
      protein9    protein9
      protein10  protein10
      
      $metadata
              sample_ID group sex treatment
      sample1   sample1     A   M   control
      sample2   sample2     B   F   control
      sample3   sample3     A   M   control
      sample4   sample4     B   F treatment
      sample5   sample5     A   M treatment
      sample6   sample6     B   F treatment
      
      $design
      $design$design_formula
      [1] "~0 + group"
      
      $design$design_matrix
              A B
      sample1 1 0
      sample2 0 1
      sample3 1 0
      sample4 0 1
      sample5 1 0
      sample6 0 1
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      
      $eBayes_fit
      An object of class "MArrayLM"
      $coefficients
                       A        B
      protein1  3.249956 4.695351
      protein2  3.617250 4.761801
      protein3  3.844930 4.824251
      protein4  4.014798 4.883235
      protein5  4.152546 4.939178
      protein6  4.269655 4.992427
      protein7  4.372277 5.043269
      protein8  4.464106 5.091945
      protein9  4.547539 5.138658
      protein10 4.624225 5.183582
      
      $rank
      [1] 2
      
      $assign
      [1] 1 1
      
      $qr
      $qr
                       A          B
      sample1 -1.7320508  0.0000000
      sample2  0.0000000 -1.7320508
      sample3  0.5773503  0.0000000
      sample4  0.0000000  0.5773503
      sample5  0.5773503  0.0000000
      sample6  0.0000000  0.5773503
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      $qraux
      [1] 1.57735 1.57735
      
      $pivot
      [1] 1 2
      
      $tol
      [1] 1e-07
      
      $rank
      [1] 2
      
      
      $df.residual
       [1] 4 4 4 4 4 4 4 4 4 4
      
      $sigma
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       2.171309  1.805053  1.596609  1.451764  1.341391  1.252654  1.178762  1.115685 
       protein9 protein10 
       1.060837  1.012458 
      
      $cov.coefficients
                A         B
      A 0.3333333 0.0000000
      B 0.0000000 0.3333333
      
      $stdev.unscaled
                        A         B
      protein1  0.5773503 0.5773503
      protein2  0.5773503 0.5773503
      protein3  0.5773503 0.5773503
      protein4  0.5773503 0.5773503
      protein5  0.5773503 0.5773503
      protein6  0.5773503 0.5773503
      protein7  0.5773503 0.5773503
      protein8  0.5773503 0.5773503
      protein9  0.5773503 0.5773503
      protein10 0.5773503 0.5773503
      
      $pivot
      [1] 1 2
      
      $Amean
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       3.972654  4.189525  4.334591  4.449017  4.545862  4.631041  4.707773  4.778025 
       protein9 protein10 
       4.843098  4.903904 
      
      $method
      [1] "ls"
      
      $design
              A B
      sample1 1 0
      sample2 0 1
      sample3 1 0
      sample4 0 1
      sample5 1 0
      sample6 0 1
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      $df.prior
       [1] Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf
      
      $s2.prior
      [1] 2.347193
      
      $var.prior
      [1] 6.816652 6.816652
      
      $proportion
      [1] 0.01
      
      $s2.post
       [1] 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193
       [9] 2.347193 2.347193
      
      $t
                       A        B
      protein1  3.674208 5.308285
      protein2  4.089448 5.383409
      protein3  4.346849 5.454012
      protein4  4.538892 5.520696
      protein5  4.694622 5.583941
      protein6  4.827018 5.644141
      protein7  4.943037 5.701620
      protein8  5.046853 5.756650
      protein9  5.141177 5.809461
      protein10 5.227874 5.860250
      
      $df.total
       [1] 40 40 40 40 40 40 40 40 40 40
      
      $p.value
                           A            B
      protein1  6.993473e-04 4.413049e-06
      protein2  2.029093e-04 3.465792e-06
      protein3  9.224506e-05 2.760859e-06
      protein4  5.078967e-05 2.226698e-06
      protein5  3.116181e-05 1.815553e-06
      protein6  2.051361e-05 1.494720e-06
      protein7  1.419448e-05 1.241299e-06
      protein8  1.019662e-05 1.038942e-06
      protein9  7.542600e-06 8.757757e-07
      protein10 5.713087e-06 7.430412e-07
      
      $lods
                        A         B
      protein1  0.3072394  7.304137
      protein2  1.8439827  7.687015
      protein3  2.8791225  8.051757
      protein4  3.6925675  8.400614
      protein5  4.3780142  8.735402
      protein6  4.9789443  9.057611
      protein7  5.5192770  9.368481
      protein8  6.0136575  9.669057
      protein9  6.4717459  9.960228
      protein10 6.9002739 10.242760
      
      $F
       [1] 20.83885 22.85234 24.32068 25.53981 26.60994 27.57822 28.47104 29.30487
       [9] 30.09077 30.83660
      
      $F.p.value
       [1] 8.908471e-10 1.189474e-10 2.739460e-11 8.094687e-12 2.776202e-12
       [6] 1.054222e-12 4.317003e-13 1.875231e-13 8.545621e-14 4.053532e-14
      
      
      $results
      NULL
      
      $tags
      $tags$normalized
      [1] TRUE
      
      $tags$norm_method
      [1] "log2"
      
      
      attr(,"class")
      [1] "DAList"

---

    Code
      fit_limma_model(b)
    Output
      $data
                 sample1  sample2  sample3  sample4  sample5  sample6
      protein1  0.000000 3.459432 4.392317 4.954196 5.357552 5.672425
      protein2  1.000000 3.584963 4.459432 5.000000 5.392317 5.700440
      protein3  1.584963 3.700440 4.523562 5.044394 5.426265 5.727920
      protein4  2.000000 3.807355 4.584963 5.087463 5.459432 5.754888
      protein5  2.321928 3.906891 4.643856 5.129283 5.491853 5.781360
      protein6  2.584963 4.000000 4.700440 5.169925 5.523562 5.807355
      protein7  2.807355 4.087463 4.754888 5.209453 5.554589 5.832890
      protein8  3.000000 4.169925 4.807355 5.247928 5.584963 5.857981
      protein9  3.169925 4.247928 4.857981 5.285402 5.614710 5.882643
      protein10 3.321928 4.321928 4.906891 5.321928 5.643856 5.906891
      
      $annotation
                uniprot_id
      protein1    protein1
      protein2    protein2
      protein3    protein3
      protein4    protein4
      protein5    protein5
      protein6    protein6
      protein7    protein7
      protein8    protein8
      protein9    protein9
      protein10  protein10
      
      $metadata
              sample_ID group sex treatment
      sample1   sample1     A   M   control
      sample2   sample2     B   F   control
      sample3   sample3     A   M   control
      sample4   sample4     B   F treatment
      sample5   sample5     A   M treatment
      sample6   sample6     B   F treatment
      
      $design
      $design$design_formula
      [1] "~group"
      
      $design$design_matrix
              Intercept B
      sample1         1 0
      sample2         1 1
      sample3         1 0
      sample4         1 1
      sample5         1 0
      sample6         1 1
      attr(,"assign")
      [1] 0 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      
      $eBayes_fit
      An object of class "MArrayLM"
      $coefficients
                Intercept         B
      protein1   3.249956 1.4453946
      protein2   3.617250 1.1445511
      protein3   3.844930 0.9793217
      protein4   4.014798 0.8684370
      protein5   4.152546 0.7866320
      protein6   4.269655 0.7227719
      protein7   4.372277 0.6709916
      protein8   4.464106 0.6278387
      protein9   4.547539 0.5911190
      protein10  4.624225 0.5593573
      
      $rank
      [1] 2
      
      $assign
      [1] 0 1
      
      $qr
      $qr
               Intercept          B
      sample1 -2.4494897 -1.2247449
      sample2  0.4082483 -1.2247449
      sample3  0.4082483 -0.2898979
      sample4  0.4082483  0.5265986
      sample5  0.4082483 -0.2898979
      sample6  0.4082483  0.5265986
      attr(,"assign")
      [1] 0 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      $qraux
      [1] 1.408248 1.526599
      
      $pivot
      [1] 1 2
      
      $tol
      [1] 1e-07
      
      $rank
      [1] 2
      
      
      $df.residual
       [1] 4 4 4 4 4 4 4 4 4 4
      
      $sigma
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       2.171309  1.805053  1.596609  1.451764  1.341391  1.252654  1.178762  1.115685 
       protein9 protein10 
       1.060837  1.012458 
      
      $cov.coefficients
                 Intercept          B
      Intercept  0.3333333 -0.3333333
      B         -0.3333333  0.6666667
      
      $stdev.unscaled
                Intercept         B
      protein1  0.5773503 0.8164966
      protein2  0.5773503 0.8164966
      protein3  0.5773503 0.8164966
      protein4  0.5773503 0.8164966
      protein5  0.5773503 0.8164966
      protein6  0.5773503 0.8164966
      protein7  0.5773503 0.8164966
      protein8  0.5773503 0.8164966
      protein9  0.5773503 0.8164966
      protein10 0.5773503 0.8164966
      
      $pivot
      [1] 1 2
      
      $Amean
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       3.972654  4.189525  4.334591  4.449017  4.545862  4.631041  4.707773  4.778025 
       protein9 protein10 
       4.843098  4.903904 
      
      $method
      [1] "ls"
      
      $design
              Intercept B
      sample1         1 0
      sample2         1 1
      sample3         1 0
      sample4         1 1
      sample5         1 0
      sample6         1 1
      attr(,"assign")
      [1] 0 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      $df.prior
       [1] Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf
      
      $s2.prior
      [1] 2.347193
      
      $var.prior
      [1] 6.816651582 0.004260407
      
      $proportion
      [1] 0.01
      
      $s2.post
       [1] 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193
       [9] 2.347193 2.347193
      
      $t
                Intercept         B
      protein1   3.674208 1.1554672
      protein2   4.089448 0.9149690
      protein3   4.346849 0.7828825
      protein4   4.538892 0.6942398
      protein5   4.694622 0.6288438
      protein6   4.827018 0.5777932
      protein7   4.943037 0.5363994
      protein8   5.046853 0.5019024
      protein9   5.141177 0.4725482
      protein10  5.227874 0.4471575
      
      $df.total
       [1] 40 40 40 40 40 40 40 40 40 40
      
      $p.value
                   Intercept         B
      protein1  6.993473e-04 0.2547547
      protein2  2.029093e-04 0.3656909
      protein3  9.224506e-05 0.4383068
      protein4  5.078967e-05 0.4915445
      protein5  3.116181e-05 0.5330262
      protein6  2.051361e-05 0.5666423
      protein7  1.419448e-05 0.5946534
      protein8  1.019662e-05 0.6184878
      protein9  7.542600e-06 0.6391038
      protein10 5.713087e-06 0.6571734
      
      $lods
                Intercept         B
      protein1  0.3072394 -4.594066
      protein2  1.8439827 -4.595647
      protein3  2.8791225 -4.596359
      protein4  3.6925675 -4.596775
      protein5  4.3780142 -4.597049
      protein6  4.9789443 -4.597245
      protein7  5.5192770 -4.597391
      protein8  6.0136575 -4.597505
      protein9  6.4717459 -4.597596
      protein10 6.9002739 -4.597670
      
      $F
       [1] 20.83885 22.85234 24.32068 25.53981 26.60994 27.57822 28.47104 29.30487
       [9] 30.09077 30.83660
      
      $F.p.value
       [1] 8.908471e-10 1.189474e-10 2.739460e-11 8.094687e-12 2.776202e-12
       [6] 1.054222e-12 4.317003e-13 1.875231e-13 8.545621e-14 4.053532e-14
      
      
      $results
      NULL
      
      $tags
      $tags$normalized
      [1] TRUE
      
      $tags$norm_method
      [1] "log2"
      
      
      attr(,"class")
      [1] "DAList"

---

    Code
      fit_limma_model(c)
    Output
      $data
                 sample1  sample2  sample3  sample4  sample5  sample6
      protein1  0.000000 3.459432 4.392317 4.954196 5.357552 5.672425
      protein2  1.000000 3.584963 4.459432 5.000000 5.392317 5.700440
      protein3  1.584963 3.700440 4.523562 5.044394 5.426265 5.727920
      protein4  2.000000 3.807355 4.584963 5.087463 5.459432 5.754888
      protein5  2.321928 3.906891 4.643856 5.129283 5.491853 5.781360
      protein6  2.584963 4.000000 4.700440 5.169925 5.523562 5.807355
      protein7  2.807355 4.087463 4.754888 5.209453 5.554589 5.832890
      protein8  3.000000 4.169925 4.807355 5.247928 5.584963 5.857981
      protein9  3.169925 4.247928 4.857981 5.285402 5.614710 5.882643
      protein10 3.321928 4.321928 4.906891 5.321928 5.643856 5.906891
      
      $annotation
                uniprot_id
      protein1    protein1
      protein2    protein2
      protein3    protein3
      protein4    protein4
      protein5    protein5
      protein6    protein6
      protein7    protein7
      protein8    protein8
      protein9    protein9
      protein10  protein10
      
      $metadata
              sample_ID group sex treatment
      sample1   sample1     A   M   control
      sample2   sample2     B   F   control
      sample3   sample3     A   M   control
      sample4   sample4     B   F treatment
      sample5   sample5     A   M treatment
      sample6   sample6     B   F treatment
      
      $design
      $design$design_formula
      [1] "~0 + treatment"
      
      $design$design_matrix
              control treatment
      sample1       1         0
      sample2       1         0
      sample3       1         0
      sample4       0         1
      sample5       0         1
      sample6       0         1
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$treatment
      [1] "contr.treatment"
      
      
      $design$contrast_matrix
                 Contrasts
      Levels      Treatment_vs_Control
        control                     -1
        treatment                    1
      
      $design$contrast_vector
      [1] "Treatment_vs_Control= treatment - control"
      
      
      $eBayes_fit
      An object of class "MArrayLM"
      $coefficients
                 Contrasts
                  Treatment_vs_Control
        protein1              2.710808
        protein2              2.349454
        protein3              2.129872
        protein4              1.969822
        protein5              1.843274
        protein6              1.738480
        protein7              1.649076
        protein8              1.571197
        protein9              1.502307
        protein10             1.440643
      
      $rank
      [1] 2
      
      $assign
      [1] 1 1
      
      $qr
      $qr
                 control  treatment
      sample1 -1.7320508  0.0000000
      sample2  0.5773503 -1.7320508
      sample3  0.5773503  0.0000000
      sample4  0.0000000  0.5773503
      sample5  0.0000000  0.5773503
      sample6  0.0000000  0.5773503
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$treatment
      [1] "contr.treatment"
      
      
      $qraux
      [1] 1.57735 1.00000
      
      $pivot
      [1] 1 2
      
      $tol
      [1] 1e-07
      
      $rank
      [1] 2
      
      
      $df.residual
       [1] 4 4 4 4 4 4 4 4 4 4
      
      $sigma
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
      1.6560014 1.2959508 1.0989453 0.9671413 0.8702035 0.7947801 0.7338403 0.6832441 
       protein9 protein10 
      0.6403590 0.6034134 
      
      $cov.coefficients
                            Contrasts
      Contrasts              Treatment_vs_Control
        Treatment_vs_Control            0.6666667
      
      $stdev.unscaled
                 Contrasts
                  Treatment_vs_Control
        protein1             0.8164966
        protein2             0.8164966
        protein3             0.8164966
        protein4             0.8164966
        protein5             0.8164966
        protein6             0.8164966
        protein7             0.8164966
        protein8             0.8164966
        protein9             0.8164966
        protein10            0.8164966
      
      $pivot
      [1] 1 2
      
      $Amean
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       3.972654  4.189525  4.334591  4.449017  4.545862  4.631041  4.707773  4.778025 
       protein9 protein10 
       4.843098  4.903904 
      
      $method
      [1] "ls"
      
      $design
              control treatment
      sample1       1         0
      sample2       1         0
      sample3       1         0
      sample4       0         1
      sample5       0         1
      sample6       0         1
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$treatment
      [1] "contr.treatment"
      
      
      $contrasts
                 Contrasts
      Levels      Treatment_vs_Control
        control                     -1
        treatment                    1
      
      $df.prior
       [1] 20.76311 39.60907      Inf      Inf      Inf      Inf      Inf      Inf
       [9]      Inf      Inf
      
      $s2.prior
      [1] 0.9922625
      
      $var.prior
      [1] 8.071417
      
      $proportion
      [1] 0.01
      
      $s2.post
       [1] 1.2749537 1.0552977 0.9922625 0.9922625 0.9922625 0.9922625 0.9922625
       [8] 0.9922625 0.9922625 0.9922625
      
      $t
                 Contrasts
                  Treatment_vs_Control
        protein1              2.940338
        protein2              2.801078
        protein3              2.618700
        protein4              2.421917
        protein5              2.266325
        protein6              2.137480
        protein7              2.027556
        protein8              1.931804
        protein9              1.847103
        protein10             1.771286
      
      $df.total
       [1] 24.76311 40.00000 40.00000 40.00000 40.00000 40.00000 40.00000 40.00000
       [9] 40.00000 40.00000
      
      $p.value
                 Contrasts
                  Treatment_vs_Control
        protein1           0.007007173
        protein2           0.007809268
        protein3           0.012406936
        protein4           0.020064186
        protein5           0.028912479
        protein6           0.038723564
        protein7           0.049304122
        protein8           0.060489317
        protein9           0.072138609
        protein10          0.084132465
      
      $lods
                 Contrasts
                  Treatment_vs_Control
        protein1             -2.362805
        protein2             -2.514481
        protein3             -2.714500
        protein4             -3.172617
        protein5             -3.509516
        protein6             -3.771576
        protein7             -3.983028
        protein8             -4.158125
        protein9             -4.305954
        protein10            -4.432656
      
      $F
       [1] 8.645586 7.846036 6.857591 5.865681 5.136228 4.568820 4.110984 3.731865
       [9] 3.411789 3.137453
      
      $F.p.value
       [1] 0.007007173 0.007564511 0.008826548 0.015438887 0.023431505 0.032558993
       [7] 0.042605557 0.053383745 0.064732233 0.076513205
      
      
      $results
      NULL
      
      $tags
      $tags$normalized
      [1] TRUE
      
      $tags$norm_method
      [1] "log2"
      
      
      attr(,"class")
      [1] "DAList"

---

    Code
      fit_limma_model(d)
    Message
      Estimated intra-block correlation = 0.549
    Output
      $data
                 sample1  sample2  sample3  sample4  sample5  sample6
      protein1  0.000000 3.459432 4.392317 4.954196 5.357552 5.672425
      protein2  1.000000 3.584963 4.459432 5.000000 5.392317 5.700440
      protein3  1.584963 3.700440 4.523562 5.044394 5.426265 5.727920
      protein4  2.000000 3.807355 4.584963 5.087463 5.459432 5.754888
      protein5  2.321928 3.906891 4.643856 5.129283 5.491853 5.781360
      protein6  2.584963 4.000000 4.700440 5.169925 5.523562 5.807355
      protein7  2.807355 4.087463 4.754888 5.209453 5.554589 5.832890
      protein8  3.000000 4.169925 4.807355 5.247928 5.584963 5.857981
      protein9  3.169925 4.247928 4.857981 5.285402 5.614710 5.882643
      protein10 3.321928 4.321928 4.906891 5.321928 5.643856 5.906891
      
      $annotation
                uniprot_id
      protein1    protein1
      protein2    protein2
      protein3    protein3
      protein4    protein4
      protein5    protein5
      protein6    protein6
      protein7    protein7
      protein8    protein8
      protein9    protein9
      protein10  protein10
      
      $metadata
              sample_ID group sex treatment
      sample1   sample1     A   M   control
      sample2   sample2     B   F   control
      sample3   sample3     A   M   control
      sample4   sample4     B   F treatment
      sample5   sample5     A   M treatment
      sample6   sample6     B   F treatment
      
      $design
      $design$design_formula
      [1] "~0 + sex + (1 | treatment)"
      
      $design$design_matrix
              F M
      sample1 0 1
      sample2 1 0
      sample3 0 1
      sample4 1 0
      sample5 0 1
      sample6 1 0
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$sex
      [1] "contr.treatment"
      
      
      $design$random_factor
      [1] "treatment"
      
      
      $eBayes_fit
      An object of class "MArrayLM"
      $coefficients
                       F        M
      protein1  4.375960 3.569347
      protein2  4.479819 3.899232
      protein3  4.565841 4.103340
      protein4  4.642462 4.255571
      protein5  4.712630 4.379093
      protein6  4.777845 4.484237
      protein7  4.839024 4.576522
      protein8  4.896798 4.659252
      protein9  4.951628 4.734568
      protein10 5.003872 4.803936
      
      $rank
      [1] 2
      
      $assign
      NULL
      
      $qr
      $qr
                       F           M
      sample1 -1.9361078  1.19726894
      sample2  0.6177629  1.52153232
      sample3 -0.2340241 -0.51247009
      sample4  0.5165002 -0.06696502
      sample5 -0.3389079 -0.74214627
      sample6  0.4265803  0.24248249
      
      $qraux
      [1] 1.000000 1.351157
      
      $pivot
      [1] 1 2
      
      $tol
      [1] 1e-07
      
      $rank
      [1] 2
      
      
      $df.residual
       [1] 4 4 4 4 4 4 4 4 4 4
      
      $sigma
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       2.625995  2.109991  1.823286  1.628639  1.483495  1.369098  1.275549  1.196999 
       protein9 protein10 
       1.129716  1.071179 
      
      $cov.coefficients
                F         M
      F 0.4319541 0.2671159
      M 0.2671159 0.4319541
      
      $stdev.unscaled
                        F         M
      protein1  0.6572322 0.6572322
      protein2  0.6572322 0.6572322
      protein3  0.6572322 0.6572322
      protein4  0.6572322 0.6572322
      protein5  0.6572322 0.6572322
      protein6  0.6572322 0.6572322
      protein7  0.6572322 0.6572322
      protein8  0.6572322 0.6572322
      protein9  0.6572322 0.6572322
      protein10 0.6572322 0.6572322
      
      $pivot
      [1] 1 2
      
      $ndups
      [1] 1
      
      $spacing
      [1] 1
      
      $block
      [1] "control"   "control"   "control"   "treatment" "treatment" "treatment"
      
      $correlation
      [1] 0.5486051
      
      $Amean
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       3.972654  4.189525  4.334591  4.449017  4.545862  4.631041  4.707773  4.778025 
       protein9 protein10 
       4.843098  4.903904 
      
      $method
      [1] "ls"
      
      $design
              F M
      sample1 0 1
      sample2 1 0
      sample3 0 1
      sample4 1 0
      sample5 0 1
      sample6 1 0
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$sex
      [1] "contr.treatment"
      
      
      $df.prior
       [1] 38.58884      Inf      Inf      Inf      Inf      Inf      Inf      Inf
       [9]      Inf      Inf
      
      $s2.prior
      [1] 2.880884
      
      $var.prior
      [1] 5.553852 5.553852
      
      $proportion
      [1] 0.01
      
      $s2.post
       [1] 3.257975 2.880884 2.880884 2.880884 2.880884 2.880884 2.880884 2.880884
       [9] 2.880884 2.880884
      
      $t
                       F        M
      protein1  3.688763 3.008819
      protein2  4.015862 3.495404
      protein3  4.092975 3.678374
      protein4  4.161661 3.814839
      protein5  4.224562 3.925568
      protein6  4.283022 4.019822
      protein7  4.337866 4.102550
      protein8  4.389656 4.176712
      protein9  4.438808 4.244228
      protein10 4.485641 4.306411
      
      $df.total
       [1] 40 40 40 40 40 40 40 40 40 40
      
      $p.value
                           F            M
      protein1  6.702158e-04 0.0045222260
      protein2  2.534988e-04 0.0011730460
      protein3  2.007491e-04 0.0006908872
      protein4  1.628886e-04 0.0004623830
      protein5  1.343898e-04 0.0003325038
      protein6  1.123071e-04 0.0002504882
      protein7  9.484042e-05 0.0001949980
      protein8  8.080341e-05 0.0001555746
      protein9  6.937542e-05 0.0001265250
      protein10 5.996918e-05 0.0001045035
      
      $lods
                         F          M
      protein1  -0.4042978 -2.0591941
      protein2   1.5721473 -0.2414474
      protein3   1.8622358  0.3674842
      protein4   2.1252654  0.8418699
      protein5   2.3699841  1.2394903
      protein6   2.6007173  1.5869121
      protein7   2.8200573  1.8986387
      protein8   3.0297491  2.1834902
      protein9   3.2310595  2.4472479
      protein10  3.4249574  2.6939172
      
      $F
       [1]  7.232234  8.892773  9.441926  9.907174 10.319479 10.694430 11.041138
       [8] 11.365461 11.671429 11.961954
      
      $F.p.value
       [1] 1.976640e-03 1.373782e-04 7.932749e-05 4.981602e-05 3.298429e-05
       [6] 2.267086e-05 1.602857e-05 1.158892e-05 8.534199e-06 6.382476e-06
      
      
      $results
      NULL
      
      $tags
      $tags$normalized
      [1] TRUE
      
      $tags$norm_method
      [1] "log2"
      
      
      attr(,"class")
      [1] "DAList"

# extract_DA_results gives consistent results

    Code
      extract_DA_results(a)
    Message
      ! Warning: more than 20% of the data is differentially adundant in some terms
      ! Criteria for DA: |logFC| > 1, p-value < 0.05, p.value adjustment = "none"
      !  Problematic terms: "A" and "B"
      ! Assumption that most proteins are not DA may be violated
      ! Warning: more than 20% of the data is differentially adundant in some terms
      ! Criteria for DA: |logFC| > 1, p-value < 0.05, p.value adjustment = "BH"
      !  Problematic terms: "A" and "B"
      ! Assumption that most proteins are not DA may be violated
    Output
      $data
                 sample1  sample2  sample3  sample4  sample5  sample6
      protein1  0.000000 3.459432 4.392317 4.954196 5.357552 5.672425
      protein2  1.000000 3.584963 4.459432 5.000000 5.392317 5.700440
      protein3  1.584963 3.700440 4.523562 5.044394 5.426265 5.727920
      protein4  2.000000 3.807355 4.584963 5.087463 5.459432 5.754888
      protein5  2.321928 3.906891 4.643856 5.129283 5.491853 5.781360
      protein6  2.584963 4.000000 4.700440 5.169925 5.523562 5.807355
      protein7  2.807355 4.087463 4.754888 5.209453 5.554589 5.832890
      protein8  3.000000 4.169925 4.807355 5.247928 5.584963 5.857981
      protein9  3.169925 4.247928 4.857981 5.285402 5.614710 5.882643
      protein10 3.321928 4.321928 4.906891 5.321928 5.643856 5.906891
      
      $annotation
                uniprot_id
      protein1    protein1
      protein2    protein2
      protein3    protein3
      protein4    protein4
      protein5    protein5
      protein6    protein6
      protein7    protein7
      protein8    protein8
      protein9    protein9
      protein10  protein10
      
      $metadata
              sample_ID group sex treatment
      sample1   sample1     A   M   control
      sample2   sample2     B   F   control
      sample3   sample3     A   M   control
      sample4   sample4     B   F treatment
      sample5   sample5     A   M treatment
      sample6   sample6     B   F treatment
      
      $design
      $design$design_formula
      [1] "~0 + group"
      
      $design$design_matrix
              A B
      sample1 1 0
      sample2 0 1
      sample3 1 0
      sample4 0 1
      sample5 1 0
      sample6 0 1
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      
      $eBayes_fit
      An object of class "MArrayLM"
      $coefficients
                       A        B
      protein1  3.249956 4.695351
      protein2  3.617250 4.761801
      protein3  3.844930 4.824251
      protein4  4.014798 4.883235
      protein5  4.152546 4.939178
      protein6  4.269655 4.992427
      protein7  4.372277 5.043269
      protein8  4.464106 5.091945
      protein9  4.547539 5.138658
      protein10 4.624225 5.183582
      
      $rank
      [1] 2
      
      $assign
      [1] 1 1
      
      $qr
      $qr
                       A          B
      sample1 -1.7320508  0.0000000
      sample2  0.0000000 -1.7320508
      sample3  0.5773503  0.0000000
      sample4  0.0000000  0.5773503
      sample5  0.5773503  0.0000000
      sample6  0.0000000  0.5773503
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      $qraux
      [1] 1.57735 1.57735
      
      $pivot
      [1] 1 2
      
      $tol
      [1] 1e-07
      
      $rank
      [1] 2
      
      
      $df.residual
       [1] 4 4 4 4 4 4 4 4 4 4
      
      $sigma
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       2.171309  1.805053  1.596609  1.451764  1.341391  1.252654  1.178762  1.115685 
       protein9 protein10 
       1.060837  1.012458 
      
      $cov.coefficients
                A         B
      A 0.3333333 0.0000000
      B 0.0000000 0.3333333
      
      $stdev.unscaled
                        A         B
      protein1  0.5773503 0.5773503
      protein2  0.5773503 0.5773503
      protein3  0.5773503 0.5773503
      protein4  0.5773503 0.5773503
      protein5  0.5773503 0.5773503
      protein6  0.5773503 0.5773503
      protein7  0.5773503 0.5773503
      protein8  0.5773503 0.5773503
      protein9  0.5773503 0.5773503
      protein10 0.5773503 0.5773503
      
      $pivot
      [1] 1 2
      
      $Amean
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       3.972654  4.189525  4.334591  4.449017  4.545862  4.631041  4.707773  4.778025 
       protein9 protein10 
       4.843098  4.903904 
      
      $method
      [1] "ls"
      
      $design
              A B
      sample1 1 0
      sample2 0 1
      sample3 1 0
      sample4 0 1
      sample5 1 0
      sample6 0 1
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      $df.prior
       [1] Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf
      
      $s2.prior
      [1] 2.347193
      
      $var.prior
      [1] 6.816652 6.816652
      
      $proportion
      [1] 0.01
      
      $s2.post
       [1] 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193
       [9] 2.347193 2.347193
      
      $t
                       A        B
      protein1  3.674208 5.308285
      protein2  4.089448 5.383409
      protein3  4.346849 5.454012
      protein4  4.538892 5.520696
      protein5  4.694622 5.583941
      protein6  4.827018 5.644141
      protein7  4.943037 5.701620
      protein8  5.046853 5.756650
      protein9  5.141177 5.809461
      protein10 5.227874 5.860250
      
      $df.total
       [1] 40 40 40 40 40 40 40 40 40 40
      
      $p.value
                           A            B
      protein1  6.993473e-04 4.413049e-06
      protein2  2.029093e-04 3.465792e-06
      protein3  9.224506e-05 2.760859e-06
      protein4  5.078967e-05 2.226698e-06
      protein5  3.116181e-05 1.815553e-06
      protein6  2.051361e-05 1.494720e-06
      protein7  1.419448e-05 1.241299e-06
      protein8  1.019662e-05 1.038942e-06
      protein9  7.542600e-06 8.757757e-07
      protein10 5.713087e-06 7.430412e-07
      
      $lods
                        A         B
      protein1  0.3072394  7.304137
      protein2  1.8439827  7.687015
      protein3  2.8791225  8.051757
      protein4  3.6925675  8.400614
      protein5  4.3780142  8.735402
      protein6  4.9789443  9.057611
      protein7  5.5192770  9.368481
      protein8  6.0136575  9.669057
      protein9  6.4717459  9.960228
      protein10 6.9002739 10.242760
      
      $F
       [1] 20.83885 22.85234 24.32068 25.53981 26.60994 27.57822 28.47104 29.30487
       [9] 30.09077 30.83660
      
      $F.p.value
       [1] 8.908471e-10 1.189474e-10 2.739460e-11 8.094687e-12 2.776202e-12
       [6] 1.054222e-12 4.317003e-13 1.875231e-13 8.545621e-14 4.053532e-14
      
      
      $results
      $results$A
                   logFC     CI.L     CI.R average_intensity        t         B
      protein1  3.249956 1.462250 5.037663          3.972654 3.674208 0.3072394
      protein2  3.617250 1.829543 5.404957          4.189525 4.089448 1.8439827
      protein3  3.844930 2.057223 5.632637          4.334591 4.346849 2.8791225
      protein4  4.014798 2.227091 5.802505          4.449017 4.538892 3.6925675
      protein5  4.152546 2.364839 5.940253          4.545862 4.694622 4.3780142
      protein6  4.269655 2.481948 6.057362          4.631041 4.827018 4.9789443
      protein7  4.372277 2.584570 6.159984          4.707773 4.943037 5.5192770
      protein8  4.464106 2.676399 6.251813          4.778025 5.046853 6.0136575
      protein9  4.547539 2.759832 6.335246          4.843098 5.141177 6.4717459
      protein10 4.624225 2.836518 6.411932          4.903904 5.227874 6.9002739
                     P.Value    adj.P.Val sig.PVal sig.FDR
      protein1  6.993473e-04 6.993473e-04        1       1
      protein2  2.029093e-04 2.254548e-04        1       1
      protein3  9.224506e-05 1.153063e-04        1       1
      protein4  5.078967e-05 7.255667e-05        1       1
      protein5  3.116181e-05 5.193635e-05        1       1
      protein6  2.051361e-05 4.102723e-05        1       1
      protein7  1.419448e-05 3.548621e-05        1       1
      protein8  1.019662e-05 3.398872e-05        1       1
      protein9  7.542600e-06 3.398872e-05        1       1
      protein10 5.713087e-06 3.398872e-05        1       1
      
      $results$B
                   logFC     CI.L     CI.R average_intensity        t         B
      protein1  4.695351 2.907644 6.483058          3.972654 5.308285  7.304137
      protein2  4.761801 2.974094 6.549508          4.189525 5.383409  7.687015
      protein3  4.824251 3.036544 6.611958          4.334591 5.454012  8.051757
      protein4  4.883235 3.095528 6.670942          4.449017 5.520696  8.400614
      protein5  4.939178 3.151471 6.726885          4.545862 5.583941  8.735402
      protein6  4.992427 3.204720 6.780134          4.631041 5.644141  9.057611
      protein7  5.043269 3.255562 6.830976          4.707773 5.701620  9.368481
      protein8  5.091945 3.304238 6.879651          4.778025 5.756650  9.669057
      protein9  5.138658 3.350951 6.926365          4.843098 5.809461  9.960228
      protein10 5.183582 3.395875 6.971289          4.903904 5.860250 10.242760
                     P.Value    adj.P.Val sig.PVal sig.FDR
      protein1  4.413049e-06 4.413049e-06        1       1
      protein2  3.465792e-06 3.850880e-06        1       1
      protein3  2.760859e-06 3.451073e-06        1       1
      protein4  2.226698e-06 3.180997e-06        1       1
      protein5  1.815553e-06 3.025921e-06        1       1
      protein6  1.494720e-06 2.989439e-06        1       1
      protein7  1.241299e-06 2.989439e-06        1       1
      protein8  1.038942e-06 2.989439e-06        1       1
      protein9  8.757757e-07 2.989439e-06        1       1
      protein10 7.430412e-07 2.989439e-06        1       1
      
      
      $tags
      $tags$normalized
      [1] TRUE
      
      $tags$norm_method
      [1] "log2"
      
      $tags$DA_criteria
      $tags$DA_criteria$pval_thresh
      [1] 0.05
      
      $tags$DA_criteria$lfc_thresh
      [1] 1
      
      $tags$DA_criteria$adj_method
      [1] "BH"
      
      
      $tags$extract_intercept
      [1] FALSE
      
      
      attr(,"class")
      [1] "DAList"

---

    Code
      extract_DA_results(b)
    Output
      $data
                 sample1  sample2  sample3  sample4  sample5  sample6
      protein1  0.000000 3.459432 4.392317 4.954196 5.357552 5.672425
      protein2  1.000000 3.584963 4.459432 5.000000 5.392317 5.700440
      protein3  1.584963 3.700440 4.523562 5.044394 5.426265 5.727920
      protein4  2.000000 3.807355 4.584963 5.087463 5.459432 5.754888
      protein5  2.321928 3.906891 4.643856 5.129283 5.491853 5.781360
      protein6  2.584963 4.000000 4.700440 5.169925 5.523562 5.807355
      protein7  2.807355 4.087463 4.754888 5.209453 5.554589 5.832890
      protein8  3.000000 4.169925 4.807355 5.247928 5.584963 5.857981
      protein9  3.169925 4.247928 4.857981 5.285402 5.614710 5.882643
      protein10 3.321928 4.321928 4.906891 5.321928 5.643856 5.906891
      
      $annotation
                uniprot_id
      protein1    protein1
      protein2    protein2
      protein3    protein3
      protein4    protein4
      protein5    protein5
      protein6    protein6
      protein7    protein7
      protein8    protein8
      protein9    protein9
      protein10  protein10
      
      $metadata
              sample_ID group sex treatment
      sample1   sample1     A   M   control
      sample2   sample2     B   F   control
      sample3   sample3     A   M   control
      sample4   sample4     B   F treatment
      sample5   sample5     A   M treatment
      sample6   sample6     B   F treatment
      
      $design
      $design$design_formula
      [1] "~group"
      
      $design$design_matrix
              Intercept B
      sample1         1 0
      sample2         1 1
      sample3         1 0
      sample4         1 1
      sample5         1 0
      sample6         1 1
      attr(,"assign")
      [1] 0 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      
      $eBayes_fit
      An object of class "MArrayLM"
      $coefficients
                Intercept         B
      protein1   3.249956 1.4453946
      protein2   3.617250 1.1445511
      protein3   3.844930 0.9793217
      protein4   4.014798 0.8684370
      protein5   4.152546 0.7866320
      protein6   4.269655 0.7227719
      protein7   4.372277 0.6709916
      protein8   4.464106 0.6278387
      protein9   4.547539 0.5911190
      protein10  4.624225 0.5593573
      
      $rank
      [1] 2
      
      $assign
      [1] 0 1
      
      $qr
      $qr
               Intercept          B
      sample1 -2.4494897 -1.2247449
      sample2  0.4082483 -1.2247449
      sample3  0.4082483 -0.2898979
      sample4  0.4082483  0.5265986
      sample5  0.4082483 -0.2898979
      sample6  0.4082483  0.5265986
      attr(,"assign")
      [1] 0 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      $qraux
      [1] 1.408248 1.526599
      
      $pivot
      [1] 1 2
      
      $tol
      [1] 1e-07
      
      $rank
      [1] 2
      
      
      $df.residual
       [1] 4 4 4 4 4 4 4 4 4 4
      
      $sigma
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       2.171309  1.805053  1.596609  1.451764  1.341391  1.252654  1.178762  1.115685 
       protein9 protein10 
       1.060837  1.012458 
      
      $cov.coefficients
                 Intercept          B
      Intercept  0.3333333 -0.3333333
      B         -0.3333333  0.6666667
      
      $stdev.unscaled
                Intercept         B
      protein1  0.5773503 0.8164966
      protein2  0.5773503 0.8164966
      protein3  0.5773503 0.8164966
      protein4  0.5773503 0.8164966
      protein5  0.5773503 0.8164966
      protein6  0.5773503 0.8164966
      protein7  0.5773503 0.8164966
      protein8  0.5773503 0.8164966
      protein9  0.5773503 0.8164966
      protein10 0.5773503 0.8164966
      
      $pivot
      [1] 1 2
      
      $Amean
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       3.972654  4.189525  4.334591  4.449017  4.545862  4.631041  4.707773  4.778025 
       protein9 protein10 
       4.843098  4.903904 
      
      $method
      [1] "ls"
      
      $design
              Intercept B
      sample1         1 0
      sample2         1 1
      sample3         1 0
      sample4         1 1
      sample5         1 0
      sample6         1 1
      attr(,"assign")
      [1] 0 1
      attr(,"contrasts")
      attr(,"contrasts")$group
      [1] "contr.treatment"
      
      
      $df.prior
       [1] Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf
      
      $s2.prior
      [1] 2.347193
      
      $var.prior
      [1] 6.816651582 0.004260407
      
      $proportion
      [1] 0.01
      
      $s2.post
       [1] 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193 2.347193
       [9] 2.347193 2.347193
      
      $t
                Intercept         B
      protein1   3.674208 1.1554672
      protein2   4.089448 0.9149690
      protein3   4.346849 0.7828825
      protein4   4.538892 0.6942398
      protein5   4.694622 0.6288438
      protein6   4.827018 0.5777932
      protein7   4.943037 0.5363994
      protein8   5.046853 0.5019024
      protein9   5.141177 0.4725482
      protein10  5.227874 0.4471575
      
      $df.total
       [1] 40 40 40 40 40 40 40 40 40 40
      
      $p.value
                   Intercept         B
      protein1  6.993473e-04 0.2547547
      protein2  2.029093e-04 0.3656909
      protein3  9.224506e-05 0.4383068
      protein4  5.078967e-05 0.4915445
      protein5  3.116181e-05 0.5330262
      protein6  2.051361e-05 0.5666423
      protein7  1.419448e-05 0.5946534
      protein8  1.019662e-05 0.6184878
      protein9  7.542600e-06 0.6391038
      protein10 5.713087e-06 0.6571734
      
      $lods
                Intercept         B
      protein1  0.3072394 -4.594066
      protein2  1.8439827 -4.595647
      protein3  2.8791225 -4.596359
      protein4  3.6925675 -4.596775
      protein5  4.3780142 -4.597049
      protein6  4.9789443 -4.597245
      protein7  5.5192770 -4.597391
      protein8  6.0136575 -4.597505
      protein9  6.4717459 -4.597596
      protein10 6.9002739 -4.597670
      
      $F
       [1] 20.83885 22.85234 24.32068 25.53981 26.60994 27.57822 28.47104 29.30487
       [9] 30.09077 30.83660
      
      $F.p.value
       [1] 8.908471e-10 1.189474e-10 2.739460e-11 8.094687e-12 2.776202e-12
       [6] 1.054222e-12 4.317003e-13 1.875231e-13 8.545621e-14 4.053532e-14
      
      
      $results
      $results$B
                    logFC      CI.L     CI.R average_intensity         t         B
      protein1  1.4453946 -1.082805 3.973594          3.972654 1.1554672 -4.594066
      protein2  1.1445511 -1.383648 3.672750          4.189525 0.9149690 -4.595647
      protein3  0.9793217 -1.548878 3.507521          4.334591 0.7828825 -4.596359
      protein4  0.8684370 -1.659762 3.396636          4.449017 0.6942398 -4.596775
      protein5  0.7866320 -1.741567 3.314831          4.545862 0.6288438 -4.597049
      protein6  0.7227719 -1.805428 3.250971          4.631041 0.5777932 -4.597245
      protein7  0.6709916 -1.857208 3.199191          4.707773 0.5363994 -4.597391
      protein8  0.6278387 -1.900361 3.156038          4.778025 0.5019024 -4.597505
      protein9  0.5911190 -1.937080 3.119318          4.843098 0.4725482 -4.597596
      protein10 0.5593573 -1.968842 3.087557          4.903904 0.4471575 -4.597670
                  P.Value adj.P.Val sig.PVal sig.FDR
      protein1  0.2547547 0.6571734        0       0
      protein2  0.3656909 0.6571734        0       0
      protein3  0.4383068 0.6571734        0       0
      protein4  0.4915445 0.6571734        0       0
      protein5  0.5330262 0.6571734        0       0
      protein6  0.5666423 0.6571734        0       0
      protein7  0.5946534 0.6571734        0       0
      protein8  0.6184878 0.6571734        0       0
      protein9  0.6391038 0.6571734        0       0
      protein10 0.6571734 0.6571734        0       0
      
      
      $tags
      $tags$normalized
      [1] TRUE
      
      $tags$norm_method
      [1] "log2"
      
      $tags$DA_criteria
      $tags$DA_criteria$pval_thresh
      [1] 0.05
      
      $tags$DA_criteria$lfc_thresh
      [1] 1
      
      $tags$DA_criteria$adj_method
      [1] "BH"
      
      
      $tags$extract_intercept
      [1] FALSE
      
      
      attr(,"class")
      [1] "DAList"

---

    Code
      extract_DA_results(c)
    Message
      ! Warning: more than 20% of the data is differentially adundant in a term
      ! Criteria for DA: |logFC| > 1, p-value < 0.05, p.value adjustment = "none"
      !  Problematic term: "Treatment_vs_Control"
      ! Assumption that most proteins are not DA may be violated
      ! Warning: more than 20% of the data is differentially adundant in a term
      ! Criteria for DA: |logFC| > 1, p-value < 0.05, p.value adjustment = "BH"
      !  Problematic term: "Treatment_vs_Control"
      ! Assumption that most proteins are not DA may be violated
    Output
      $data
                 sample1  sample2  sample3  sample4  sample5  sample6
      protein1  0.000000 3.459432 4.392317 4.954196 5.357552 5.672425
      protein2  1.000000 3.584963 4.459432 5.000000 5.392317 5.700440
      protein3  1.584963 3.700440 4.523562 5.044394 5.426265 5.727920
      protein4  2.000000 3.807355 4.584963 5.087463 5.459432 5.754888
      protein5  2.321928 3.906891 4.643856 5.129283 5.491853 5.781360
      protein6  2.584963 4.000000 4.700440 5.169925 5.523562 5.807355
      protein7  2.807355 4.087463 4.754888 5.209453 5.554589 5.832890
      protein8  3.000000 4.169925 4.807355 5.247928 5.584963 5.857981
      protein9  3.169925 4.247928 4.857981 5.285402 5.614710 5.882643
      protein10 3.321928 4.321928 4.906891 5.321928 5.643856 5.906891
      
      $annotation
                uniprot_id
      protein1    protein1
      protein2    protein2
      protein3    protein3
      protein4    protein4
      protein5    protein5
      protein6    protein6
      protein7    protein7
      protein8    protein8
      protein9    protein9
      protein10  protein10
      
      $metadata
              sample_ID group sex treatment
      sample1   sample1     A   M   control
      sample2   sample2     B   F   control
      sample3   sample3     A   M   control
      sample4   sample4     B   F treatment
      sample5   sample5     A   M treatment
      sample6   sample6     B   F treatment
      
      $design
      $design$design_formula
      [1] "~0 + treatment"
      
      $design$design_matrix
              control treatment
      sample1       1         0
      sample2       1         0
      sample3       1         0
      sample4       0         1
      sample5       0         1
      sample6       0         1
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$treatment
      [1] "contr.treatment"
      
      
      $design$contrast_matrix
                 Contrasts
      Levels      Treatment_vs_Control
        control                     -1
        treatment                    1
      
      $design$contrast_vector
      [1] "Treatment_vs_Control= treatment - control"
      
      
      $eBayes_fit
      An object of class "MArrayLM"
      $coefficients
                 Contrasts
                  Treatment_vs_Control
        protein1              2.710808
        protein2              2.349454
        protein3              2.129872
        protein4              1.969822
        protein5              1.843274
        protein6              1.738480
        protein7              1.649076
        protein8              1.571197
        protein9              1.502307
        protein10             1.440643
      
      $rank
      [1] 2
      
      $assign
      [1] 1 1
      
      $qr
      $qr
                 control  treatment
      sample1 -1.7320508  0.0000000
      sample2  0.5773503 -1.7320508
      sample3  0.5773503  0.0000000
      sample4  0.0000000  0.5773503
      sample5  0.0000000  0.5773503
      sample6  0.0000000  0.5773503
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$treatment
      [1] "contr.treatment"
      
      
      $qraux
      [1] 1.57735 1.00000
      
      $pivot
      [1] 1 2
      
      $tol
      [1] 1e-07
      
      $rank
      [1] 2
      
      
      $df.residual
       [1] 4 4 4 4 4 4 4 4 4 4
      
      $sigma
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
      1.6560014 1.2959508 1.0989453 0.9671413 0.8702035 0.7947801 0.7338403 0.6832441 
       protein9 protein10 
      0.6403590 0.6034134 
      
      $cov.coefficients
                            Contrasts
      Contrasts              Treatment_vs_Control
        Treatment_vs_Control            0.6666667
      
      $stdev.unscaled
                 Contrasts
                  Treatment_vs_Control
        protein1             0.8164966
        protein2             0.8164966
        protein3             0.8164966
        protein4             0.8164966
        protein5             0.8164966
        protein6             0.8164966
        protein7             0.8164966
        protein8             0.8164966
        protein9             0.8164966
        protein10            0.8164966
      
      $pivot
      [1] 1 2
      
      $Amean
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       3.972654  4.189525  4.334591  4.449017  4.545862  4.631041  4.707773  4.778025 
       protein9 protein10 
       4.843098  4.903904 
      
      $method
      [1] "ls"
      
      $design
              control treatment
      sample1       1         0
      sample2       1         0
      sample3       1         0
      sample4       0         1
      sample5       0         1
      sample6       0         1
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$treatment
      [1] "contr.treatment"
      
      
      $contrasts
                 Contrasts
      Levels      Treatment_vs_Control
        control                     -1
        treatment                    1
      
      $df.prior
       [1] 20.76311 39.60907      Inf      Inf      Inf      Inf      Inf      Inf
       [9]      Inf      Inf
      
      $s2.prior
      [1] 0.9922625
      
      $var.prior
      [1] 8.071417
      
      $proportion
      [1] 0.01
      
      $s2.post
       [1] 1.2749537 1.0552977 0.9922625 0.9922625 0.9922625 0.9922625 0.9922625
       [8] 0.9922625 0.9922625 0.9922625
      
      $t
                 Contrasts
                  Treatment_vs_Control
        protein1              2.940338
        protein2              2.801078
        protein3              2.618700
        protein4              2.421917
        protein5              2.266325
        protein6              2.137480
        protein7              2.027556
        protein8              1.931804
        protein9              1.847103
        protein10             1.771286
      
      $df.total
       [1] 24.76311 40.00000 40.00000 40.00000 40.00000 40.00000 40.00000 40.00000
       [9] 40.00000 40.00000
      
      $p.value
                 Contrasts
                  Treatment_vs_Control
        protein1           0.007007173
        protein2           0.007809268
        protein3           0.012406936
        protein4           0.020064186
        protein5           0.028912479
        protein6           0.038723564
        protein7           0.049304122
        protein8           0.060489317
        protein9           0.072138609
        protein10          0.084132465
      
      $lods
                 Contrasts
                  Treatment_vs_Control
        protein1             -2.362805
        protein2             -2.514481
        protein3             -2.714500
        protein4             -3.172617
        protein5             -3.509516
        protein6             -3.771576
        protein7             -3.983028
        protein8             -4.158125
        protein9             -4.305954
        protein10            -4.432656
      
      $F
       [1] 8.645586 7.846036 6.857591 5.865681 5.136228 4.568820 4.110984 3.731865
       [9] 3.411789 3.137453
      
      $F.p.value
       [1] 0.007007173 0.007564511 0.008826548 0.015438887 0.023431505 0.032558993
       [7] 0.042605557 0.053383745 0.064732233 0.076513205
      
      
      $results
      $results$Treatment_vs_Control
                   logFC         CI.L     CI.R average_intensity        t         B
      protein1  2.710808  0.811120551 4.610496          3.972654 2.940338 -2.362805
      protein2  2.349454  0.654240940 4.044668          4.189525 2.801078 -2.514481
      protein3  2.129872  0.486067158 3.773676          4.334591 2.618700 -2.714500
      protein4  1.969822  0.326016953 3.613626          4.449017 2.421917 -3.172617
      protein5  1.843274  0.199469089 3.487078          4.545862 2.266325 -3.509516
      protein6  1.738480  0.094675327 3.382284          4.631041 2.137480 -3.771576
      protein7  1.649076  0.005271096 3.292880          4.707773 2.027556 -3.983028
      protein8  1.571197 -0.072607531 3.215002          4.778025 1.931804 -4.158125
      protein9  1.502307 -0.141497359 3.146112          4.843098 1.847103 -4.305954
      protein10 1.440643 -0.203161861 3.084447          4.903904 1.771286 -4.432656
                    P.Value  adj.P.Val sig.PVal sig.FDR
      protein1  0.007007173 0.03904634        1       1
      protein2  0.007809268 0.03904634        1       1
      protein3  0.012406936 0.04135645        1       1
      protein4  0.020064186 0.05016047        1       0
      protein5  0.028912479 0.05782496        1       0
      protein6  0.038723564 0.06453927        1       0
      protein7  0.049304122 0.07043446        1       0
      protein8  0.060489317 0.07561165        0       0
      protein9  0.072138609 0.08015401        0       0
      protein10 0.084132465 0.08413247        0       0
      
      
      $tags
      $tags$normalized
      [1] TRUE
      
      $tags$norm_method
      [1] "log2"
      
      $tags$DA_criteria
      $tags$DA_criteria$pval_thresh
      [1] 0.05
      
      $tags$DA_criteria$lfc_thresh
      [1] 1
      
      $tags$DA_criteria$adj_method
      [1] "BH"
      
      
      $tags$extract_intercept
      [1] FALSE
      
      
      attr(,"class")
      [1] "DAList"

---

    Code
      extract_DA_results(d)
    Message
      ! Warning: more than 20% of the data is differentially adundant in some terms
      ! Criteria for DA: |logFC| > 1, p-value < 0.05, p.value adjustment = "none"
      !  Problematic terms: "F" and "M"
      ! Assumption that most proteins are not DA may be violated
      ! Warning: more than 20% of the data is differentially adundant in some terms
      ! Criteria for DA: |logFC| > 1, p-value < 0.05, p.value adjustment = "BH"
      !  Problematic terms: "F" and "M"
      ! Assumption that most proteins are not DA may be violated
    Output
      $data
                 sample1  sample2  sample3  sample4  sample5  sample6
      protein1  0.000000 3.459432 4.392317 4.954196 5.357552 5.672425
      protein2  1.000000 3.584963 4.459432 5.000000 5.392317 5.700440
      protein3  1.584963 3.700440 4.523562 5.044394 5.426265 5.727920
      protein4  2.000000 3.807355 4.584963 5.087463 5.459432 5.754888
      protein5  2.321928 3.906891 4.643856 5.129283 5.491853 5.781360
      protein6  2.584963 4.000000 4.700440 5.169925 5.523562 5.807355
      protein7  2.807355 4.087463 4.754888 5.209453 5.554589 5.832890
      protein8  3.000000 4.169925 4.807355 5.247928 5.584963 5.857981
      protein9  3.169925 4.247928 4.857981 5.285402 5.614710 5.882643
      protein10 3.321928 4.321928 4.906891 5.321928 5.643856 5.906891
      
      $annotation
                uniprot_id
      protein1    protein1
      protein2    protein2
      protein3    protein3
      protein4    protein4
      protein5    protein5
      protein6    protein6
      protein7    protein7
      protein8    protein8
      protein9    protein9
      protein10  protein10
      
      $metadata
              sample_ID group sex treatment
      sample1   sample1     A   M   control
      sample2   sample2     B   F   control
      sample3   sample3     A   M   control
      sample4   sample4     B   F treatment
      sample5   sample5     A   M treatment
      sample6   sample6     B   F treatment
      
      $design
      $design$design_formula
      [1] "~0 + sex + (1 | treatment)"
      
      $design$design_matrix
              F M
      sample1 0 1
      sample2 1 0
      sample3 0 1
      sample4 1 0
      sample5 0 1
      sample6 1 0
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$sex
      [1] "contr.treatment"
      
      
      $design$random_factor
      [1] "treatment"
      
      
      $eBayes_fit
      An object of class "MArrayLM"
      $coefficients
                       F        M
      protein1  4.375960 3.569347
      protein2  4.479819 3.899232
      protein3  4.565841 4.103340
      protein4  4.642462 4.255571
      protein5  4.712630 4.379093
      protein6  4.777845 4.484237
      protein7  4.839024 4.576522
      protein8  4.896798 4.659252
      protein9  4.951628 4.734568
      protein10 5.003872 4.803936
      
      $rank
      [1] 2
      
      $assign
      NULL
      
      $qr
      $qr
                       F           M
      sample1 -1.9361078  1.19726894
      sample2  0.6177629  1.52153232
      sample3 -0.2340241 -0.51247009
      sample4  0.5165002 -0.06696502
      sample5 -0.3389079 -0.74214627
      sample6  0.4265803  0.24248249
      
      $qraux
      [1] 1.000000 1.351157
      
      $pivot
      [1] 1 2
      
      $tol
      [1] 1e-07
      
      $rank
      [1] 2
      
      
      $df.residual
       [1] 4 4 4 4 4 4 4 4 4 4
      
      $sigma
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       2.625995  2.109991  1.823286  1.628639  1.483495  1.369098  1.275549  1.196999 
       protein9 protein10 
       1.129716  1.071179 
      
      $cov.coefficients
                F         M
      F 0.4319541 0.2671159
      M 0.2671159 0.4319541
      
      $stdev.unscaled
                        F         M
      protein1  0.6572322 0.6572322
      protein2  0.6572322 0.6572322
      protein3  0.6572322 0.6572322
      protein4  0.6572322 0.6572322
      protein5  0.6572322 0.6572322
      protein6  0.6572322 0.6572322
      protein7  0.6572322 0.6572322
      protein8  0.6572322 0.6572322
      protein9  0.6572322 0.6572322
      protein10 0.6572322 0.6572322
      
      $pivot
      [1] 1 2
      
      $ndups
      [1] 1
      
      $spacing
      [1] 1
      
      $block
      [1] "control"   "control"   "control"   "treatment" "treatment" "treatment"
      
      $correlation
      [1] 0.5486051
      
      $Amean
       protein1  protein2  protein3  protein4  protein5  protein6  protein7  protein8 
       3.972654  4.189525  4.334591  4.449017  4.545862  4.631041  4.707773  4.778025 
       protein9 protein10 
       4.843098  4.903904 
      
      $method
      [1] "ls"
      
      $design
              F M
      sample1 0 1
      sample2 1 0
      sample3 0 1
      sample4 1 0
      sample5 0 1
      sample6 1 0
      attr(,"assign")
      [1] 1 1
      attr(,"contrasts")
      attr(,"contrasts")$sex
      [1] "contr.treatment"
      
      
      $df.prior
       [1] 38.58884      Inf      Inf      Inf      Inf      Inf      Inf      Inf
       [9]      Inf      Inf
      
      $s2.prior
      [1] 2.880884
      
      $var.prior
      [1] 5.553852 5.553852
      
      $proportion
      [1] 0.01
      
      $s2.post
       [1] 3.257975 2.880884 2.880884 2.880884 2.880884 2.880884 2.880884 2.880884
       [9] 2.880884 2.880884
      
      $t
                       F        M
      protein1  3.688763 3.008819
      protein2  4.015862 3.495404
      protein3  4.092975 3.678374
      protein4  4.161661 3.814839
      protein5  4.224562 3.925568
      protein6  4.283022 4.019822
      protein7  4.337866 4.102550
      protein8  4.389656 4.176712
      protein9  4.438808 4.244228
      protein10 4.485641 4.306411
      
      $df.total
       [1] 40 40 40 40 40 40 40 40 40 40
      
      $p.value
                           F            M
      protein1  6.702158e-04 0.0045222260
      protein2  2.534988e-04 0.0011730460
      protein3  2.007491e-04 0.0006908872
      protein4  1.628886e-04 0.0004623830
      protein5  1.343898e-04 0.0003325038
      protein6  1.123071e-04 0.0002504882
      protein7  9.484042e-05 0.0001949980
      protein8  8.080341e-05 0.0001555746
      protein9  6.937542e-05 0.0001265250
      protein10 5.996918e-05 0.0001045035
      
      $lods
                         F          M
      protein1  -0.4042978 -2.0591941
      protein2   1.5721473 -0.2414474
      protein3   1.8622358  0.3674842
      protein4   2.1252654  0.8418699
      protein5   2.3699841  1.2394903
      protein6   2.6007173  1.5869121
      protein7   2.8200573  1.8986387
      protein8   3.0297491  2.1834902
      protein9   3.2310595  2.4472479
      protein10  3.4249574  2.6939172
      
      $F
       [1]  7.232234  8.892773  9.441926  9.907174 10.319479 10.694430 11.041138
       [8] 11.365461 11.671429 11.961954
      
      $F.p.value
       [1] 1.976640e-03 1.373782e-04 7.932749e-05 4.981602e-05 3.298429e-05
       [6] 2.267086e-05 1.602857e-05 1.158892e-05 8.534199e-06 6.382476e-06
      
      
      $results
      $results$F
                   logFC     CI.L     CI.R average_intensity        t          B
      protein1  4.375960 1.978369 6.773552          3.972654 3.688763 -0.4042978
      protein2  4.479819 2.225246 6.734391          4.189525 4.015862  1.5721473
      protein3  4.565841 2.311269 6.820413          4.334591 4.092975  1.8622358
      protein4  4.642462 2.387889 6.897034          4.449017 4.161661  2.1252654
      protein5  4.712630 2.458058 6.967203          4.545862 4.224562  2.3699841
      protein6  4.777845 2.523272 7.032417          4.631041 4.283022  2.6007173
      protein7  4.839024 2.584452 7.093596          4.707773 4.337866  2.8200573
      protein8  4.896798 2.642225 7.151370          4.778025 4.389656  3.0297491
      protein9  4.951628 2.697056 7.206201          4.843098 4.438808  3.2310595
      protein10 5.003872 2.749299 7.258444          4.903904 4.485641  3.4249574
                     P.Value    adj.P.Val sig.PVal sig.FDR
      protein1  6.702158e-04 0.0006702158        1       1
      protein2  2.534988e-04 0.0002816653        1       1
      protein3  2.007491e-04 0.0002509364        1       1
      protein4  1.628886e-04 0.0002326980        1       1
      protein5  1.343898e-04 0.0002239831        1       1
      protein6  1.123071e-04 0.0002239831        1       1
      protein7  9.484042e-05 0.0002239831        1       1
      protein8  8.080341e-05 0.0002239831        1       1
      protein9  6.937542e-05 0.0002239831        1       1
      protein10 5.996918e-05 0.0002239831        1       1
      
      $results$M
                   logFC     CI.L     CI.R average_intensity        t          B
      protein1  3.569347 1.171756 5.966939          3.972654 3.008819 -2.0591941
      protein2  3.899232 1.644659 6.153804          4.189525 3.495404 -0.2414474
      protein3  4.103340 1.848768 6.357913          4.334591 3.678374  0.3674842
      protein4  4.255571 2.000999 6.510144          4.449017 3.814839  0.8418699
      protein5  4.379093 2.124521 6.633666          4.545862 3.925568  1.2394903
      protein6  4.484237 2.229664 6.738809          4.631041 4.019822  1.5869121
      protein7  4.576522 2.321949 6.831094          4.707773 4.102550  1.8986387
      protein8  4.659252 2.404680 6.913825          4.778025 4.176712  2.1834902
      protein9  4.734568 2.479996 6.989140          4.843098 4.244228  2.4472479
      protein10 4.803936 2.549363 7.058508          4.903904 4.306411  2.6939172
                     P.Value    adj.P.Val sig.PVal sig.FDR
      protein1  0.0045222260 0.0045222260        1       1
      protein2  0.0011730460 0.0013033844        1       1
      protein3  0.0006908872 0.0008636089        1       1
      protein4  0.0004623830 0.0006605472        1       1
      protein5  0.0003325038 0.0005541729        1       1
      protein6  0.0002504882 0.0005009763        1       1
      protein7  0.0001949980 0.0004874949        1       1
      protein8  0.0001555746 0.0004874949        1       1
      protein9  0.0001265250 0.0004874949        1       1
      protein10 0.0001045035 0.0004874949        1       1
      
      
      $tags
      $tags$normalized
      [1] TRUE
      
      $tags$norm_method
      [1] "log2"
      
      $tags$DA_criteria
      $tags$DA_criteria$pval_thresh
      [1] 0.05
      
      $tags$DA_criteria$lfc_thresh
      [1] 1
      
      $tags$DA_criteria$adj_method
      [1] "BH"
      
      
      $tags$extract_intercept
      [1] FALSE
      
      
      attr(,"class")
      [1] "DAList"

