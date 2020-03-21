Sequences = dict()

Sequences["DIF"] = """\chemfig{
               \color{<O1>}{O}% 9
         =[:60]\color{<C14>}{C}% 8
                  (
        -[:120,,,2]\color{<O2>}{O}H% 10
                  )
              -\mcfabove{\color{<C13>}{C}}{\mcfright{H}{_2}}% 7
        -[:300]\color{<C7>}{C}% 5
              -\color{<C8>}{C}% 4
                  (
       =_[:300,,,1]\color{<C9>}{CH}% 3
         -[:240,,1]\mcfbelow{\color{<C10>}{C}}{H}% 2
           =_[:180]\mcfbelow{\color{<C11>}{C}}{H}% 1
        -[:120,,,2]\color{<C12>}{C}H% 6
         =_[:60,,2]\phantom{C}% -> 5
                  )
         -[:60]\mcfabove{\color{<N1>}{N}}{H}% 11
              -\color{<C3>}{C}% 12
       =^[:300]\color{<C4>}{C}% 13
                  (
            -[:240]\color{<Cl4>}{Cl}% 19
                  )
              -\mcfbelow{\color{<C5>}{C}}{H}% 14
    =^[:60,,,1]\color{<C6>}{C}H% 15
     -[:120,,1]\mcfabove{\color{<C1>}{C}}{H}% 16
       =^[:180]\color{<C2>}{C}% 17
                  (
            -[:240]\phantom{C}% -> 12
                  )
        -[:120]\color{<Cl2>}{Cl}% 18
}"""

Sequences["ALA"] = """\chemfig{
                H_3\color{<CB>}{C}% 1
    :[:60,,2,2]H\color{<CA>}{C}% 2
                   (
        -[:120,,2,2]H_2N% 6
                   )
          -[,,2]C% 3
                   (
             =[:300]O% 4
                   )
      -[:60,,,1]OH% 5
}"""